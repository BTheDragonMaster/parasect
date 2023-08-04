from scipy.spatial import distance
import subprocess
from subprocess import CalledProcessError
import os
import shutil
import argparse


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', type=float, default=3.4, help="Minimum sphere size.")
    parser.add_argument('-c', type=float, default=3.0, help="Cluster distance.")
    parser.add_argument('-d', type=float, default=4.0, help="Minimum distance between pocket and substrate.")
    parser.add_argument('-m', type=float, default=6.2, help="Maximum sphere size")
    parser.add_argument('-i', type=str, required=True, help="Input directory containing .pdb files")
    parser.add_argument('-o', type=str, required=True, help="Output directory")
    args = parser.parse_args()
    return args


def run():

    args = parse_arguments()

    if not os.path.exists(args.o):
        os.mkdir(args.o)
    counter = 0
    for pdb_file in os.listdir(args.i):
        if pdb_file.endswith('.pdb'):
            counter += 1
            pdb_name = pdb_file.split('.pdb')[0]
            try:
                run_fpocket(args.i, pdb_file, args.s, args.m, args.c, args.o)
                pdb_out = os.path.join(args.o, f"{pdb_name}_out/{pdb_name}_out.pdb")
                get_nearby_pockets(pdb_out, pdb_name, os.path.join(args.o, f"{pdb_name}_out/"), args.d)

                if counter % 100 == 0:
                    print(f"Processed {counter} structures.")
            except CalledProcessError:
                print(f"Failed to run fPocket on {pdb_name.split('_0')[0]}.")


def get_nearby_pockets(pdb_file, pdb_name, out_dir, min_distance_to_substrate):
    min_dist = 10000
    beta_carbon_pos = (32.833,  97.884,  32.711)
    alpha_carbon_pos = (32.341,  99.049,  31.853)
    pocket_to_coords = parse_pocket_from_pdb(pdb_file)
    nearby_pockets = []
    pocket_coords = []
    for pocket, coords in pocket_to_coords.items():
        for coord in coords:
            dist_1 = distance.euclidean(alpha_carbon_pos, coord)
            dist_2 = distance.euclidean(beta_carbon_pos, coord)
            if dist_1 < min_distance_to_substrate or dist_2 < min_distance_to_substrate:
                nearby_pockets.append(pocket)
                pocket_coords.append(coords)
                break

            if dist_1 < min_dist:
                min_dist = dist_1
            if dist_2 < min_dist:
                min_dist = dist_2

    summary_file = os.path.join(out_dir, f'substrate_pockets_{pdb_name.split("_0")[0]}.txt')
    with open(summary_file, 'w') as pockets:
        for i, pocket in enumerate(nearby_pockets):
            for x, y, z in pocket_coords[i]:
                pockets.write(f"{pocket}\t{x}\t{y}\t{z}\n")

    if len(nearby_pockets) == 0:
        print(f"No pockets found for {pdb_name.split('_0')[0]}. Closest pocket: {min_dist}")


def run_fpocket(in_dir, pdb_file, sphere_size, max_sphere, cluster_distance, out_dir):
    in_file = os.path.join(in_dir, pdb_file)
    pdb_name = pdb_file.split('.pdb')[0]
    tmp_dir = os.path.join(out_dir, "tmp_fpocket.tmp")
    with open(tmp_dir, 'w') as tmp:
        subprocess.check_call(["fpocket", '-f', in_file, "-m", str(sphere_size), "-D", str(cluster_distance), "-M",
                               str(max_sphere)], stdout=tmp)

    for pdb in os.listdir(in_dir):
        full_path = os.path.join(in_dir, pdb)
        if pdb == f"{pdb_name}_out" and os.path.isdir(full_path):
            destination = os.path.join(out_dir, pdb)
            if os.path.exists(destination):
                shutil.rmtree(destination)

            shutil.move(full_path, destination)
    os.remove(tmp_dir)


def parse_pocket_from_pdb(pdb_file):
    pocket_to_coords = {}
    with open(pdb_file, 'r') as pdb:
        for line in pdb:
            if line.startswith("HETATM"):
                atom_name = line[12:16].strip()

                if atom_name in {"POL", "APOL"}:
                    pocket_nr = int(line[22:26].strip())
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    if pocket_nr not in pocket_to_coords:
                        pocket_to_coords[pocket_nr] = []
                    pocket_to_coords[pocket_nr].append((x, y, z))
    return pocket_to_coords


if __name__ == "__main__":
    run()
