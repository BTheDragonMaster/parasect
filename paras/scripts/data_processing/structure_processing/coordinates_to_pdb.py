import argparse


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=str, required=True, help="Path containing coordinate file")
    parser.add_argument('-o', type=str, required=True, help="Output path.")
    args = parser.parse_args()
    return args


def run():
    args = parse_arguments()
    with open(args.i, 'r') as coordinate_file:
        with open(args.o, 'w') as out:
            counter = 0
            for line in coordinate_file:
                counter += 1
                coordinates = []
                nrs = line.strip().split()
                for nr in nrs:
                    coordinates.append(float(nr))

                pdb_line = f"HETATM{counter:5d} FILT STP E{2:4d}    {coordinates[0]:8.3f}{coordinates[1]:8.3f}{coordinates[2]:8.3f}  0.00  0.00          Ve\n"
                out.write(pdb_line)
            out.write("TER")


if __name__ == "__main__":
    run()
