import os
import argparse


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=str, required=True, help="Input directory containing pdb files of A-domains.")
    parser.add_argument('-o', type=str, required=True, help="Output directory.")
    parser.add_argument('-d', type=str, required=True, help="Path to dummy substrate file.")
    args = parser.parse_args()
    return args


def run():
    args = parse_arguments()
    if not os.path.exists(args.o):
        os.mkdir(args.o)

    for pdb_file in os.listdir(args.i):
        pdb_path = os.path.join(args.i, pdb_file)
        if pdb_file.endswith('.pdb') and os.path.isfile(pdb_path):
            out_path = os.path.join(args.o, pdb_file)
            add_dummy_substrate(pdb_path, args.d, out_path)


def get_highest_serial_number(pdb_file):
    highest_serial_number = 0
    with open(pdb_file, 'r') as pdb:
        for line in pdb:
            serial_number = get_atom_serial_number(line)
            if serial_number is not None:
                serial_number = int(serial_number.strip())
                if serial_number > highest_serial_number:
                    highest_serial_number = serial_number
    return highest_serial_number


def add_dummy_substrate(pdb_file, dummy_substrate_file, out_file):
    highest_serial_number = 0
    with open(out_file, 'w') as out:
        with open(pdb_file, 'r') as pdb:
            for line in pdb:
                if line.strip() != "END":
                    out.write(line)
                serial_number = get_atom_serial_number(line)
                if serial_number != None:
                    serial_number = int(serial_number.strip())
                    if serial_number > highest_serial_number:
                        highest_serial_number = serial_number
        highest_serial_number = get_highest_serial_number(pdb_file)

        replacement_serial_numbers = ["{:5d}".format(highest_serial_number + 1),
                                      "{:5d}".format(highest_serial_number + 2),
                                      "{:5d}".format(highest_serial_number + 3),
                                      "{:5d}".format(highest_serial_number + 4)]

        with open(dummy_substrate_file, 'r') as dummy:
            index = 0
            for line in dummy:
                if line.strip() and line.strip() != "TER":
                    new_line = line[:6] + replacement_serial_numbers[index] + line[11:]
                    out.write(new_line)
                elif line.strip() == "TER":
                    out.write("TER\n")
            out.write("END")


def get_atom_serial_number(line):
    if line.startswith("ATOM") or line.startswith("HETATM"):
        atom_serial_number = line[6:11]
        return atom_serial_number

    return None


if __name__ == "__main__":
    run()
