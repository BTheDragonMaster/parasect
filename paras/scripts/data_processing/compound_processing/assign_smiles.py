from paras.scripts.parsers.smiles import parse_smiles
from sys import argv


def get_unique_substrate_names(domain_to_substrate):
    pass


def get_unique_substrates(substrate_file):
    unique_substrates = set()
    with open(substrate_file, 'r') as substrate_names:
        for line in substrate_names:
            line = line.strip()
            if line:
                substrates = line.split('|')
                for substrate in substrates:
                    unique_substrates.add(substrate)

    return unique_substrates


def assign_smiles(substrate_file, smiles_file, out_file):
    unique_substrates = get_unique_substrates(substrate_file)
    substrate_to_smiles = parse_smiles(smiles_file, lowercase=True)
    with open(out_file, 'w') as out:
        out.write("substrate\tsmiles\n")
        for substrate in unique_substrates:
            if substrate.lower() in substrate_to_smiles:
                smiles = substrate_to_smiles[substrate.lower()]
            else:
                smiles = ''

            out.write(f"{substrate}\t{smiles}\n")


if __name__ == "__main__":
    assign_smiles(argv[1], argv[2], argv[3])
