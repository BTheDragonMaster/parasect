def parse_smiles(smiles_file, lowercase=False):
    substrate_to_smiles = {}
    with open(smiles_file, 'r') as names_and_smiles:
        for line in names_and_smiles:
            line = line.strip()
            name, smiles = line.split()
            if lowercase:
                substrate_to_smiles[name.lower()] = smiles
            else:
                substrate_to_smiles[name] = smiles

    return substrate_to_smiles
