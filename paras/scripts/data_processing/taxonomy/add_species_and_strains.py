from parasect.core.parsing import parse_taxonomy_file
from parasect.core.tabular import Tabular
from sys import argv

def parse_species(species_file):
    syn_to_species = {}
    with open(species_file, 'r') as species_info:
        species_info.readline()
        for line in species_info:
            line = line.strip()
            if line:
                tax_info = line.split('\t')
                if len(tax_info) == 2:
                    syn_to_species[tax_info[0]] = (tax_info[1], "")
                else:
                    syn_to_species[tax_info[0]] = (tax_info[1], tax_info[2])


    return syn_to_species

taxonomy = parse_taxonomy_file(argv[1])
syn_to_species = parse_species(argv[2])
out_file = argv[3]

protein_to_species = {}
issues = False
for protein, tax in taxonomy.items():
    protein_synonyms = protein.split('|')
    match_found = False
    for synonym in protein_synonyms:
        if synonym in syn_to_species:
            species, strain = syn_to_species[synonym]
            species = species.strip()
            strain = strain.strip()
            if strain and species.endswith(strain) and len(species.split()) > 2:
                print(species, strain)
                species = species[:-len(strain)].strip()
                print(species)

            if not strain:
                strain = "Unknown"

            if species.split()[0] != tax.genus:
                print(f"Mismatching genera for {protein}: {species.split()[0]}, {tax.genus}")

            if not match_found:
                protein_to_species[protein] = (species, strain)
            else:
                if protein_to_species[protein] != (species, strain):
                    print(f"Mismatching species for protein {protein}: ")
                    print(species, strain)
                    issues = True
            match_found = True


    if not match_found:
        print(f"No species found for {protein}")

if issues:
    raise ValueError

with open(out_file, 'w') as out:
    out.write("protein_id\tdomain\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\tstrain\n")
    for protein, tax in taxonomy.items():
        species, strain = protein_to_species[protein]
        out.write(f"{protein}\t{tax.domain}\t{tax.kingdom}\t{tax.phylum}\t{tax.cls}\t{tax.order}\t{tax.family}\t{tax.genus}\t{species}\t{strain}\n")
