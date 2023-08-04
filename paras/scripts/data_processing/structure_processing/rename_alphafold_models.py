import os
from shutil import copy
from sys import argv


def rename_alphafold_models(structure_dir, alignment_dir, output_dir):

    aln_id_to_prot_id = {}

    for file_name in os.listdir(alignment_dir):

        if file_name.endswith('a3m'):
            aln_id = file_name.split('.a3m')[0]
            with open(os.path.join(alignment_dir, file_name)) as alignment:
                prot_id = alignment.readline().strip()[1:]
                aln_id_to_prot_id[aln_id] = prot_id

    for file_name in os.listdir(structure_dir):

        if file_name.endswith('.pdb'):
            file_path = os.path.join(structure_dir, file_name)
            aln_id = file_name.split('_')[0].strip()
            rank = file_name.split('_')[3]
            if rank == '001':
                prot_id = aln_id_to_prot_id[aln_id]
                new_file_path = os.path.join(output_dir, f"{prot_id}.pdb")
                copy(file_path, new_file_path)


if __name__ == "__main__":
    rename_alphafold_models(argv[1], argv[2], argv[3])
