from sys import argv
import os


def process_nrpspredictor(nrps_file, out_dir):
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    out_file = os.path.join(out_dir, 'specificities.txt')
    fasta_file = os.path.join(out_dir, 'domains.fasta')

    prot_id_to_nr = {}
    with open(nrps_file, 'r') as nrps:
        with open(fasta_file, 'w') as fasta:
            with open(out_file, 'w') as out:
                out.write("ID\tspec\n")
                nrps.readline()

                for line in nrps:
                    line = line.strip()
                    line_info = line.split('\t')
                    prot_id = line_info[0]
                    sequence = line_info[1]
                    spec = line_info[4]

                    if prot_id not in prot_id_to_nr:
                        prot_id_to_nr[prot_id] = 0

                    prot_id_to_nr[prot_id] += 1
                    identifier = f"{prot_id}.A{prot_id_to_nr[prot_id]}"

                    out.write(f"{identifier}\t{spec}\n")
                    fasta.write(f">{identifier}\n{sequence}\n")


if __name__ == "__main__":
    process_nrpspredictor(argv[1], argv[2])
