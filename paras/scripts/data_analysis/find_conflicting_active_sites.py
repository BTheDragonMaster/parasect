import os
from sys import argv
from paras.scripts.parsers.fasta import read_fasta
from paras.scripts.parsers.parsers import parse_specificities


def find_conflicting_active_sites(signature_file, specificity_file, out_folder, strict=1):
    id_to_sig = read_fasta(signature_file)
    id_to_spec = parse_specificities(specificity_file)

    sig_to_ids = {}
    sig_to_spec = {}

    suspicious_sigs = set()

    for seq_id, sig in id_to_sig.items():
        specs = set(id_to_spec[seq_id])
        if sig not in sig_to_ids:
            sig_to_ids[sig] = [seq_id]
            sig_to_spec[sig] = specs

        else:
            sig_to_ids[sig].append(seq_id)
            if strict:
                if specs != sig_to_spec[sig]:
                    suspicious_sigs.add(sig)
            else:
                if not specs.intersection(sig_to_spec[sig]):
                    suspicious_sigs.add(sig)
            sig_to_spec[sig] = sig_to_spec[sig].union(specs)

    write_suspicious_id_groups(suspicious_sigs, sig_to_ids, id_to_spec, out_folder)


def write_suspicious_id_groups(suspicious_sigs, sig_to_ids, id_to_specs, out_folder):
    if not os.path.exists(out_folder):
        os.mkdir(out_folder)

    for sig in suspicious_sigs:
        out_file = os.path.join(out_folder, f"{sig}.txt")
        with open(out_file, 'w') as out:
            out.write(f"domain_name\tspec\n")
            for seq_id in sig_to_ids[sig]:
                spec_string = '|'.join(sorted(id_to_specs[seq_id]))
                out.write(f"{seq_id}\t{spec_string}\n")


if __name__ == "__main__":
    find_conflicting_active_sites(argv[1], argv[2], argv[3], int(argv[4]))




