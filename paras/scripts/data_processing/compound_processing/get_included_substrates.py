from paras.scripts.data_analysis.count_substrates import count_substrates_from_file
from sys import argv


def get_included_substrates(specificities_file, domain_list, limit, out_file):
    spec_to_count = count_substrates_from_file(domain_list, specificities_file)
    with open(out_file, 'w') as out:
        for spec, count in spec_to_count.items():
            if count >= limit:
                out.write(f"{spec}\n")


if __name__ == "__main__":
    get_included_substrates(argv[1], argv[2], int(argv[3]), argv[4])

