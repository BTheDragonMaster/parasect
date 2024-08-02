from paras.scripts.parsers.parsers import parse_substrate_list
from paras.scripts.parsers.tabular import Tabular
import paras.data.compound_data

import os
from sys import argv


INCLUDED_SUBSTRATES = [x.lower() for x in parse_substrate_list(os.path.join(os.path.dirname(paras.data.compound_data.__file__), 'all_substrates.txt'))]


def check_mapping_file(mapping_file):
    mapping_data = Tabular(mapping_file, [0])

    for datapoint in mapping_data.data:
        full_name = mapping_data.get_value(datapoint, "full name").lower()
        if full_name not in INCLUDED_SUBSTRATES:
            print(full_name)
        

if __name__ == "__main__":
    check_mapping_file(argv[1])