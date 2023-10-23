from sys import argv

from paras.scripts.parsers.parsers import parse_domain_list
from paras.scripts.parsers.iterate_over_dir import iterate_over_dir


def remove_data(domains_to_remove, crossval_folder):
    for file_name, domain_file in iterate_over_dir(crossval_folder, extension='.txt'):
        domains = parse_domain_list(domain_file)

        with open(domain_file, 'w') as out:
            for domain in domains:
                if domain not in domains_to_remove:
                    out.write(f"{domain}\n")


if __name__ == "__main__":
    remove_data(['WP_010639243.1.A1', 'EAU39346.1.A1', 'UHJ79951.1.A4', 'OAQ83772.1.A8', 'antaF.A1', 'BBD17766.1.A1', 'AXM43064.1.A1', 'CBW75453.1.A3', 'CBJ79718.1.A1', 'A05163.A1', 'QVQ62857.1.A1', 'ClsI.A1', 'THA23581.1.A1', 'AAC44128.1.A1'], argv[1])