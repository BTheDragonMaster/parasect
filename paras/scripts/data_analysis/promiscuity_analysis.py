from paras.scripts.parsers.parsers import parse_specificities, parse_domain_list
from paras.scripts.parsers.tabular import Tabular
from sys import argv


def get_promiscuous(domain_list, parasect_data):
    domains = parse_domain_list(domain_list)
    domain_to_spec = parse_specificities(parasect_data)

    promiscuous_domain_to_spec = {}
    for domain in domains:
        if domain in domain_to_spec and len(domain_to_spec[domain]) > 1:
            promiscuous_domain_to_spec[domain] = domain_to_spec[domain]
    return promiscuous_domain_to_spec


def get_promiscuity_type(domain_list, parasect_data, evidence_file, out_file):
    domains = parse_domain_list(domain_list)
    promiscuous_domain_to_spec = get_promiscuous(domain_list, parasect_data)
    evidence_data = Tabular(evidence_file, [2, 4])
    promiscuous_to_evidence = {}
    mildly_promiscuous = 0
    for datapoint in evidence_data.data:
        domain_name = f"{evidence_data.get_value(datapoint, 'protein_name')}.A{evidence_data.get_value(datapoint, 'domain_nr')}"
        for domain in promiscuous_domain_to_spec:
            if domain_name in domain.split('|'):
                evidence_1 = evidence_data.get_value(datapoint, 'evidence_1')
                evidence_2 = evidence_data.get_value(datapoint, 'evidence_2')
                evidence_3 = evidence_data.get_value(datapoint, 'evidence_3')

                if 'atp_ppi_exchange' in [evidence_1, evidence_2, evidence_3]:
                    promiscuous_to_evidence[domain_name] = 'atp_ppi_exchange'
                elif not evidence_1.strip():
                    promiscuous_to_evidence[domain_name] = 'no data'
                elif 'structure_inference' in [evidence_1, evidence_2, evidence_3]:
                    promiscuous_to_evidence[domain_name] = 'structure_inference'
                else:
                    promiscuous_to_evidence[domain_name] = 'other'
                break

    atp_ppi_exchange = list(promiscuous_to_evidence.values()).count('atp_ppi_exchange')
    structure_inference = list(promiscuous_to_evidence.values()).count('structure_inference')
    other = list(promiscuous_to_evidence.values()).count('other')
    no_data = list(promiscuous_to_evidence.values()).count('no data')
    print(promiscuous_to_evidence)
    print("ATP-PPi:", atp_ppi_exchange, "Structure inference:", structure_inference, "Other:", other, "No data:",
          no_data)

    domain_to_spec = parse_specificities(parasect_data)
    for domain in list(domain_to_spec.keys()):
        if domain not in domains:
            del domain_to_spec[domain]

    domain_to_evidence = {}

    mildy_promiscuous_with_atp_ppi = 0
    mildly_promiscuous_no_data = 0

    domain_to_data = {}
    evidence_types = set()
    mildy_promiscuous_structure_inference = 0

    for datapoint in evidence_data.data:
        domain_name = f"{evidence_data.get_value(datapoint, 'protein_name')}.A{evidence_data.get_value(datapoint, 'domain_nr')}"
        for domain in domain_to_spec:
            if domain_name in domain.split('|'):
                evidence_1 = evidence_data.get_value(datapoint, 'evidence_1').strip()
                evidence_2 = evidence_data.get_value(datapoint, 'evidence_2').strip()
                evidence_3 = evidence_data.get_value(datapoint, 'evidence_3').strip()
                if evidence_1:
                    evidence_types.add(evidence_1)
                if evidence_2:
                    evidence_types.add(evidence_2)
                if evidence_3:
                    evidence_types.add(evidence_3)

                minor_spec = evidence_data.get_value(datapoint, "minor specificity").strip()

                if minor_spec:
                    mildly_promiscuous += 1
                    if 'atp_ppi_exchange' in [evidence_1, evidence_2, evidence_3]:
                        mildy_promiscuous_with_atp_ppi += 1
                    elif not evidence_1.strip() and not evidence_2.strip() and not evidence_3.strip():
                        mildly_promiscuous_no_data += 1
                    elif 'structure_inference' in [evidence_1, evidence_2, evidence_3]:
                        mildy_promiscuous_structure_inference += 1
                    else:
                        print(evidence_1, evidence_2, evidence_3)

                if 'atp_ppi_exchange' in [evidence_1, evidence_2, evidence_3]:
                    domain_to_evidence[domain_name] = 'atp_ppi_exchange'
                elif not evidence_1.strip() and not evidence_2.strip() and not evidence_3.strip():
                    domain_to_evidence[domain_name] = 'no data'
                elif 'structure_inference' in [evidence_1, evidence_2, evidence_3]:
                    domain_to_evidence[domain_name] = 'structure_inference'
                else:
                    domain_to_evidence[domain_name] = 'other'
                break

    for datapoint in evidence_data.data:
        domain_name = f"{evidence_data.get_value(datapoint, 'protein_name')}.A{evidence_data.get_value(datapoint, 'domain_nr')}"
        for domain in domain_to_spec:
            if domain_name in domain.split('|'):
                evidence_1 = evidence_data.get_value(datapoint, 'evidence_1').strip()
                evidence_2 = evidence_data.get_value(datapoint, 'evidence_2').strip()
                evidence_3 = evidence_data.get_value(datapoint, 'evidence_3').strip()

                domain_to_data[domain] = {}
                for evidence in evidence_types:
                    domain_to_data[domain][evidence] = 0

                if evidence_1:
                    domain_to_data[domain][evidence_1] = 1
                if evidence_2:
                    domain_to_data[domain][evidence_2] = 1
                if evidence_3:
                    domain_to_data[domain][evidence_3] = 1

                domain_to_data[domain]["Minor specificity"] = 0
                domain_to_data[domain]["Promiscuous"] = 0

                minor_spec = evidence_data.get_value(datapoint, "minor specificity").strip()

                if minor_spec:
                    domain_to_data[domain]["Minor specificity"] = 1

                if domain in promiscuous_domain_to_spec:
                    domain_to_data[domain]["Promiscuous"] = 1

    atp_ppi_exchange = list(domain_to_evidence.values()).count('atp_ppi_exchange')
    structure_inference = list(domain_to_evidence.values()).count('structure_inference')
    other = list(domain_to_evidence.values()).count('other')
    no_data = list(domain_to_evidence.values()).count('no data')
    print("ATP-PPi:", atp_ppi_exchange, "Structure inference:", structure_inference, "Other:", other, "No data:", no_data)
    print(len(domains))
    print(len(promiscuous_domain_to_spec))
    print(mildly_promiscuous)
    print(mildy_promiscuous_with_atp_ppi)
    print(mildly_promiscuous_no_data)
    print(mildy_promiscuous_structure_inference)

    with open(out_file, 'w') as out:
        out.write("Domain\tPromiscuous\tMinor specificity")
        for evidence in sorted(evidence_types):
            out.write(f'\t{evidence}')
        out.write('\n')
        for domain, data in domain_to_data.items():
            out.write(f"{domain}\t{data['Promiscuous']}\t{data['Minor specificity']}")
            for evidence in sorted(evidence_types):
                out.write(f"\t{data[evidence]}")
            out.write('\n')


if __name__ == "__main__":
    get_promiscuity_type(argv[1], argv[2], argv[3], argv[4])


