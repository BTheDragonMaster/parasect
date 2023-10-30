import os
from argparse import ArgumentParser

from pymol import cmd


from paras.scripts.data_processing.structure_processing.add_dummy_substrate import add_dummy_substrate
from paras.scripts.parsers.parsers import parse_sequence_feature_importances
from paras.scripts.parsers.pdb import sequence_from_pdb
from paras.scripts.machine_learning.feature_inference.plot_importances_pymol import importance_to_colour
from paras.scripts.feature_extraction.sequence_feature_extraction.hmm.run_hmmscan import run_hmmpfam2_single, HMM2_FILE
from paras.scripts.feature_extraction.sequence_feature_extraction.extract_domains import parse_hmm2_results
from paras.scripts.feature_extraction.sequence_feature_extraction.sequence_labels import parse_domain_id
from paras.scripts.feature_extraction.sequence_feature_extraction.read_positions import \
    HMM2_POSITIONS_EXTENDED_SIGNATURE


def parse_arguments():
    parser = ArgumentParser(description="Visualise sequence feature importance in pymol.")
    parser.add_argument('-i', required=True, help="Feature importance file.")
    parser.add_argument('-o', required=True, help="Output directory.")
    parser.add_argument('-pdb', type=str, required=True, help="PDB file to project features on.")
    parser.add_argument('-d', type=str, default=None, help="Link to file with PDB dummy substrate")
    parser.add_argument('-header', action='store_true', help="If feature importance files contains header, give this argument")
    args = parser.parse_args()

    return args


def parse_importance_label(label):
    feature_name, _, feature_residue = label.split('_')
    return feature_name, int(feature_residue)


def residues_from_pdb(pdb_file):
    sequence, residue_numbers = sequence_from_pdb(pdb_file)
    hmm_output = run_hmmpfam2_single(HMM2_FILE, f">query\n{sequence}")
    positions = None
    for hit_key, hit in parse_hmm2_results(hmm_output).items():
        seq_id, hit_id, hit_start, hit_end = parse_domain_id(hit_key)
        if seq_id == 'query' and hit_id == 'AMP-binding':
            hmm_signature_positions = HMM2_POSITIONS_EXTENDED_SIGNATURE

            profile = hit.aln[1].seq
            query = hit.aln[0].seq
            offset = hit.hit_start
            query_offset = hit.query_start

            hmm_positions = [p - offset for p in hmm_signature_positions]

            positions = []
            query_position_skipping_gaps = 0
            hmm_position_skipping_gaps = 0

            for i, amino in enumerate(query):
                hmm_amino = profile[i]
                if hmm_amino in '.-':
                    if amino != '-':
                        query_position_skipping_gaps += 1
                    continue
                if hmm_position_skipping_gaps in hmm_positions:
                    positions.append(query_position_skipping_gaps + query_offset)

                hmm_position_skipping_gaps += 1
                if amino != '-':
                    query_position_skipping_gaps += 1

            print(positions)
            print("".join([sequence[i] for i in positions]))

    return positions


def get_importance_per_residue(feature_to_importance):
    residue_to_importance = {}
    for feature, importance in feature_to_importance.items():
        name, residue = parse_importance_label(feature)
        if residue not in residue_to_importance:
            residue_to_importance[residue] = 0
        residue_to_importance[residue] += importance

    return residue_to_importance


def sort_by_feature_type(feature_to_importance):
    feature_type_to_feature_to_importance = {}
    for feature, importance in feature_to_importance.items():
        feature_type, residue = parse_importance_label(feature)
        if feature_type not in feature_type_to_feature_to_importance:
            feature_type_to_feature_to_importance[feature_type] = {}

        feature_type_to_feature_to_importance[feature_type][residue] = importance
    return feature_type_to_feature_to_importance


def normalise_importances(residue_to_importance):
    max_importance = 0.0

    for residue, importance in residue_to_importance.items():
        if abs(importance) > max_importance:
            max_importance = abs(importance)

    scaling_factor = 1.0 / max_importance

    for residue, importance in residue_to_importance.items():
        residue_to_importance[residue] = importance * scaling_factor


def run():
    args = parse_arguments()

    if not os.path.exists(args.o):
        os.mkdir(args.o)

    signature_positions = residues_from_pdb(args.pdb)
    signature_residue_nrs = [nr + 1 for nr in signature_positions]
    feature_to_importance = parse_sequence_feature_importances(args.i, args.header)
    residue_to_importance = get_importance_per_residue(feature_to_importance)

    normalise_importances(residue_to_importance)
    pdb_dummy = os.path.join(args.o, "domain_with_substrate.pdb")
    if args.d:
        add_dummy_substrate(args.pdb, args.d, pdb_dummy)
        pdb_in = pdb_dummy
    else:
        pdb_in = args.pdb

    residue_importance_file = os.path.join(args.o, 'residue_importances.pse')

    visualise_in_pymol(pdb_in, signature_residue_nrs, residue_to_importance, residue_importance_file)

    feature_type_to_importances = sort_by_feature_type(feature_to_importance)

    for feature_type, residue_to_importance in feature_type_to_importances.items():
        residue_importance_file = os.path.join(args.o, f'{feature_type}_importances.pse')
        normalise_importances(residue_to_importance)

        visualise_in_pymol(pdb_in, signature_residue_nrs, residue_to_importance, residue_importance_file)


def visualise_in_pymol(pdb_in, positions, residue_to_importance, pse_out):
    cmd.load(pdb_in, object="a_domain")
    cmd.select('chA', 'chain A')
    cmd.select('chB', 'chain B')
    cmd.color('grey40', 'chA')
    cmd.color('grey40', 'chB')
    cmd.delete('chA')
    cmd.delete('chB')

    for i, position in enumerate(positions):
        importance = residue_to_importance[i + 1]
        colour = importance_to_colour(importance)
        selection_name = f"residue_{position}"

        cmd.select(selection_name, f"resi {position}")
        cmd.set_color(selection_name, list(colour))
        cmd.color(selection_name, selection_name)
        cmd.set("sphere_transparency", f"{1 - abs(importance)}", selection_name)
        cmd.show("sticks", selection_name)
        cmd.hide("cartoon", selection_name)
        # cmd.delete(selection_name)

    cmd.save(pse_out)
    cmd.reinitialize()


if __name__ == "__main__":
    run()