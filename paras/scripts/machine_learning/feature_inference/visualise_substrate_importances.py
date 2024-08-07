import os
from argparse import ArgumentParser

from pymol import cmd


from paras.scripts.data_processing.structure_processing.add_dummy_substrate import add_dummy_substrate
from paras.scripts.parsers.parsers import parse_residue_feature_importance_matrix

from paras.scripts.machine_learning.feature_inference.plot_importances_pymol import importance_to_colour

from paras.scripts.machine_learning.feature_inference.importance_per_sequence_property import residues_from_pdb


def parse_arguments():
    parser = ArgumentParser(description="Visualise sequence feature importance in pymol.")
    parser.add_argument('-i', required=True, help="Feature importance file.")
    parser.add_argument('-o', required=True, help="Output directory.")
    parser.add_argument('-pdb', type=str, required=True, help="Folder containing representative .pdb files for each substrate.")
    parser.add_argument('-d', type=str, default=None, help="Link to file with PDB dummy substrate")
    args = parser.parse_args()

    return args

SUBSTRATE_TO_PYMOL = {"alanine": "ALA",
                      "cysteine": "CYS",
                      "aspartic acid": "ASP",
                      "glutamic acid": "GLU",
                      "phenylalanine": "PHE",
                      "glycine": "GLY",
                      "histidine": "HIS",
                      "isoleucine": "ILE",
                      "lysine": "LYS",
                      "leucine": "LEU",
                      "methionine": "MET",
                      "asparagine": "ASN",
                      "proline": "PRO",
                      "glutamine": "GLN",
                      "arginine": "ARG",
                      "serine": "SER",
                      "threonine": "THR",
                      "valine": "VAL",
                      "tryptophan": "TRP",
                      "tyrosine": "TYR"}
def run():
    args = parse_arguments()

    if not os.path.exists(args.o):
        os.mkdir(args.o)

    substrate_to_importances = parse_residue_feature_importance_matrix(args.i)

    for substrate, importances in substrate_to_importances.items():
        print(substrate)
        pdb_in = os.path.join(args.pdb, f"{substrate}.pdb")

        signature_positions = residues_from_pdb(pdb_in)
        signature_residue_nrs = [nr + 1 for nr in signature_positions]
        position_to_importance = {}
        for i, importance in enumerate(importances):
            position = signature_residue_nrs[i]
            position_to_importance[position] = importance

        pdb_dummy = os.path.join(args.o, f"domain_with_{substrate}.pdb")

        if args.d:
            add_dummy_substrate(pdb_in, args.d, pdb_dummy)
            pdb_in = pdb_dummy

        residue_importance_file = os.path.join(args.o, f'residue_importances_{substrate}.pse')
        visualise_in_pymol(pdb_in, substrate, position_to_importance, residue_importance_file)


def visualise_in_pymol(pdb_in, substrate, position_to_importance, pse_out):
    cmd.load(pdb_in, object="a_domain")
    cmd.select('chA', 'chain A')
    cmd.select('chB', 'chain B')
    cmd.color('grey40', 'chA')
    cmd.color('grey40', 'chB')
    cmd.delete('chA')
    cmd.delete('chB')

    for position, importance in position_to_importance.items():
        colour = importance_to_colour(importance)
        selection_name = f"residue_{position}"

        cmd.select(selection_name, f"resi {position}")
        cmd.set_color(selection_name, list(colour))
        cmd.color(selection_name, selection_name)
        cmd.set("sphere_transparency", f"{1 - abs(importance)}", selection_name)
        cmd.show("sticks", selection_name)
        cmd.hide("cartoon", selection_name)

        # cmd.delete(selection_name)

    cmd.select('oxygens', 'name oe1+od1+od2+oe2+og1+og+oh')
    cmd.color('firebrick', 'oxygens')
    cmd.delete('oxygens')

    cmd.select('nitrogens', 'name nd2+ne2+nz+nh1+nh2+ne+nd1+ne2+ne1')
    cmd.color('deepblue', 'nitrogens')
    cmd.delete('nitrogens')

    cmd.select('sulphurs', 'name sd+sg')
    cmd.color('brightorange', 'sulphurs')
    cmd.delete('sulphurs')

    if substrate in SUBSTRATE_TO_PYMOL:

        cmd.wizard("mutagenesis")
        cmd.do("refresh_wizard")

        # Mutate
        cmd.get_wizard().set_dep('dep')
        cmd.get_wizard().set_mode(SUBSTRATE_TO_PYMOL[substrate])
        print(cmd.get_wizard().dep)
        cmd.get_wizard().do_select("/a_domain//C/1")

        print("Selected")

        # cmd.frame("1")
        cmd.get_wizard().apply()
        cmd.set_wizard("done")

    cmd.save(pse_out)
    cmd.reinitialize()


if __name__ == "__main__":
    run()