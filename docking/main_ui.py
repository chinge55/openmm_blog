import streamlit as st 
import ui_utils
import core
from docking_utils import conduct_docking, create_vina_config

st.title("Docking Experiment")
pdb_id = ui_utils.get_and_display_protein()
protein_chain = st.text_input("Enter Protein Chain", placeholder = "A")
molecule_path = ui_utils.get_and_display_molecule()
config = create_vina_config()
save_poses = True

st.title("Docking Experiment Options")
nonstandard_residue = st.text_input("Enter Nonstandard Residue",placeholder="PTR", help = "There can be nonstandard residue, please check rcsb website")
ligand_residue_name = st.text_input("Enter Ligand Residue Name",placeholder = "IZA", help = "Enter the residue name of the ligand")
if st.button("Conduct Docking"):
    with st.spinner("Conducting Docking"):
        conduct_docking(pdb_id, molecule_path, config, nonstandard_residue, ligand_residue_name, protein_chain, save_poses)

    with open("temp/out_pdbqt.txt") as file: 
        st.download_button(
            label = "Download PDBQT File",
            file_name = "out.pdbqt",
            data = file
        )
