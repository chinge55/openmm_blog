import streamlit as st
from streamlit_molstar import st_molstar_remote
from .molecule_image import get_molecule_image
from .pubchem_handler import download_sdf
from rdkit import Chem
def get_and_display_protein():
    pdb_id = st.text_input("Enter the protein ID", placeholder = "2b7a")
    pdb_id = pdb_id.upper()
    st_molstar_remote(f"https://files.rcsb.org/view/{pdb_id}.cif", key='sds')
    return pdb_id

def get_and_display_molecule():
    molecule_name = st.text_input("Enter Ligand Pubchem ID", placeholder="164885969")
    cond, path = download_sdf(molecule_name)
    if not cond:
        st.write("Error downloading molecule")
        return
    m = Chem.MolFromMolFile(path)
    st.image(get_molecule_image(m), caption = f"2D Image of {molecule_name}")
    return path

