from core import VinaConfig
import streamlit as st 
import pandas as pd 
import numpy as np 
from molecule_image import get_molecule_image, display_protein
# Displaying Protein Reference: https://discuss.streamlit.io/t/new-component-streamlit-molstar-for-visualisation-and-analysis-of-large-scale-molecular-data/67631 

st.title("Docking Experiment")
pdb_id = st.text_input("Enter the protein ID")
display_protein(pdb_id)
st.write("This is a simple web app to conduct docking experiment")

molecule_name = st.text_input("Enter the molecule name")
st.image(get_molecule_image(molecule_name), caption = f"2D Image of {molecule_name}")