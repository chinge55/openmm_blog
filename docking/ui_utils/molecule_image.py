import rdkit 
from rdkit import Chem
from rdkit.Chem import Draw
import numpy as np
from streamlit_molstar import st_molstar, st_molstar_rcsb, st_molstar_remote



def get_molecule_image(molecule:str):
    """
    Returns 2D image of the molecule as numpy array

    Params:
        molecule: str: SMILES representation of the molecule
    
    Returns:
        molecule_img: np.array: 2D image of the molecule
    """
    molecule_img = Draw.MolToImage(molecule)
    molecule_img_np = np.array(molecule_img)
    return molecule_img_np


if __name__ == "__main__":
    get_molecule_image("CCO")