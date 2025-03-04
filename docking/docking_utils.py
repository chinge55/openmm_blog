import core 
import logging
import prody
import matplotlib.pyplot as plt
import streamlit as st
import pandas as pd
import os

def create_vina_config():
    exhaustiveness = st.slider("Exhaustiveness", 1, 32, 8)
    n_poses = st.slider("Number of Poses", 1, 20, 10)
    box_size = st.slider("Box Size", 1, 100, 60)
    config = core.VinaConfig()
    config.exhaustiveness = exhaustiveness
    config.n_poses = n_poses
    config.box_size = [box_size, box_size, box_size]
    return config

def conduct_docking(protein_name, ligand_sdf_path, config, nonstandard_residue, ligand_residue_name, protein_chain, save_poses):
    fname = prody.fetchPDB(protein_name, folder = "temp/", tp = "http")
    protein = prody.parsePDB(fname)
    docking_experiment = core.ProteinPreparation(protein, protein_chain, nonstandard_residue)
    binding_pocket_center = docking_experiment.get_binding_pocket(ligand_residue_name, False)
    condition = docking_experiment.fix_protein()
    if not condition:
        raise Exception("Protein Preparation Failed")
    logging.info("Creating pdbqt for protein")
    protein_pdb_path = f"temp/{docking_experiment.random_path_number}H.pdb"
    protein_pdbqt_path = protein_pdb_path.replace(".pdb", ".pdbqt")
    core.convert_pdb_pdbqt_protein(protein_pdb_path, protein_pdbqt_path)


    ligand_pdbqt_path = ligand_sdf_path.replace(".sdf", ".pdbqt")
    core.convert_sdf_pdbqt_ligand(ligand_sdf_path, ligand_pdbqt_path)
    logging.info(f"Conducting Docking for CID: {ligand_sdf_path}")
    affinity_arr, rmsd_lb_arr, rmsd_ub_arr = core.conduct_docking(protein_pdbqt_path, ligand_pdbqt_path, binding_pocket_center, config, save_poses)
    create_table(affinity_arr, rmsd_lb_arr, rmsd_ub_arr)
    # delete everything in temp folder
    files = os.listdir("temp/")
    for file in files:
        if not file.endswith(".txt"):
            os.remove(f"temp/{file}")



def create_table(affinity_arr, rmsd_lb_arr, rmsd_ub_arr):
    df = pd.DataFrame({'Affinity': affinity_arr, 'RMSD Lower Bound': rmsd_lb_arr, 'RMSD Upper Bound': rmsd_ub_arr})
    st.table(df)