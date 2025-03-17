import os
import core
import prody
import pandas as pd
import logging
import requests
import pubchempy as pcp
logging.basicConfig(level=logging.INFO)
def download_sdf(cid: int, save = True, path = "temp/") -> str:
    """Download the 3D structure of a compound in SDF format from PubChem and return the SDF text."""
    try:
        url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF'
        response = requests.get(url)
        response.raise_for_status()
        sdf =  response.text
    except (requests.exceptions.RequestException, pcp.PubChemHTTPError) as e:
        print(f"Error downloading SDF for CID {cid}: {e}")
        return False, None
    if sdf is None:
        print(f"empty SDF for CID {cid}")
        return False, None 
    #get Smiles of the compound
    download_path = f"{path}{cid}_pubchem.sdf"
    if save:
        with open(download_path, "w") as f:
            f.write(sdf)
    return True, download_path

if __name__ == "__main__":
    pubchem_entries_filepath = "/home/sangam/workspace/sangam/openmm_practice/docking/pubchem_entries.csv"
    temp_path = "temp/"
    os.makedirs(temp_path, exist_ok = True)
    fname = prody.fetchPDB("1P4Q", folder = "temp/")
    pdb = prody.parsePDB(fname)
    docking_experiment = core.ProteinPreparation(pdb, "B", "")
    binding_pocket_center = [-5.0, -6.0, -5.0]
    binding_pocket_center = docking_experiment.get_binding_pocket("IZA", False, binding_pocket_center)
    condition = docking_experiment.fix_protein()
    if not condition:
        raise Exception("Protein Preparation Failed")
    logging.info("Creating pdbqt for protein")
    protein_pdb_path = f"temp/{docking_experiment.random_path_number}H.pdb"
    protein_pdbqt_path = protein_pdb_path.replace(".pdb", ".pdbqt")
    print(f"Protein PDB Path: {protein_pdb_path}")
    core.convert_pdb_pdbqt_protein(protein_pdb_path, protein_pdbqt_path)
    df = pd.read_csv(pubchem_entries_filepath)
    docking_config = core.VinaConfig()
    for pin in range(len(df)):
        try:
            print(f"Working for: {pin}/ {len(df)}")
            cid = df['CID'][pin]
            cid = int(cid)
            logging.info(f"Downloading SDF for CID: {cid}")
            condition, download_path = download_sdf(cid)
            if not condition:
                raise Exception(f"Download Failed for CID: {cid}")
            logging.info(f"Conversion of SDF to PDBQT for CID: {cid}")
            ligand_pdbqt_path = download_path.replace(".sdf", ".pdbqt")
            core.convert_sdf_pdbqt_ligand(download_path, ligand_pdbqt_path)
            logging.info(f"Conducting Docking for CID: {cid}")
            affinity_arr, rmsd_lb_arr, rmsd_ub_arr = core.conduct_docking(protein_pdbqt_path, ligand_pdbqt_path, binding_pocket_center,docking_config, True)
            for a, b, c in zip(affinity_arr, rmsd_lb_arr, rmsd_ub_arr):
                logging.info(f"Affinity: {a}, RMSD Lower Bound: {b}, RMSD Upper Bound: {c}")
            # open a file and save the results in a csv file, only the first index of all required. 
            with open("results.csv", "a") as f:
                f.write(f"{cid},{affinity_arr}\n")
        except Exception as e:
            logging.error(f"Error: {e}")
            continue
