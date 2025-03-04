import subprocess


def convert_sdf_pdbqt_ligand(sdf_path:str, pdbqt_path:str):
    """
    Openbabel Command:
    $ obabel <ligand_name.sdf> -O <ligand_name.pdbqt>
    """
    command = ["obabel", sdf_path, "-O", pdbqt_path]
    try: 
        result = subprocess.run(command, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error: {e}")
def convert_pdb_pdbqt_protein(pdb_path:str, pdbqt_path:str):
    """
    Openbabel Command:
    $ obabel <protein_name.pdb> -xr -O <protein_name.pdbqt>
    too much work. Just use Subprocess
    """
    command = ["obabel", pdb_path, "-xr", "-O", pdbqt_path]
    try: 
        result = subprocess.run(command, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error: {e}")


if __name__ == "__main__":
    pdb_filepath = "../temp/191H.pdb"
    ligand_filepath = "../temp/5330790_pubchem.sdf"
    convert_sdf_pdbqt_ligand(ligand_filepath, ligand_filepath.replace(".sdf", ".pdbqt"))
    # convert_pdb_pdbqt_protein(pdb_filepath, pdb_filepath.replace(".pdb", ".pdbqt"))