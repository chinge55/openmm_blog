from vina import Vina
from typing import List
from dataclasses import dataclass
import logging 
import io 
import contextlib
logging.basicConfig(level=logging.INFO)

@dataclass
class VinaConfig:
    seed:int = 100
    exhaustiveness:int = 32
    n_poses:int = 20
    box_size = [60, 60, 60]

def seperate_poses(poses:str):
    models = []
    current_model = []
    for line in poses.splitlines():
        if line.startswith("MODEL"):
            if current_model:
                models.append("\n".join(current_model))
                current_model = []
        current_model.append(line)

    if current_model:
        models.append("\n".join(current_model))

    return models

def get_score(curr_model:str):
    for line in curr_model.splitlines():
        if line.startswith("REMARK VINA RESULT"):
            parts = line.split()
            binding_affinity = float(parts[3])
            rmsd_lb = float(parts[4])
            rmsd_ub = float(parts[5])
            return binding_affinity, rmsd_lb, rmsd_ub
    return None, None, None
def conduct_docking(receptor_path:str, ligand_path:str, box_center:List, config, write_poses = False):
    """
    Ref: https://github.com/ccsb-scripps/AutoDock-Vina/issues/112 
    """
    # config = VinaConfig()
    v = Vina(sf_name = 'vina', seed = config.seed)
    v.set_receptor(receptor_path)
    v.set_ligand_from_file(ligand_path)
    v.compute_vina_maps(center=box_center, box_size=config.box_size)
    # energy = v.score()
    # logging.info(f"Score before minimization: {energy[0]} kcal/mol")
    # energy_minimized = v.optimize()
    # logging.info(f"Score after minimization: {energy_minimized[0]} kcal/mol")
    # if write_poses:
    #     pass
    v.dock(exhaustiveness=config.exhaustiveness, n_poses=config.n_poses)
    poses = v.poses()
    models = seperate_poses(poses)
    affinity_arr = []
    rmsd_lb_arr = []
    rmsd_ub_arr = []
    for model in models:
        affinity, rmsd_lb, rmsd_ub = get_score(model)
        affinity_arr.append(affinity)
        rmsd_lb_arr.append(rmsd_lb)
        rmsd_ub_arr.append(rmsd_ub)
    if write_poses:
        out_filename = "out.pdbqt"
        logging.info(f"Writing Autodock Vina Poses on File: {out_filename}")
        v.write_pose(f"temp/{out_filename}")
    return affinity_arr, rmsd_lb_arr, rmsd_ub_arr

if __name__ == "__main__":
    ligand_path ="/home/ubuntu/analysis/docking/src/temp/5330790_pubchem.pdbqt"
    protein_path ="/home/ubuntu/analysis/docking/src/temp/696H.pdbqt"
    affinity_arr, rmsd_lb_arr, rmsd_ub_arr = conduct_docking(protein_path, ligand_path, [114.23, 65.5, 9.55])
    for a,b,c in zip(affinity_arr, rmsd_lb_arr, rmsd_ub_arr):
        logging.info(f"Affinity: {a} RMSD Lower Bound: {b} RMSD Upper Bound: {c}")