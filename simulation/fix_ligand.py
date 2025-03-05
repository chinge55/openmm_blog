from rdkit import Chem 
from rdkit.Chem import AllChem


# original_smiles = "CS(=O)(=O)NC1=C(C(=C(C=C1)F)C(=O)C2=CNC3=C2C=C(C=N3)Cl)F"
template_mol = Chem.MolFromMolFile("/home/sangam/workspace/sangam/openmm_practice/temp/ligand_original.sdf")
pose_mol = Chem.MolFromPDBFile("/home/sangam/workspace/sangam/openmm_practice/results/ligand.pdb", sanitize = False)
# Chem.MolToPDBFile(pose_mol, "pose.pdb")
editable_mol = Chem.EditableMol(pose_mol)
for bond in pose_mol.GetBonds():
    editable_mol.RemoveBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())

pose_mol = editable_mol.GetMol()

corrected_mol = AllChem.AssignBondOrdersFromTemplate(template_mol, pose_mol)
corrected_mol = Chem.AddHs(corrected_mol, addCoords = True)
Chem.MolToPDBFile(corrected_mol, "fixed.pdb")
