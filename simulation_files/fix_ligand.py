from rdkit import Chem
from rdkit.Chem import AllChem

# Load the original molecule (which has correct bonds)
original_mol = Chem.MolFromMolFile("ligand_original.sdf")

# Load the docked molecule (without sanitization)
docked_mol = Chem.MolFromPDBFile("ligand_docked.pdb", sanitize=False)

if original_mol.GetNumAtoms() != docked_mol.GetNumAtoms():
    raise ValueError("Mismatch in the number of atoms between original and docked ligand.")

conf = docked_mol.GetConformer()
new_conf = Chem.Conformer(original_mol.GetNumAtoms())

for i, atom in enumerate(original_mol.GetAtoms()):
    pos = conf.GetAtomPosition(i)
    new_conf.SetAtomPosition(i, pos)

original_mol.RemoveAllConformers()
original_mol.AddConformer(new_conf)

Chem.MolToMolFile(original_mol, "ligand_final.sdf")

print("Successfully transferred coordinates! Saved as ligand_final.sdf")