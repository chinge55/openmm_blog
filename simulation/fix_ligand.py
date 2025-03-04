from rdkit import Chem
from rdkit.Chem import AllChem

# Load the ligand
ligand = Chem.SDMolSupplier("../results/ligand.sdf")[0]

# Assign stereochemistry
Chem.AssignStereochemistry(ligand, force=True, cleanIt=True)

# Save the updated molecule
w = Chem.SDWriter("../results/ligand_fixed.sdf")
w.write(ligand)
w.close()

print("✅ Stereochemistry fixed and saved as 'ligand_fixed.sdf'")
