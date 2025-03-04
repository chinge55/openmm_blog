from simtk.openmm.app import PDBFile, Modeller, ForceField
from simtk.openmm import unit

# Load prepared protein & ligand
protein_pdb = PDBFile("../results/protein.pdb")  # Your prepared protein
ligand_pdb = PDBFile("../results/ligand.pdb")    # Docked ligand

# Create a Modeller object
modeller = Modeller(protein_pdb.topology, protein_pdb.positions)
modeller.add(ligand_pdb.topology, ligand_pdb.positions)  # Merge protein and ligand

# Save merged complex
with open("../results/protein_ligand_complex.pdb", "w") as f:
    PDBFile.writeFile(modeller.topology, modeller.positions, f)
print("✅ Merged protein and ligand successfully.")
