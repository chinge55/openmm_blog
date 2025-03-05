import numpy as np
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
import mdtraj as md
import pdbfixer
import openmm as mm
import openmm.app as app
from openmm import unit
from openff.toolkit.topology import Molecule, Topology
from openmmforcefields.generators import GAFFTemplateGenerator
from rdkit import Chem
from rdkit.Chem import Draw
import time
import sys

# --- Input Files ---
ligand_file = "/home/sangam/workspace/sangam/openmm_practice/simulation_files/ligand_final.sdf"     # Molecule file for ligand (SDF or MOL2 recommended)
# ligand_file = "/home/sangam/workspace/sangam/openmm_practice/temp/ligand_original.sdf"     # Molecule file for ligand (SDF or MOL2 recommended)
protein_pdb_file = "/home/sangam/workspace/sangam/openmm_practice/simulation_files/prepared_protein.pdb"     # Molecule file for protein (PDB format)

def sanitize_ligand(original_sdf_path, docked_pdb_path):
    # Load the original molecule (which has correct bonds)
    original_mol = Chem.MolFromMolFile(original_sdf_path)

    # Load the docked molecule (without sanitization)
    docked_mol = Chem.MolFromPDBFile(docked_pdb_path, sanitize=False)

    if original_mol.GetNumAtoms() != docked_mol.GetNumAtoms():
        raise ValueError("Mismatch in the number of atoms between original and docked ligand.")

    conf = docked_mol.GetConformer()
    new_conf = Chem.Conformer(original_mol.GetNumAtoms())

    for i, atom in enumerate(original_mol.GetAtoms()):
        pos = conf.GetAtomPosition(i)
        new_conf.SetAtomPosition(i, pos)

    original_mol.RemoveAllConformers()
    original_mol.AddConformer(new_conf)
    return original_mol

def prepare_ligand(ligand_file):
    """
    Load the ligand molecule and add hydrogens.
    """
    # split molecule
    rdkit_mol = Chem.MolFromMolFile(ligand_file)
    if rdkit_mol is None:
        raise ValueError(f"Could not read molecule from {pdb_file}")
    print(f"Loaded Molecule")
    prepared_ligand = Chem.rdmolops.AddHs(rdkit_mol, addCoords = True)
    prepared_ligand = Chem.MolFromMolBlock(Chem.MolToMolBlock(prepared_ligand))
    return rdkit_mol 

def save_mol_as_image(mol, filename="molecule.png", size=(300, 300)):
    """Save an RDKit molecule as an image."""
    img = Draw.MolToImage(mol, size=size)
    img.save(filename)
    print(f"Saved molecule image as {filename}")

def rdkit_to_openmm(rdkit_mol, name="LIG"):
    """
    Convert an RDKit molecule to an OpenMM molecule.
    Inspired by @hannahbrucemcdonald and @glass-w.

    Parameters
    ----------
    rdkit_mol: rdkit.Chem.rdchem.Mol
        RDKit molecule to convert.
    name: str
        Molecule name.

    Returns
    -------
    omm_molecule: openmm.app.Modeller
        OpenMM modeller object holding the molecule of interest.
    """
    # convert RDKit to OpenFF
    off_mol = Molecule.from_rdkit(rdkit_mol)

    # add name for molecule
    off_mol.name = name

    # add names for atoms
    element_counter_dict = {}
    for off_atom, rdkit_atom in zip(off_mol.atoms, rdkit_mol.GetAtoms()):
        element = rdkit_atom.GetSymbol()
        if element in element_counter_dict.keys():
            element_counter_dict[element] += 1
        else:
            element_counter_dict[element] = 1
        off_atom.name = element + str(element_counter_dict[element])

    # convert from OpenFF to OpenMM
    off_mol_topology = off_mol.to_topology()
    mol_topology = off_mol_topology.to_openmm()
    mol_positions = off_mol.conformers[0]

    # convert units from Ångström to nanometers
    # since OpenMM works in nm
    mol_positions = mol_positions.to("nanometers")

    # combine topology and positions in modeller object
    omm_mol = app.Modeller(mol_topology, mol_positions)

    return omm_mol

def prepare_protein(
    pdb_file, ignore_missing_residues=True, ignore_terminal_missing_residues=True, ph=7.0
):
    """
    Use pdbfixer to prepare the protein from a PDB file. Hetero atoms such as ligands are
    removed and non-standard residues replaced. Missing atoms to existing residues are added.
    Missing residues are ignored by default, but can be included.

    Parameters
    ----------
    pdb_file: pathlib.Path or str
        PDB file containing the system to simulate.
    ignore_missing_residues: bool, optional
        If missing residues should be ignored or built.
    ignore_terminal_missing_residues: bool, optional
        If missing residues at the beginning and the end of a chain should be ignored or built.
    ph: float, optional
        pH value used to determine protonation state of residues

    Returns
    -------
    fixer: pdbfixer.pdbfixer.PDBFixer
        Prepared protein system.
    """
    fixer = pdbfixer.PDBFixer(str(pdb_file))
    fixer.removeHeterogens()  # co-crystallized ligands are unknown to PDBFixer
    fixer.findMissingResidues()  # identify missing residues, needed for identification of missing atoms

    # if missing terminal residues shall be ignored, remove them from the dictionary
    if ignore_terminal_missing_residues:
        chains = list(fixer.topology.chains())
        keys = fixer.missingResidues.keys()
        for key in list(keys):
            chain = chains[key[0]]
            if key[1] == 0 or key[1] == len(list(chain.residues())):
                del fixer.missingResidues[key]

    # if all missing residues shall be ignored ignored, clear the dictionary
    if ignore_missing_residues:
        fixer.missingResidues = {}

    fixer.findNonstandardResidues()  # find non-standard residue
    fixer.replaceNonstandardResidues()  # replace non-standard residues with standard one
    fixer.findMissingAtoms()  # find missing heavy atoms
    fixer.addMissingAtoms()  # add missing atoms and residues
    fixer.addMissingHydrogens(ph)  # add missing hydrogens
    return fixer

def merge_protein_and_ligand(protein, ligand):
    """
    Merge two OpenMM objects.

    Parameters
    ----------
    protein: pdbfixer.pdbfixer.PDBFixer
        Protein to merge.
    ligand: openmm.app.Modeller
        Ligand to merge.

    Returns
    -------
    complex_topology: openmm.app.topology.Topology
        The merged topology.
    complex_positions: openmm.unit.quantity.Quantity
        The merged positions.
    """
    # combine topologies
    md_protein_topology = md.Topology.from_openmm(protein.topology)  # using mdtraj for protein top
    md_ligand_topology = md.Topology.from_openmm(ligand.topology)  # using mdtraj for ligand top
    md_complex_topology = md_protein_topology.join(md_ligand_topology)  # add them together
    complex_topology = md_complex_topology.to_openmm()

    # combine positions
    total_atoms = len(protein.positions) + len(ligand.positions)

    # create an array for storing all atom positions as tupels containing a value and a unit
    # called OpenMM Quantities
    complex_positions = unit.Quantity(np.zeros([total_atoms, 3]), unit=unit.nanometers)
    complex_positions[: len(protein.positions)] = protein.positions  # add protein positions
    complex_positions[len(protein.positions) :] = ligand.positions  # add ligand positions

    return complex_topology, complex_positions

def generate_forcefield(
    rdkit_mol=None, protein_ff="amber14-all.xml", solvent_ff="amber14/tip3pfb.xml"
):
    """
    Generate an OpenMM Forcefield object and register a small molecule.

    Parameters
    ----------
    rdkit_mol: rdkit.Chem.rdchem.Mol
        Small molecule to register in the force field.
    protein_ff: string
        Name of the force field.
    solvent_ff: string
        Name of the solvent force field.

    Returns
    -------
    forcefield: openmm.app.Forcefield
        Forcefield with registered small molecule.
    """
    forcefield = app.ForceField(protein_ff, solvent_ff)

    if rdkit_mol is not None:
        gaff = GAFFTemplateGenerator(
            molecules=Molecule.from_rdkit(rdkit_mol, allow_undefined_stereo=True)
        )
        forcefield.registerTemplateGenerator(gaff.generator)

    return forcefield

if __name__ == "__main__":
    ligand_name = "UNL"
    prepared_ligand =prepare_ligand(ligand_file)
    # omm_ligand = rdkit_to_openmm(prepared_ligand, ligand_name)
    # prepared_protein = prepare_protein(protein_pdb_file)
    # complex_topology, complex_positions = merge_protein_and_ligand(prepared_protein, omm_ligand)
    # print("Complex topology has", complex_topology.getNumAtoms(), "atoms.")
    # from openmm.app import PDBFile
    # with open("complex.pdb", "w") as pdb_file:
    #     PDBFile.writeFile(complex_topology, complex_positions, pdb_file)
    # forcefield = generate_forcefield(prepared_ligand)
    # print(f'Forcefield Generated')
    # modeller = app.Modeller(complex_topology, complex_positions)
    # start = time.time()
    # modeller.addSolvent(forcefield, padding=1 * unit.nanometers, ionicStrength=0.15 * unit.molar)
    # end = time.time()
    # print(f"Added solvent in {end-start} seconds")
    # print(f"Added solvent to the system. Creating OpenMM system...")
    # system = forcefield.createSystem(modeller.topology, nonbondedMethod=app.PME)
    # integrator = mm.LangevinIntegrator(
    #     300 * unit.kelvin, 1.0 / unit.picoseconds, 2.0 * unit.femtoseconds
    # )
    # simulation = app.Simulation(modeller.topology, system, integrator)
    # simulation.context.setPositions(modeller.positions)
    # print(f"Minimizing energy...")
    # start_time = time.time()
    # simulation.minimizeEnergy()
    # end_time = time.time()
    # print(f"Minimized energy in {end_time-start_time} seconds")
    # with open("topology.pdb", "w") as pdb_file:
    #     app.PDBFile.writeFile(
    #         simulation.topology,
    #         simulation.context.getState(getPositions=True, enforcePeriodicBox=True).getPositions(),
    #         file=pdb_file,
    #         keepIds=True,
    #     )
    # steps = 50000  # corresponds to 20 fs
    # write_interval = 5000  # write every 2 fs
    # log_interval = 2500  # log progress to stdout every 2 fs
    # simulation.reporters.append(
    # md.reporters.XTCReporter(file=str("trajectory.xtc"), reportInterval=write_interval)
    # )
    # simulation.reporters.append(
    #     app.StateDataReporter(
    #         sys.stdout,
    #         log_interval,
    #         step=True,
    #         potentialEnergy=True,
    #         temperature=True,
    #         progress=True,
    #         remainingTime=True,
    #         speed=True,
    #         totalSteps=steps,
    #         separator="\t",
    #     )
    # )
    # simulation.context.setVelocitiesToTemperature(300 * unit.kelvin)
    # simulation.step(steps)  # perform the simulation