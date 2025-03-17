import pdbfixer 
import prody
import random
import logging 
import openmm 
from openmm import app, unit
logging.basicConfig(level=logging.INFO)

class ProteinPreparation:
    allowed_chains = ["A","B"]
    def __init__(self, receptor, chain:str, nonstandard_residue:str ):
        self.receptor = receptor
        if chain not in self.allowed_chains:
            raise ValueError(f"Chain {chain} not allowed. Choose from {self.allowed_chains}")
        self.chain = chain
        if nonstandard_residue == "":
            nonstandard_residue = None
        self.nonstandard_residue = nonstandard_residue
        self.get_chain()
    
    def get_chain(self):
        """Get the protein chain from the receptor. If nonstandard_residue is provided, remove it from the protein
        Raises: 
            Error if chain is not found in the receptor or non-standard residue is not found in the chain
        """
        print(f"Nonstandard Residue: {self.nonstandard_residue}, {type(self.nonstandard_residue)}")
        if self.nonstandard_residue is not None:
            self.protein = self.receptor.select(f"protein and chain {self.chain} and not resname {self.nonstandard_residue}")
        else:
            self.protein = self.receptor.select(f"protein and chain {self.chain}")
    
    def get_binding_pocket(self, ligand_residue_name:str, save_ligand = True, pocket_center = None):
        """Get the binding pocket of the protein. Currently, binding pocket is defined as the center of mass of the ligand.
        The residue needs to be on the chain that was previously selected.

        Args:
            ligand_residue_name (str): Residue name of the ligand
            save_ligand (bool, optional): Save the ligand as a separate PDB file. Defaults to True.

        Returns:
            ligand_center (np.array): Center of mass of the ligand (x,y,z) coordinates
        """
        if not pocket_center:
            selection_command = f"resname {ligand_residue_name} and chain {self.chain}"
            ligand = self.receptor.select(selection_command)
            self.default_ligand = ligand
            self.ligand_center = ligand.getCoords().mean(axis=0)
            logging.info(f"Center of Mass of Ligand:{self.ligand_center}")
        else:
            self.ligand_center = pocket_center
        self.random_path_number = random.randint(0,1000)
        pdb_path = f"temp/{self.random_path_number}.pdb"
        logging.info(f"Saving the protein only at: {pdb_path}")
        prody.writePDB(pdb_path, self.protein)
        if save_ligand:
            pdb_path = f"temp/{self.random_path_number}L.pdb"
            logging.info(f"Saving the ligand only at: {pdb_path}")
            prody.writePDB(pdb_path, ligand)
        return self.ligand_center
    
    @staticmethod
    def minimize_energy(fixer):
        forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

        system = forcefield.createSystem(fixer.topology, nonbondedMethod=app.CutoffNonPeriodic,
                nonbondedCutoff=0.9*unit.nanometer)

        atom_elements = [atom.element.name for atom in fixer.topology.atoms()]
        for i in range(system.getNumParticles()):
            if atom_elements[i]!='hydrogen':
                system.setParticleMass( i, 0.0 )
        integrator = openmm.LangevinIntegrator(298*unit.kelvin, 1/unit.picosecond, 1*unit.femtosecond)
        platform = openmm.Platform.getPlatformByName('CPU')

        simulation = app.Simulation(fixer.topology, system, integrator, platform)
        simulation.context.setPositions(fixer.positions)
        simulation.minimizeEnergy()
        positions = simulation.context.getState(getPositions=True).getPositions()
        return fixer, positions

    def fix_protein(self):
        pdb_path = f"temp/{self.random_path_number}.pdb"
        fixer = pdbfixer.PDBFixer(pdb_path)
        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(7.0)
        fixer, positions = ProteinPreparation.minimize_energy(fixer)
        pdb_path = f"temp/{self.random_path_number}H.pdb"
        app.PDBFile.writeFile(fixer.topology, positions, open(f'{pdb_path}', 'w'))
        logging.info(f"Fixed Protein saved at: {pdb_path}")
        return True



if __name__ == "__main__":
    fname = prody.fetchPDB("2B7A", folder = "temp/")
    pdb = prody.parsePDB(fname)
    docking_experiment = ProteinPreparation(pdb, "A", "PTR")
    binding_pocket_center = docking_experiment.get_binding_pocket("IZA")
    condition = docking_experiment.fix_protein()
    print(condition)