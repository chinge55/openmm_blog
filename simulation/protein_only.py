from openmm.app import *
from openmm import *
from openmm import unit

# or more explicitly:
# from openmm.app import Simulation, PDBFile, ForceField
# from openmm import LangevinIntegrator, Platform
# import openmm.unit as unit

pdb = PDBFile('../results/protein.pdb')
forcefield = ForceField('amber14/protein.ff14SB.xml')
modeller = Modeller(pdb.topology, pdb.positions)

system = forcefield.createSystem(modeller.topology)
integrator = LangevinIntegrator(300*unit.kelvin, 1.0/unit.picosecond, 0.004*unit.picoseconds)
simulation = Simulation(modeller.topology, system, integrator)
simulation.context.setPositions(modeller.positions)

print("Setup complete!")
