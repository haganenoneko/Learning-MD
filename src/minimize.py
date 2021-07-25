from simtk.openmm.app import *
from simtk.openmm import * 
from simtk.unit import * 
from sys import stdout 

import glob 
import os 

def MinimizePDB(input_pdb, out_path, iters=100):
    """"
    https://www.sciencedirect.com/science/article/pii/S0969212620301726
    CHARMM36 force field with CHARMM 36m modifications 
    1. 5000 steps of minimization with steepest descent
    2. Heat to 300K with random velocities using Berendsen thermostat = 1fs for 50ps 
    3. Couple semi-isotropically to Berendsen barostat 
    4. Equilibrate at 1fs for 25ps, and 2fs for 200ps
        Gradually relax restrains on protein heavy atoms: backbone (4000 to 500), side chains (2000 to 200; kJ/mol)
    5. Production = NVT, 2fs for ?, 300K. 
        Restraints on protein heavy atoms (backbone=500, side chains=200; kJ/mol)
        
    LINCS used to constrain H bonds 
    PME for long range elecotrstatics
    """"
    fname = os.path.basename(input_pdb)
    print('Minimizing... < %s >' % fname)
    
    # create system 
    system = forcefield.createSystem(modeller.topology, nonbondedMethod=PME)
    integrator = VerletIntegrator(0.001*picoseconds)

    # create simulation 
    simulation = Simulation(modeller.topology, system, integrator)
    simulation.context.setPositions(modeller.positions)
    
    # assess initial energy
    state = simulation.context.getState(getEnergy=True)
    print("Initial potential energy = ", state.getPotentialEnergy())
    
    simulation.minimizeEnergy(maxIterations=iters)
    print("Potential energy after %d iterations: %.3e" % (iters, state.getPotentialEnergy()))
    
    print('Saving...')
    positions = simulation.context.getState(getPositions=True).getPositions()
    PDBFile.writeFile(simulation.topology, positions, open(out_path * 'output.pdb', 'w'))
    print('Done')
    