#!/usr/bin/env python
import psi4,os,importlib
import numpy as np
mem="1000mb"
cpu=4
psi4.set_num_threads(cpu)
psi4.set_memory(mem)

psi4.core.set_output_file(os.path.dirname(os.path.abspath(__file__))+'/'+'.'.join(os.path.basename(__file__).split('.')[:-1])+'.out', False)




opttheory="mp2"
optbasis="6-31g*"
options={'scf_type':'df','g_convergence':'gau','freeze_core':'true','opt_coordinates':'both',"frozen_dihedral":"20 10 4 1"}
psi4.set_options(options)


crdlist=importlib.import_module('.'.join(os.path.basename(__file__).split('.')[:-1])+'_xyz')
basename=crdlist.basename
scanrange=crdlist.scanrange
psi4.core.print_out('SCAN ATOMS: 20 10 4 1\n')
for i in scanrange:
    psi4.core.print_out('STARTING '+i+'\n')
    importname=basename+str(i)
    psi4_xyz=getattr(crdlist,importname)
    mol=psi4.geometry(psi4_xyz)
    mol.update_geometry() # This update is required for psi4 to load the molecule
    psi4.set_options(options)
    try:
        psi4.optimize(opttheory+'/'+optbasis,molecule=mol)
    except psi4.OptimizationConvergenceError:
        tmpopt=options
        tmpopt['g_convergence']='gau_loose'
        psi4.set_options(tmpopt)
        try:
           psi4.optimize(opttheory+'/'+optbasis,molecule=mol)
        except psi4.OptimizationConvergenceError:
           psi4.core.print_out('ENDING ERROR '+i+'\n')
        else:
           E=psi4.energy(opttheory+'/'+optbasis)
           psi4.core.print_out('ENERGY is: '+str(E)+'\n')
           psi4.core.print_out('ENDING SUCCESS '+i+'\n')
    else:
        E=psi4.energy(opttheory+'/'+optbasis)
        psi4.core.print_out('ENERGY is: '+str(E)+'\n')
        psi4.core.print_out('ENDING SUCCESS '+i+'\n')
