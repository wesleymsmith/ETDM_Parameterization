#!/usr/bin/env python
import psi4,os,importlib
import numpy as np
mem="1000mb"
cpu=4
psi4.set_num_threads(cpu)
psi4.set_memory(mem)

psi4.core.set_output_file(os.path.dirname(os.path.abspath(__file__))+'/'+'.'.join(os.path.basename(__file__).split('.')[:-1])+'.out', False)




opttheory="hf"
optbasis="6-31g*"
options={'scf_type': 'df', 'g_convergence':'gau','freeze_core':'true'}
psi4.set_options(options)


psi4_xyz="""
0 1
h    	-4.858650	-2.687319	0.927118
o    	-4.556517	-3.595501	0.914798
h    	-3.627853	-3.536450	0.690459

"""
mol=psi4.geometry(psi4_xyz)
mol.update_geometry() # This update is required for psi4 to load the molecule
E=psi4.energy(opttheory+"/"+optbasis)
psi4.core.print_out('INTERACTION MONOMER ENERGY is : '+str(E) )