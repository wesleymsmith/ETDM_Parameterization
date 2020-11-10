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


bsse=False
bssebasis="None"

if bsse:
    psi4.core.print_out('INTERACTION ENERGY correction BSSE\n')
else:
    psi4.core.print_out('INTERACTION ENERGY correction None\n')

crdlist=importlib.import_module('.'.join(os.path.basename(__file__).split('.')[:-1])+'_xyz')
basename=crdlist.basename
intrange=crdlist.intrange
basecoor=getattr(crdlist,'basecoor')
elist={}
for i in intrange:
    importname=basename+str(i)
    intecoor=getattr(crdlist,importname)
    mol=psi4.geometry(basecoor+"--"+intecoor)
    mol.update_geometry() # This update is required for psi4 to load the molecule
    E=psi4.energy(opttheory+'/'+optbasis)
    key=float('.'.join(i.split('_')))
    elist[key]=E
mindist=min(elist, key=elist.get)
mine=elist[mindist]

if bsse and optbasis != bssebasis:
    importname=basename+'_'.join(str(mindist).split('.'))
    intecoor=getattr(crdlist,importname)
    mol=psi4.geometry(basecoor+"--"+intecoor)
    mol.update_geometry()
    mine=psi4.energy(opttheory+'/'+bssebasis,bsse_type='cp')
    
psi4.core.print_out('INTERACTION DISTANCE and ENERGY are: '+str(mindist)+' '+str(mine)+'\n')

psi4.core.print_out('INTERACTION TABLE NOBSSE START\n')
for d in elist.keys():
    psi4.core.print_out(str(d)+' '+str(elist[d])+'\n')
psi4.core.print_out('INTERACTION TABLE NOBSSE END\n')

