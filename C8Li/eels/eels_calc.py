from ase import Atoms
from gpaw import GPAW, PW, FermiDirac
from gpaw.response.df import DielectricFunction
from gpaw.bztools import find_high_symmetry_monkhorst_pack
from gpaw.mpi import world
import numpy as np
from ase.parallel import paropen
from ase.io import read
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.use('Agg')
import os

atoms = read('../relax/POSCAR')
def calc_scf(encut=600,k=10,width=0.01):
    calc = GPAW(
             mode=PW(encut),
             xc='LDA',
             kpts={'density':k,'gamma':True},
             random=True,
             occupations=FermiDirac(width),
             txt='gs.out'
               )
    return calc

def gs_response(calc,k=30,nbands=90,width=0.001):
#    kpts = find_high_symmetry_monkhorst_pack(calc,density=k, pbc=True)
    responseGS = calc.fixed_density(
#                kpts=kpts,
                kpts={'density':k},
                parallel={'band':1},
                nbands=nbands,
                occupations=FermiDirac(width),
                convergence={'bands':60},
                txt='gsresponse.out'
                    )
    return responseGS

def eels_calc(calc, direction='M', eta=25e-3, pbc=True, ecut=100, domega0=0.01, xc='RPA'):
    f = paropen('%s/%s/q_list.dat' % (xc, direction),'w')
    Q={'M':42, 'K':42, 'A':52}
    for i in range(0,3*int(Q[direction])+1,2):
        if os.path.isfile('%s/%s/EELS_%d' % (xc, direction, i)): 
            continue
        df = DielectricFunction(calc,eta=eta,pbc=pbc,ecut=ecut,
                            domega0=domega0,
    #                        integrationmode='tetrahedron_integration',
                            txt='%s/%s/df_%d.out' % (xc, direction, i))  
        if direction == 'M':
            q_c = np.array([i/42.0, 0.0, 0.0])
        elif direction == 'K':
            q_c = np.array([i/42.0, i/42.0, 0.0])
        elif direction == 'A':
            q_c = np.array([0.0, 0.0, i/52.0])
        df.get_dielectric_function(q_c=q_c,xc=xc, filename='%s/%s/epsilon_%d.csv' % (xc, direction, i))
        df.get_eels_spectrum(q_c=q_c,xc=xc, filename='%s/%s/EELS_%d' % (xc, direction, i))

        # Calculate cartesian momentum vector:
        cell_cv = atoms.cell
        bcell_cv = 2*np.pi*np.linalg.inv(cell_cv).T
        q_v = np.dot(q_c, bcell_cv)
        print(np.sqrt(np.inner(q_v, q_v)),file=f)
    f.close()

def main():
#    calc = calc_scf()
#    atoms.calc = calc
#    atoms.get_potential_energy()
#    responseGS = gs_response(calc)
#    responseGS.write('gsresponse.gpw','all')
    
    for xc in ['RPA', 'ALDA']:
        for direction in ['A', 'M', 'K']:
            if not os.path.isdir('%s/%s' %(xc, direction)):
                os.system('mkdir -p %s/%s' % (xc, direction))
            eels_calc(calc='/opt/gpw_files/C8Li/gsresponse.gpw', xc=xc, direction=direction)

if __name__=="__main__":
    main()
