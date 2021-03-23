import numpy as np
from gpaw import GPAW, PW, FermiDirac
from ase.parallel import paropen
from ase.io import read

atoms = read('../relax/POSCAR')
calc = GPAW(mode=PW(600),
            xc='LDA',
            kpts={'density':10,'gamma':True},
            occupations = FermiDirac(0.01),
            txt='gs.out'
            )
atoms.calc = calc
atoms.get_potential_energy()
#calc.write('gs.gpw','all')

path = atoms.cell.bandpath('GKMGAHLA,KH,ML',npoints=600)
calc = calc.fixed_density(symmetry='off',
                        kpts=path.kpts,
                        txt='band.out')
calc.write('band.gpw')
x, X, _ = path.get_linear_kpoint_axis()
with paropen('kpath.dat','w') as f:
    for k in x:
        print(k,file=f)

with paropen('highk.dat','w') as f:
    for k in X:
        print(k,file=f)
ef = calc.get_fermi_level()
e_kn = calc.band_structure().energies[0]
e_nk = e_kn.T - ef
np.savetxt('e_nk.dat',e_nk)
