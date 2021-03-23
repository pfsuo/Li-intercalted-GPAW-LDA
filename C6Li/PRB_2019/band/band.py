import numpy as np
from gpaw import GPAW
from ase.parallel import paropen

atoms = GPAW('gs.gpw').atoms
path = atoms.cell.bandpath('GKMGAHLA,KH,ML',npoints=600)
calc = GPAW('gs.gpw').fixed_density(symmetry='off',
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
