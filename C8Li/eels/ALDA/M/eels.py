from gpaw import GPAW, PW, FermiDirac
from gpaw.response.df import DielectricFunction
from gpaw.mpi import world
import numpy as np
from ase.parallel import paropen
from ase.io import read

f = paropen('C8Li_q_list','w')
for i in range(23,42,2):
    df = DielectricFunction('../../gsresponse.gpw',eta=25e-3,pbc=True,ecut=100,
                        domega0=0.01,
#                        integrationmode='tetrahedron_integration',
                        txt='df_%d.out' % i)  
    q_c = [i/42.0, 0.0, 0.0]  # Γ-A excitation
    df.get_eels_spectrum(q_c=q_c, xc='ALDA', filename='C8Li_EELS_%d' % i)

    # Calculate cartesian momentum vector:
    cell_cv = read('../../gsresponse.out').cell
    bcell_cv = 2*np.pi*np.linalg.inv(cell_cv).T
    q_v = np.dot(q_c, bcell_cv)
    print(np.sqrt(np.inner(q_v, q_v)),file=f)
f.close()
