from ase import Atoms
from gpaw import GPAW, PW, FermiDirac
from ase.eos import EquationOfState as EOS
from ase.build import graphene, add_adsorbate
from ase.build.supercells import make_supercell
from ase.optimize.sciopt import SciPyFminBFGS as BFGS
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.use('Agg')

natom = 8

def my_model(a,c=3.7):
    C2 = graphene(a=a)
    slab = make_supercell(C2,[[2,1,0],[-1,1,0],[0,0,1]])
    x1, y1 = np.dot([1/3, 1/3],slab.cell[:2,:2])
    x2, y2 = np.dot([2/3, 0],slab.cell[:2,:2])
    add_adsorbate(slab,'Li',position=(x1,y1),height=c/2)
    add_adsorbate(slab,'Li',position=(x2,y2),height=c/2)
    slab.cell[2] = [0, 0, c]
    slab.pbc = True
    return slab

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

def a_fit(c=3.7):
    a_s = np.arange(2.44, 2.63, 0.02)
    Vs = []
    Es = []
    for a in a_s:
        atoms = my_model(a, c)
        atoms.calc = calc_scf()
        Vs.append(atoms.get_volume())
        relax = BFGS(atoms, logfile='relax_fit.log')
        relax.run(fmax=0.01)
        Es.append(atoms.get_potential_energy())
    eos = EOS(Vs, Es, eos='birchmurnaghan')
    v0, e0, B = eos.fit()    
    eos.plot(filename='a_fit.svg')
    plt.close()
    return np.sqrt(v0*2/np.sqrt(3)/c)/np.sqrt(3)

def main():
    a = a_fit()
    atoms = my_model(a, 3.7)
    atoms.calc = calc_scf(600,10,0.01)
    relax = BFGS(atoms,logfile='relax.log')
    relax.run(fmax=0.01)
    atoms.write('POSCAR',vasp5=True)
    atoms.get_potential_energy()

if __name__=="__main__":
    main()
