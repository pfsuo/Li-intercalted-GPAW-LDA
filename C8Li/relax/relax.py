from ase import Atoms
from gpaw import GPAW, PW, FermiDirac
from ase.eos import EquationOfState as EOS
from ase.build import graphene, add_adsorbate
from ase.build.supercells import make_supercell
from ase.optimize.sciopt import SciPyFminBFGS as BFGS
import numpy as np
import matplotlib.pyplot as plt

natom = 9

def my_model(a,c):
    C2 = graphene(a=a)
    slab = make_supercell(C2,[[2,0,0],[0,2,0],[0,0,1]])
    x, y = np.dot([1/6, 1/3],slab.cell[:2,:2])
    add_adsorbate(slab,'Li',position=(x,y),height=c/2)
    slab.cell[2] = [0, 0, c]
    slab.pbc = True
    return slab

def calc_scf(encut,k,width):
    calc = GPAW(
             mode=PW(encut),
             xc='LDA',
             kpts={'density':k,'gamma':True},
             random=True,
             occupations=FermiDirac(width),
             txt='gs.out'
               )
    return calc

def a_fit(c):
    a_s = np.arange(2.40,2.61,0.02)
    Es = []
    Vs = []
    for a in a_s:
        atoms = my_model(a,c)
        atoms.calc = calc_scf(600,11,0.01)
        relax = BFGS(atoms)
        relax.run(fmax=0.01)
        Vs.append(atoms.get_volume())
        Es.append(atoms.get_potential_energy())
    eos = EOS(Vs, Es, eos='birchmurnaghan')
    v0, e0, B = eos.fit()
    eos.plot(filename='a_fit.svg')
    plt.close()
    return np.sqrt(v0/3*2/np.sqrt(3)/c)

def main():
    a = a_fit(3.7)
    atoms = my_model(a, 3.7)
    calc = calc_scf(600,11,0.01)
    atoms.calc = calc 
    relax = BFGS(atoms,logfile='relax.log')
    relax.run(fmax=0.01)
    atoms.write('POSCAR',vasp5=True)
    atoms.get_potential_energy()
    calc.write('gs.gpw','all')

if __name__=="__main__":
    main()
