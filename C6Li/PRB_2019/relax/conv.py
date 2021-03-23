from ase import Atoms
from gpaw import GPAW, PW, FermiDirac
from ase.eos import EquationOfState as EOS
from ase.build import graphene, add_adsorbate
from ase.build.supercells import make_supercell
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.use('Agg')

natom = 7

def my_model(a,c):
    C2 = graphene(a=a)
    slab = make_supercell(C2,[[2,1,0],[-1,1,0],[0,0,1]])
    x, y = np.dot([1/3, 1/3],slab.cell[:2,:2])
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

def encut_conv():
    atoms = my_model(2.471,3.7)
    encuts = np.arange(300, 860, 50)
    energies = []
    for encut in encuts:
        atoms.calc = calc_scf(encut,10,0.01)
        energies.append(atoms.get_potential_energy()/natom)
    plt.plot(encuts, energies,'o-', lw=2)
    plt.xlabel('wavefunction cutoff (eV)', size = 15)
    plt.ylabel('total energy per atom (eV)', size = 15)
    plt.tight_layout()
    plt.savefig('Encut.svg')
    plt.close()

def k_conv():
    atoms = my_model(2.471,3.7)
    ks = np.arange(5.5, 12.1, 0.5)
    energies = []
    for k in ks:
        atoms.calc = calc_scf(600,k,0.01)
        energies.append(atoms.get_potential_energy()/natom)
    plt.plot(ks, energies,'o-', lw=2)
    plt.xlabel(r'Kmesh density (1/$\mathrm{\AA}$)', size = 15)
    plt.ylabel('total energy per atom (eV)', size = 15)
    plt.tight_layout()
    plt.savefig('k.svg')
    plt.close()

def width_conv():
    atoms = my_model(2.471,3.7)
    widths = np.arange(0.001, 0.052, 0.005)
    energies = []
    for width in widths:
        atoms.calc = calc_scf(600,10,width)
        atoms.get_potential_energy()
        gs_contents = [line for line in open('gs.out') if line.strip()]
        for line in gs_contents:
            if 'Entropy ' in line:
                energy = float(line.split()[-1])/natom
        energies.append(energy)
    plt.plot(widths, energies,'o-', lw=2)
    plt.xlabel('MP smearing width (eV)', size = 15)
    plt.ylabel('Entropy (-ST) per atom (eV)', size = 15)
    plt.tight_layout()
    plt.savefig('width.svg')
    plt.close()

def main():
    encut_conv()
    k_conv()
    width_conv()
    atoms = my_model(2.471,3.70)
    atoms.write('POSCAR',vasp5=True)

if __name__=="__main__":
    main()
