# Creates: EELS.png
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.use('Agg')
from ase.io import read
import sys

prefix = sys.argv[1]
def get_data():
    f = open('q_list.dat','w')
    for i in range(1,42,2):
        q_c = [i/42.0, i/42.0, 0.0]  # Γ-A excitation
        cell_cv = read('../../gsresponse.out').cell
        bcell_cv = 2*np.pi*np.linalg.inv(cell_cv).T
        q_v = np.dot(q_c, bcell_cv)
        print(np.sqrt(np.inner(q_v, q_v)),file=f)
    f.close()
    energy = np.loadtxt(prefix + '_EELS_1', delimiter=',')[:,0]
    np.savetxt('energy.dat', energy)
    Q = np.loadtxt('q_list.dat')
    eels_LFE = np.zeros([int(len(energy)),int(len(Q))])
    eels_NLFE = np.zeros([int(len(energy)),int(len(Q))])
    for i, q in enumerate(Q):
        filename = prefix + '_EELS_' + str(i*2 + 1)
        d = np.loadtxt(filename, delimiter=',')
        eels_NLFE[:,i] = d[:,1]
        eels_LFE[:,i] = d[:,2]
    np.savetxt('eels_NLFE.dat', eels_NLFE)
    np.savetxt('eels_LFE.dat', eels_LFE)

def plot_eels(figsize=(6,5)):
    fontsize=10
    energy = np.loadtxt('energy.dat')[:803]
    q_list = np.loadtxt('q_list.dat')
    xticks = [q_list[0],q_list[int(int(len(q_list)/2)-1)],q_list[-1]]
    xticklabels = ['', r'$\mathrm{2\pi/a\sqrt{3}}$', r'$\mathrm{4\pi/a\sqrt{3}}$']
    eels_LFE = np.loadtxt('eels_LFE.dat')[:803,:]
    eels_LFE /= eels_LFE.max()
    eels_NLFE = np.loadtxt('eels_NLFE.dat')[:803,:]
    eels_NLFE /= eels_LFE.max()
    x, y = np.meshgrid(q_list, energy)
    
    norm = mpl.colors.Normalize(0,1)
    s_m = mpl.cm.ScalarMappable(cmap='turbo', norm=norm)
    s_m.set_array([eels_LFE])
    from mpl_toolkits.axes_grid1 import host_subplot
    ax = host_subplot(111)
    plt.contourf(x,y,eels_LFE,np.arange(0,1.01,.002), cmap='turbo',norm=norm)
    cbar = plt.colorbar()
    cbar.set_ticks([0,1])
    plt.xlabel(r'Momentum ($\mathrm{\AA}^{-1}$)', size=15)
    plt.ylabel('Energy (eV)', size=15)
    plt.xticks(size=fontsize)
    plt.yticks(size=fontsize)
    ax2 = ax.twin()
    ax2.set_xticks(xticks)
    ax2.set_xticklabels(xticklabels)
    ax2.axis["right"].major_ticklabels.set_visible(False)
    ax2.axis["right"].major_ticks.set_visible(False)
    ax2.axis['top']
    plt.savefig('EELS_LFE.png',dpi=600)
    plt.close()
    
    ax = host_subplot(111)
    plt.contourf(x,y,eels_NLFE,np.arange(0,1.01,.002), cmap='turbo',norm=norm)
    cbar = plt.colorbar()
    cbar.set_ticks([0,1])
    plt.xlabel(r'Momentum ($\mathrm{\AA}^{-1}$)', size=15)
    plt.ylabel('Energy (eV)', size=15)
    plt.xticks(size=fontsize)
    plt.yticks(size=fontsize)
    ax2 = ax.twin()
    ax2.set_xticks(xticks)
    ax2.set_xticklabels(xticklabels)
    ax2.axis["right"].major_ticklabels.set_visible(False)
    ax2.axis["right"].major_ticks.set_visible(False)
    ax2.axis['top']
    plt.savefig('EELS_NLFE.png',dpi=600)
    plt.close()

def main():
    get_data()
    plot_eels()

if __name__=="__main__":
    main()
