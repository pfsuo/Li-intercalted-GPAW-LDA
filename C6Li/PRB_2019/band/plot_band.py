import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.use('Agg')
from gpaw import GPAW

title = 'Bandstructure of C$_6$Li'
calc = GPAW('band.gpw')
ef = calc.get_fermi_level()
nbnd = calc.get_number_of_bands()
natom = calc.atoms.get_global_number_of_atoms()
NEDOS = 3001
kx = np.loadtxt('kpath.dat')
X = np.loadtxt('highk.dat')
e_nk = np.loadtxt('e_nk.dat')
klabels = ['Γ','K','M','Γ','A','H','L','','A|K','','H|M','L']
w_kni = abs(calc.get_projections(locfun='projectors'))
energy, dos = calc.get_dos(npts=NEDOS,width=0.1)

def plotband(figsize=(5,5)):
    plt.figure(figsize=figsize)
    for e_k in e_nk:
        plt.plot(kx,e_k,c='r',lw=2)
    for k in X[1:-1]:
        plt.axvline(x=k,c='k',lw=0.5,ls='--')
    plt.axhline(y=0,c='k',lw=0.5,ls='--')
    plt.xticks(X,klabels,size=15)
    plt.yticks(size=15)
    plt.ylabel(r'$\varepsilon_n(k) - \varepsilon_{\mathrm{F}}$ (eV)',size=15)
    plt.title(title, size=18)
    plt.axis([0,kx[-1],-5,5])
    plt.tight_layout()
    plt.savefig('band.svg' )
    plt.close()

def plotdos(figsize=(5,5)):
    plt.figure(figsize=figsize)
    plt.plot(energy - ef, dos, c = 'r', lw=2)
    plt.xlabel(r'$\varepsilon_n(k) - \varepsilon_{\mathrm{F}}$ (eV)',size=15)
    plt.ylabel('Density of States (1/eV)',size=15)
    plt.xlim([-8,8])
    plt.ylim(ymin=0)
    plt.tight_layout()
    plt.savefig('DOS.svg' )
    plt.close()

def plotbanddos(figsize=(9,5)):
    plt.figure(figsize=figsize)
    grid = plt.GridSpec(1,5)
    p1 = plt.subplot(grid[0,0:4])
    for e_k in e_nk:
        plt.plot(kx,e_k,c='r',lw=2)
    for k in X[1:-1]:
        plt.axvline(x=k,c='k',lw=0.5,ls='--')
    plt.axhline(y=0,c='k',lw=0.5,ls='--')
    plt.xticks(X,klabels,size=15)
    plt.yticks(size=15)
    plt.ylabel(r'$\varepsilon_n(k) - \varepsilon_{\mathrm{F}}$ (eV)',size=15)
    plt.title(title, size=18)
    plt.axis([0,kx[-1],-5,5])
    p2 = plt.subplot(grid[0,4])
    plt.plot(dos, energy - ef, c='r',lw=2)
    plt.axhline(y=0, lw=0.5, c='k',ls='--')
    plt.axvline(x=0, lw=0.5, c='k',ls='--')
    plt.fill_between(dos,energy - ef, 0, where=dos>=0, facecolor='silver',interpolate=True)
    plt.xlim(xmin=0)
    plt.ylim([-5,5])
    plt.xlabel('DOS (a.u.)',size=15)
    plt.xticks([])
    plt.ylabel('')
    plt.yticks([])
    plt.tight_layout()
    plt.savefig('banddos.svg' )
    plt.close()

def plotpband(figsize=(8,5),fmt='png'):
    plt.figure(figsize=figsize)
    for e_k in e_nk:
        plt.plot(kx,e_k,c='0.5',lw=0.5)
    for k in X[1:-1]:
        plt.axvline(x=k,c='k',lw=0.5,ls='--')
    w = np.zeros([len(kx),nbnd,2])
    w[:,:,0] = w_kni[:,:,2:23:4].sum(axis=2)
    w[:,:,1] = w_kni[:,:,24]
    scale = 60.0
    colors = ['cyan','r']
    labels = ['C (pz)', 'Li (s)']
    st = []
    for i in range(len(colors)):
        st.append(plt.scatter(-1,-1,20,c=colors[i],alpha=0.5,label=labels[i],marker='.',edgecolor='none'))
        for n in range(nbnd):
            st.append(plt.scatter(kx,e_nk[n,],w[:,n,i].T*scale,c=colors[i],alpha=0.5,marker='.',edgecolor='none'))
    plt.axhline(y=0,c='k',lw=0.5,ls='--')
    plt.xticks(X,klabels,size=15)
    plt.yticks(size=15)
    plt.ylabel(r'$\varepsilon_n(k) - \varepsilon_{\mathrm{F}}$ (eV)',size=15)
    plt.title(title, size=18)
    plt.legend(scatterpoints =1, numpoints=1,markerscale=2.0,fontsize=15,loc='best') 
    plt.xlim([0,kx[-1]])
    plt.ylim([-10,10])
    plt.tight_layout()
    plt.savefig('pband_whole.' + fmt )
    plt.ylim([-5,5])
    plt.tight_layout()
    plt.savefig('pband.' + fmt)
    plt.close()

def plotpbanddos(figsize=(10,5),fmt='png'):
    plt.figure(figsize=figsize)
    w = np.zeros([len(kx),nbnd,4])
    w[:,:,0] = w_kni[:,:,0:21:4].sum(axis=2)
    for i in range(0,21,4):
        w[:,:,1] += w_kni[:,:,i+1:i+4].sum(axis=2)
    w[:,:,2] = w_kni[:,:,24] 
    w[:,:,3] = w_kni[:,:,25:28].sum(axis=2)
    pdos_whole = np.zeros([NEDOS,2*natom])
    count = 0
    for i in 'sp':
        for j in range(natom):
            _, pdos_whole[:,count] = calc.get_orbital_ldos(a=j,angular=i,npts=NEDOS)
            count += 1
    pdos = np.zeros([NEDOS,4])
    pdos[:,0] = pdos_whole[:,:6].sum(axis=1)
    pdos[:,1] = pdos_whole[:,7:13].sum(axis=1)
    pdos[:,2] = pdos_whole[:,6]
    pdos[:,3] = pdos_whole[:,13]
    scale = 60.0
    colors = ['cyan','r']
    labels = ['C (p)', 'Li (s)']
    st = []
    grid = plt.GridSpec(1,5)
    p1 = plt.subplot(grid[0,0:4])
    for e_k in e_nk:
        plt.plot(kx,e_k,c='0.5',lw=0.5)
    for k in X[1:-1]:
        plt.axvline(x=k,c='k',lw=0.5,ls='--')
    for i in range(len(colors)):
        st.append(plt.scatter(-1,-1,20,c=colors[i],alpha=0.5,label=labels[i],marker='.',edgecolor='none'))
        for n in range(nbnd):
            st.append(plt.scatter(kx,e_nk[n,],w[:,n,i+1].T*scale,c=colors[i],alpha=0.5,marker='.',edgecolor='none'))
    plt.axhline(y=0,c='k',lw=0.5,ls='--')
    plt.xticks(X,klabels,size=15)
    plt.yticks(size=15)
    plt.ylabel(r'$\varepsilon_n(k) - \varepsilon_{\mathrm{F}}$ (eV)',size=15)
    plt.title(title, size=18)
#    plt.legend(scatterpoints =1, numpoints=1,markerscale=2.0,fontsize=15)
    plt.axis([0,kx[-1],-5,5])

    p2 = plt.subplot(grid[0,4])
    for i in range(len(colors)):
        plt.plot(pdos[:,i+1],energy - ef, c=colors[i],lw=2, label=labels[i])
    plt.xlim(xmin=0)
    plt.ylim([-5,5])
    plt.xlabel('DOS (a.u.)',size = 15)
    plt.xticks([])
    plt.ylabel('')
    plt.yticks([])
    plt.legend(fontsize=15,loc='best')
    plt.tight_layout()
    plt.savefig('pbanddos.' + fmt)
    plt.close()

def main():
    plotband((8,5))
    plotdos()
    plotbanddos((10,5))
    plotpband()
    plotpbanddos((10,5))

if __name__=='__main__':
    main()
