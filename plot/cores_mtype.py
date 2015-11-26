import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

colortipo = ['brown', 'red', 'orange', 'green', '#00D0C9', '#0076C9', 'blue']

N = 1000
Ntype = len(colortipo)
mtypes = [ -1, 0, 1, 2, 3, 4, 5 ]
mtype_labels = [ 'E', 'S0', 'Sa', 'Sb', 'Sbc', 'Sc', 'Sd' ]

x = np.random.rand(N)
y = np.random.rand(N)
# mtype below will simul a random -1 ... >5 array
z_mtype = np.random.randint(Ntype, size = N) - 1 
f = plt.figure()
ax = f.gca()
cmap = mpl.colors.ListedColormap(colortipo)
sc = ax.scatter(x, y, c = z_mtype, cmap = cmap, marker = 'o', s = 10, edgecolor = 'none', label = '')
ax.set_xlabel(r'x', fontsize = 15) 
ax.set_ylabel(r'y', fontsize = 15)
tickpos = np.linspace(mtypes[0] + .5, mtypes[-1] - .5, Ntype)
cb = f.colorbar(sc, ticks = tickpos)
cb.ax.set_yticklabels(mtype_labels)
f.savefig('cores_mtype.png')
