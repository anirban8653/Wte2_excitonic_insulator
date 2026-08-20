import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
# import matplotlib.cm as cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

mu_data = np.array([0.29,0.28,0.27,0.26,0.25,0.24,0.23,0.22,0.21,0.2,0.19,0.18])
eg_data = []
for j in range(len(mu_data)):
    en, ldos = np.loadtxt(f'DOS-{mu_data[j]}.dat',unpack=True)
    eg_data.append(ldos[::-1])

eg_data = np.array(eg_data).T





interp = 'bilinear'


def forceaspect(ax, aspect=1):
    im = ax.get_images()
    extent = im[0].get_extent()
    ax.set_aspect(abs((extent[1] - extent[0]) / (extent[3] - extent[2])) / aspect)







colors = ['navy', 'dodgerblue', 'red','white']
cmap_name = 'my_list'
cm = LinearSegmentedColormap.from_list(
    cmap_name, colors, N=1000)
norm = mpl.colors.Normalize(vmin=np.min(eg_data), vmax=np.max(eg_data))

plt.rcParams["figure.figsize"] = [6, 6]
plt.rcParams["figure.autolayout"] = True
plt.rcParams["font.family"] = "serif"
# plt.rcParams["font.serif"] = ["Times New Roman"]
plt.rcParams["mathtext.fontset"] = "dejavuserif"  # or "stix"
fig = plt.figure()


x = -mu_data
y = en


ax = fig.add_subplot(111)
im = plt.imshow(eg_data, interpolation=interp, extent=(np.amin(x), np.amax(x), np.amin(y), np.amax(y)),cmap='plasma', norm=norm)
forceaspect(ax, aspect = 1)
plt.xlabel('$\mu$ (eV)', fontsize=30)
plt.ylabel('$\omega$ (eV)', fontsize=30)
plt.xticks([-0.28,-0.25,-0.22,-0.19],fontsize=20)
plt.yticks(fontsize=20)
plt.tick_params(direction='in', length=5, width=1, colors='k', bottom=True, top=False, left=True, right=False)


cbar =plt.colorbar(im,shrink=0.77,  location = 'top')
cbar.ax.tick_params(labelsize=15)
# plt.savefig(f"topo_ldos_sum_all_density_xe_OP_pi_sm_2pi6_mu{mu}.png", dpi=600)
plt.savefig("topo_DOS.png", dpi=300)
plt.savefig("topo_DOS.pdf", dpi =300)


# Show the plot
plt.show()
