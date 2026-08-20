import numpy as np
import matplotlib.pyplot as plt

plt.rcParams["figure.figsize"] = [8, 12]
plt.rcParams["figure.autolayout"] = True
plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"  # or "stix"

mu_data = np.array([0.29,0.28,0.27,0.26,0.25,0.24,0.23,0.22,0.21,0.2,0.19,0.18])

for j in range(len(mu_data)):
    en, ldos = np.loadtxt(f'DOS-{mu_data[j]}.dat',unpack=True)

    # for i in range(1,7):

    plt.plot(en, ldos + j*0.15, color = 'k',linewidth=2)
    plt.text(-0.22, -0.03 + j*0.15, f' -{mu_data[j]} eV', fontsize=20)
    plt.hlines(0 + j*0.15, -0.03 - j * 0.01 ,0.13- j * 0.01 ,color='black', lw=1, ls ='dotted')

plt.axvline(x=0, color='black', linewidth=1, ls= '--')
# plt.ylim(0, 0.2)
# plt.text(-0.25, 0.18, f'W = {W} eV', fontsize=15)
plt.xlim(-0.22,0.25)
# plt.legend(loc=(0.22,0.52),frameon=False,fontsize=15)
plt.xlabel('$\omega$ (eV)', fontsize=35)
plt.ylabel('DOS', fontsize=35)
plt.tick_params(direction='in', length=7, width=1, colors='k', bottom=True,
                top=True, left=True, right=True)
plt.yticks(fontsize=30)
plt.xticks(fontsize=30)
plt.savefig("DOS_fountain_plot.png", dpi=300)
plt.savefig("DOS_fountain_plot.pdf", dpi = 300)
plt.show()