import numpy as np
import matplotlib.pyplot as plt

mu_data = np.array([0.29,0.28,0.27,0.26,0.25,0.24,0.23,0.22,0.21,0.2])
eg_data = []
for j in range(len(mu_data)):
    en, ldos = np.loadtxt(f'DOS-{mu_data[j]}.dat',unpack=True)
    data = []
    for i in range(len(ldos)):
        if ldos[i] < 0.02:
            data.append(en[i])

    data = np.array(data)
    eg = (data[-1]-data[0])

    eg_data.append([-mu_data[j], data[-1], data[0], eg])

eg_data = np.array(eg_data)

plt.rcParams["figure.figsize"] = [6, 6]
plt.rcParams["figure.autolayout"] = True
plt.rcParams["font.family"] = "serif"
# plt.rcParams["font.serif"] = ["Times New Roman"]
plt.rcParams["mathtext.fontset"] = "dejavuserif"  # or "stix"

plt.plot(eg_data[:,0], eg_data[:,1], '-^', color='blue', label="CBM", markersize=10, markerfacecolor='none')
plt.plot(eg_data[:,0], eg_data[:,2], '-v', color='green', label="VBM", markersize=10, markerfacecolor='none')
plt.plot(eg_data[:,0], eg_data[:,3], '-o', color='magenta', label="Gap", markersize=10, markerfacecolor='none')
plt.axhline(y=0, ls='--', color='k', lw='1')

plt.legend(loc=(0.05, 0.0), frameon=False, fontsize=20)
plt.ylim(-0.15,0.15)
plt.xlabel("$\mu$ (eV)", fontsize=30)
plt.ylabel("$\omega$ (eV)", fontsize=30)

plt.tick_params(direction='in', length=7, width=1, colors='k', bottom=True,
                top=True, left=True, right=True)
plt.yticks(fontsize=20)
plt.xticks([-0.29,-0.26,-0.23,-0.2],fontsize=20)
plt.savefig('gap_variation.png', dpi=300)
plt.savefig('gap_variation.pdf', dpi = 300)
plt.show()
