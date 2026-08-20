import numpy as np
import matplotlib.pyplot as plt

mu_data = np.array([-0.29,-0.28,-0.27,-0.26,-0.25,-0.24,-0.23,-0.22,-0.21,-0.2,-0.19, -0.18])
W_data = np.array([-0.64,-0.64,-0.64,-0.64,-0.64,-0.64,-0.64,-0.64,-0.64,-0.64,-0.39, -0.32])

for j in range(len(mu_data)):
    en, ldos = np.loadtxt(f'bulk_mu{mu_data[j]}_wmax0.25Ny12Nx12eta0.0007Ns50Ldos_wte2.dat', unpack=True)
    
    plt.rcParams["figure.figsize"] = [6, 6]
    plt.rcParams["figure.autolayout"] = True
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["mathtext.fontset"] = "dejavuserif"  # or "stix"
    
    # plt.subplot(1,2,2)
    # for i in range(0,48,2):
    ny = 12  # y length
    nx = 12
    
    lyr = ny  # number of layers for the topograph
    ival = len(np.array([i for i in range(int(nx * ny * 4/2 ))]))
    nw = 401
    x = en[0::ival] 
    
    y_all_cell_spin_up = []
    for k in range(0, ival, 4):
        y_cell_spin_up = sum([ldos[i::ival] for i in range(k + 0, k + 4)]) / 4
        # y_cell_spin_up = ldos[k+1::ival]
        y_all_cell_spin_up.append(y_cell_spin_up)
    y_all_cell_spin_up = np.array(y_all_cell_spin_up)
    
    for i in range(0, len(y_all_cell_spin_up)):
        x = np.append(x, y_all_cell_spin_up[i])
    # print(len(y))
    topo_data = x.reshape(int(nx * ny/2 ) + 1, nw).T
    # print(topo_data)
    np.savetxt(f"unitcell_avg_ldos_data_sm_2pi6_mu{mu_data[j]}.dat", topo_data, fmt="%12.6f")
    # %%
    
    # data = np.loadtxt(f'unitcell_avg_ldos_data_sm_2pi6_mu{mu}.dat')
    # for i in range(1,7):
    # # for i in range(1,7):
    #     plt.plot(data[:, 0], data[:, i],linewidth=2,label=f'{i}')
    #     plt.axvline(x=0, color='black', linewidth=0.8)
    #     # plt.ylim(0, 0.2)
    #     # plt.xlim(0.05,0.15)
    #     plt.legend(loc=(0.22,0.52),frameon=False,fontsize=15)
    #     plt.xlabel('$\omega$ (meV)', fontsize=15)
    #     plt.ylabel('DOS', fontsize=15)
    #     plt.tick_params(direction='in', length=7, width=1, colors='k', bottom=True,
    #                     top=True, left=True, right=True)
    #     plt.yticks(fontsize=15)
    #     plt.xticks(fontsize=15)
    # plt.savefig(f"xe_dos_OP_sm_mu{mu}.png", dpi=300)
    # plt.show()
    
    #%%
    
    data = np.loadtxt(f'unitcell_avg_ldos_data_sm_2pi6_mu{mu_data[j]}.dat')
    
    # for i in range(1,7):
    data_avg = (data[:, 1]+data[:, 2]+data[:, 3]+data[:, 4]+data[:, 5]+data[:,6])/6
    plt.plot(data[:, 0], data_avg,color = 'blue',linewidth=3)
    plt.axvline(x=0, color='black', linewidth=0.8)
    plt.ylim(0, 0.25)
    plt.text(-0.25, 0.23, f'W = {W_data[j]} eV', fontsize=15)
    plt.text(-0.25, 0.21, f'$\mu$ = {mu_data[j]} eV', fontsize=15)
    # plt.xlim(0.05,0.15)
    plt.legend(loc=(0.22,0.52),frameon=False,fontsize=15)
    plt.xlabel('$\omega$ (eV)', fontsize=15)
    plt.ylabel('DOS', fontsize=15)
    plt.tick_params(direction='in', length=7, width=1, colors='k', bottom=True,
                    top=True, left=True, right=True)
    plt.yticks(fontsize=15)
    plt.xticks(fontsize=15)
    plt.savefig(f"DOS{mu_data[j]}.png", dpi=300)
    plt.show()
    plot_data = []
    for i in range(0,len(data_avg)):
        plot_data.append([data[:,0][i],data_avg[i]])
    
    np.savetxt(f'DOS{mu_data[j]}.dat',np.array(plot_data))
