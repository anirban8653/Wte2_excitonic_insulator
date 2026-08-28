import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
from pathlib import Path
import numpy as np
# from scipy.linalg import eigvalsh
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.pyplot as plt
from bulk_hamiltonian import Hamiltonian_batch
np.set_printoptions(linewidth=200, suppress=True, precision=5)


# =============================================================================
# Numerical parameters
# =============================================================================

num_k       = 100 # number of k-space points
num_q       = 30  # number of q-space points
n_orb       = 8  # orbital number
kT          = 0.001 # temperature scale
mu          = 0.06   # chemical potential
omega       = 0.00   # energy 
eta         = 0.005 # broadening factor
a           = 3.477 # Lattice constants
b           = 6.249 # Lattice constants


output_file = Path(f"wte2_sus_single_core_mu_{mu*1000}_k_{num_k}_q_{num_q}.dat")


# def cp(mu):
#     return mu * np.eye(n_orb)

def rectangular_mesh(num_divisions):
    """This function creates a two-dimensional square momentum mesh"""
    axis = np.linspace(-np.pi, np.pi, num_divisions + 1, endpoint=True)
    x_grid, y_grid = np.meshgrid(axis, axis, indexing="ij")
    return x_grid.ravel(order="C"), y_grid.ravel(order="C")


def fermi_function(x):
    return 1.0/(np.exp(x/kT)+1.0)

# def eigenvalues_at_points(kx, ky):
#     """Diagonalize the Hamiltonian at paired arrays of physical k points."""

#     energies = np.empty((kx.size, n_orb), dtype=np.float64)
#     for index, (kx_value, ky_value) in enumerate(zip(kx, ky)):
#         energies[index] = eigvalsh(Hamiltonian(kx_value/a, ky_value/b)- cp(mu))
#     return energies


kx, ky = rectangular_mesh(num_k)
qx, qy = rectangular_mesh(num_q)

num_k_points = kx.size
num_q_points = qx.size

H_k = Hamiltonian_batch(kx/a, ky/b)
energies_k = np.linalg.eigvalsh(H_k) - mu
occupations_k = fermi_function(energies_k)

results = []
for q_index, (qx_value, qy_value) in enumerate(zip(qx, qy)):

    # susceptibility = 0.0 + 0.0j

    kx_q = kx + qx_value
    ky_q = ky + qy_value

    # All k+q Hamiltonians
    H_kq = Hamiltonian_batch(kx_q/a, ky_q/b)

    # All eigenvalues
    energies_kq = np.linalg.eigvalsh(H_kq) - mu

    # All occupations
    occupations_kq = fermi_function(energies_kq)

    # All k, i, j combinations
    numerator = (
        occupations_kq[:, :, None]
        - occupations_k[:, None, :]
    )

    denominator = (
        omega
        + energies_kq[:, :, None]
        - energies_k[:, None, :]
        + 1j*eta
    )

    susceptibility = -np.sum(numerator / denominator)/num_k_points
    results.append([susceptibility.real, susceptibility.imag])
    print(q_index, np.round(susceptibility.real,6), np.round(susceptibility.imag,6))

results = np.array(results)
full_results = np.column_stack((qx/a, qy/b, results))

header = (
        "qx_over_a qy_over_b Re_chi Im_chi\n"
        f"num_k={num_k}, "
        f"num_q={num_q}, n_orb={n_orb}, "
        f"kT={kT}, mu={mu}, omega={omega}, eta={eta}"
    )
np.savetxt(output_file, full_results, fmt="%.16e", header=header)
print(f"Data saved to {output_file.resolve()}")