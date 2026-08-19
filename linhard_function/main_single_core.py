import os

# Limit threaded BLAS because this version parallelizes neither the q loop nor
# the diagonalization.  These variables must be set before importing NumPy and
# SciPy.
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

from pathlib import Path

import numpy as np
from scipy.linalg import eigvalsh
from scipy.special import expit

from bulk_hamiltonian import Hamiltonian

np.set_printoptions(linewidth=200, suppress=True, precision=5)


# =============================================================================
# Numerical parameters
# =============================================================================

num_k       = 20 # number of k-space points
num_q       = 9  # number of q-space points
n_orb       = 8  # orbital number
kT          = 0.001 # temperature scale
mu          = 0.056   # chemical potential
omega       = 0.0   # energy 
eta         = 0.002 # broadening factor
a           = 3.477 # Lattice constants
b           = 6.249 # Lattice constants


K_CHUNK_SIZE = 512

output_file = Path(f"wte2_sus_single_core_mu_{mu*1000}_k_{num_k}_q_{num_q}.dat")


def cp(mu):
    return np.diag([mu] * n_orb)

def rectangular_mesh(num_divisions):
    """This function creates a two-dimensional square momentum mesh"""
    axis = np.linspace(-np.pi, np.pi, num_divisions + 1, endpoint=True)
    x_grid, y_grid = np.meshgrid(axis, axis, indexing="ij")
    return x_grid.ravel(order="C"), y_grid.ravel(order="C")


def fermi_function(energies) :
    """This function calculates fermi energy"""
    return expit(-(energies - mu) / kT)


def eigenvalues_at_points(kx, ky):
    """Diagonalize the Hamiltonian at paired arrays of physical k points."""

    energies = np.empty((kx.size, n_orb), dtype=np.float64)
    for index, (kx_value, ky_value) in enumerate(zip(kx, ky)):
        energies[index] = eigvalsh(
            Hamiltonian(kx_value, ky_value), #- cp(mu),
            lower=False,
            overwrite_a=True,
            check_finite=False,
        )
    return energies


def lindhard_at_q(qx, qy, kx, ky, energies_k, occupations_k):
    """-1/Nk * sum_{k,i,j}[f(E_i(k+q)) - f(E_j(k))]/ [omega + E_i(k+q) - E_j(k) + i*eta]"""

    num_k_points = kx.size
    susceptibility = 0.0j

    for start in range(0, num_k_points, K_CHUNK_SIZE):
        stop = min(start + K_CHUNK_SIZE, num_k_points)
        selection = slice(start, stop)

        kx_q = (kx[selection] + qx) / a
        ky_q = (ky[selection] + qy) / b

        energies_kq = eigenvalues_at_points(kx_q, ky_q)
        occupations_kq = fermi_function(energies_kq)

        numerator = (occupations_kq[:, :, np.newaxis] - occupations_k[selection, np.newaxis, :])
        denominator = ( omega + energies_kq[:, :, np.newaxis] - energies_k[selection, np.newaxis, :] + 1j * eta)

        susceptibility -= np.sum(numerator / denominator)

    return susceptibility / num_k_points



kx, ky = rectangular_mesh(num_k)
qx, qy = rectangular_mesh(num_q)

num_k_points = kx.size
num_q_points = qx.size

print(f"k mesh: {num_k + 1} x {num_k + 1} "
        f"= {num_k_points} points")
print(f"q mesh: {num_q + 1} x {num_q + 1} "
        f"= {num_q_points} points")

print("Precomputing E_n(k) ...")

energies_k = eigenvalues_at_points(kx / a, ky / b)
occupations_k = fermi_function(energies_k)

susceptibility = np.empty(num_q_points, dtype=np.complex128)


print("Calculating chi(q) ...")
for q_index, (qx_value, qy_value) in enumerate(zip(qx, qy)):
    susceptibility[q_index] = lindhard_at_q(
        qx_value,
        qy_value,
        kx,
        ky,
        energies_k,
        occupations_k,
    )

    print(q_index)

# Same four columns written by the Fortran program.
result = np.column_stack(
    (qx / a, qy / b, susceptibility.real, susceptibility.imag)
)

header = (
        "qx_over_a qy_over_b Re_chi Im_chi\n"
        f"num_k={num_k}, "
        f"num_q={num_q}, n_orb={n_orb}, "
        f"kT={kT}, mu={mu}, omega={omega}, eta={eta}"
    )
np.savetxt(output_file, result, fmt="%.16e", header=header)
print(f"Data saved to {output_file.resolve()}")



