import os

# Each worker performs serial BLAS/LAPACK work. Keeping one BLAS thread per
# worker prevents CPU oversubscription. These variables must be set before
# importing NumPy and SciPy.
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

from pathlib import Path
from multiprocessing import Pool, freeze_support

import numpy as np
from scipy.linalg import eigvalsh
from scipy.special import expit
from tqdm import tqdm

from bulk_hamiltonian import Hamiltonian

np.set_printoptions(linewidth=200, suppress=True, precision=5)


# =============================================================================
# Numerical parameters
# =============================================================================

num_k       = 50 # number of k-space points
num_q       = 19  # number of q-space points
n_orb       = 8  # orbital number
kT          = 0.001 # temperature scale
mu          = 0.056   # chemical potential
omega       = 0.0   # energy 
eta         = 0.002 # broadening factor
a           = 3.477 # Lattice constants
b           = 6.249 # Lattice constants


K_CHUNK_SIZE = 512  # this is used to save RAM, 
                    # the code sums 512 k point once, then again 512,... 
                    # not all the k points at once

# Number of worker processes. You can instead use a fixed value, for example 8.
N_PROCESSES = 8 #max(1, (os.cpu_count() or 1) - 1)

# Number of q points sent to one worker in each Pool request.
POOL_CHUNKSIZE = 1

output_file = Path(f"wte2_sus_multi_core_mu_{mu*1000}_k_{num_k}_q_{num_q}.dat")


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
            Hamiltonian(kx_value, ky_value),
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


# =============================================================================
# Multiprocessing helpers
# =============================================================================

# The q-independent arrays are stored once inside every worker process.
worker_kx = None
worker_ky = None
worker_energies_k = None
worker_occupations_k = None


def initialize_worker(kx, ky, energies_k, occupations_k):
    """Store the q-independent arrays once in every worker process."""

    global worker_kx, worker_ky
    global worker_energies_k, worker_occupations_k

    worker_kx = kx
    worker_ky = ky
    worker_energies_k = energies_k
    worker_occupations_k = occupations_k


def calculate_one_q(task):
    """Calculate one q point and return its original array index."""

    q_index, qx_value, qy_value = task

    chi_q = lindhard_at_q(
        qx_value,
        qy_value,
        worker_kx,
        worker_ky,
        worker_energies_k,
        worker_occupations_k,
    )

    return q_index, chi_q


def main():
    kx, ky = rectangular_mesh(num_k)
    qx, qy = rectangular_mesh(num_q)

    num_k_points = kx.size
    num_q_points = qx.size
    n_workers = min(N_PROCESSES, num_q_points)

    print(f"k mesh: {num_k + 1} x {num_k + 1} "
          f"= {num_k_points} points")
    print(f"q mesh: {num_q + 1} x {num_q + 1} "
          f"= {num_q_points} points")
    print(f"Worker processes: {n_workers}")

    print("Precomputing E_n(k) ...")
    energies_k = eigenvalues_at_points(kx / a, ky / b)
    occupations_k = fermi_function(energies_k)

    susceptibility = np.empty(num_q_points, dtype=np.complex128)

    # Each task contains only the index and one q point. The larger k arrays
    # are sent once to each worker through initialize_worker.
    tasks = (
        (q_index, qx_value, qy_value)
        for q_index, (qx_value, qy_value) in enumerate(zip(qx, qy))
    )

    print("Calculating chi(q) ...")
    with Pool(
        processes=n_workers,
        initializer=initialize_worker,
        initargs=(kx, ky, energies_k, occupations_k),
    ) as pool:
        completed_results = pool.imap_unordered(
            calculate_one_q,
            tasks,
            chunksize=POOL_CHUNKSIZE,
        )

        for q_index, chi_q in tqdm(
            completed_results,
            total=num_q_points,
            desc="q points",
            unit="q",
        ):
            # Results finish in an arbitrary order. q_index restores the
            # original q-mesh ordering in the output file.
            susceptibility[q_index] = chi_q

    result = np.column_stack(
        (qx / a, qy / b, susceptibility.real, susceptibility.imag)
    )

    header = (
        "qx_over_a qy_over_b Re_chi Im_chi\n"
        f"num_k={num_k}, "
        f"num_q={num_q}, n_orb={n_orb}, "
        f"kT={kT}, mu={mu}, omega={omega}, eta={eta}, "
        f"N_PROCESSES={n_workers}"
    )
    np.savetxt(output_file, result, fmt="%.16e", header=header)
    print(f"Data saved to {output_file.resolve()}")


if __name__ == "__main__":
    # Essential on Windows; harmless on Linux and macOS.
    freeze_support()
    main()


