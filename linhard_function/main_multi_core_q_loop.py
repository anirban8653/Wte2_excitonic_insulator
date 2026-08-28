import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

from pathlib import Path
from multiprocessing import Pool, cpu_count
import numpy as np
from scipy.linalg import eigvalsh
from tqdm import tqdm

from bulk_hamiltonian import Hamiltonian_batch

np.set_printoptions(linewidth=200, suppress=True, precision=5)

# =============================================================================
# Numerical parameters
# =============================================================================

num_k       = 100   # number of k-space divisions
num_q       = 50    # number of q-space divisions
n_orb       = 8     # orbital number
kT          = 0.001 # temperature scale
mu          = 0.06  # chemical potential
omega       = 0.00  # energy
eta         = 0.005 # broadening factor
a           = 3.477 # lattice constant
b           = 6.249 # lattice constant

# Number of worker processes
# You can change this manually if you want
n_processes = 50 #cpu_count()

output_file = Path(f"wte2_sus_parallel_mu_{mu*1000}_k_{num_k}_q_{num_q}.dat")

# =============================================================================
# Helper functions
# =============================================================================

def cp(mu):
    return mu * np.eye(n_orb)

def rectangular_mesh(num_divisions):
    """Create a two-dimensional square momentum mesh."""
    axis = np.linspace(-np.pi, np.pi, num_divisions + 1, endpoint=True)
    x_grid, y_grid = np.meshgrid(axis, axis, indexing="ij")
    return x_grid.ravel(order="C"), y_grid.ravel(order="C")

def fermi_function(x):
    return 1.0 / (np.exp(x / kT) + 1.0)


# def eigenvalues_at_points(kx, ky):
#     """Diagonalize the Hamiltonian at paired arrays of physical k points."""
#     energies = np.empty((kx.size, n_orb), dtype=np.float64)
#     for index, (kx_value, ky_value) in enumerate(zip(kx, ky)):
#         energies[index] = eigvalsh(Hamiltonian(kx_value / a, ky_value / b) - cp(mu))
#     return energies

# =============================================================================
# Worker globals
# These will be initialized once per worker process
# =============================================================================

worker_kx = None
worker_ky = None
worker_energies_k = None
worker_occupations_k = None
worker_num_k_points = None

def initialize_worker(kx, ky, energies_k, occupations_k):
    """Store q-independent arrays once inside each worker process."""
    global worker_kx, worker_ky
    global worker_energies_k, worker_occupations_k
    global worker_num_k_points

    worker_kx = kx
    worker_ky = ky
    worker_energies_k = energies_k
    worker_occupations_k = occupations_k
    worker_num_k_points = kx.size

def calculate_one_q(task):
    """
    Compute susceptibility for one q-point.
    task = (q_index, qx_value, qy_value)
    """

    q_index, qx_value, qy_value = task

    # All k+q points
    kx_q = worker_kx + qx_value
    ky_q = worker_ky + qy_value

    # All H(k+q)
    H_kq = Hamiltonian_batch(kx_q / a, ky_q / b)

    # All eigenvalues
    energies_kq = np.linalg.eigvalsh(H_kq) - mu

    # All occupations
    occupations_kq = fermi_function(energies_kq)

    # All k, i, j combinations
    numerator = (
        occupations_kq[:, :, None]
        - worker_occupations_k[:, None, :]
    )

    denominator = (
        omega
        + energies_kq[:, :, None]
        - worker_energies_k[:, None, :]
        + 1j * eta
    )

    susceptibility = (
        -np.sum(numerator / denominator)
        / worker_num_k_points
    )

    return (
        q_index,
        qx_value / a,
        qy_value / b,
        susceptibility.real,
        susceptibility.imag,
    )
# =============================================================================
# Main
# =============================================================================

if __name__ == "__main__":

    # Build meshes
    kx, ky = rectangular_mesh(num_k)
    qx, qy = rectangular_mesh(num_q)

    num_k_points = kx.size
    num_q_points = qx.size

    print(f"Number of k-points = {num_k_points}")
    print(f"Number of q-points = {num_q_points}")
    print(f"Using {n_processes} worker processes")

    # Compute q-independent quantities once
    H_k = Hamiltonian_batch(kx/a, ky/b)
    energies_k = np.linalg.eigvalsh(H_k) - mu
    occupations_k = fermi_function(energies_k)

    # Prepare q tasks
    q_tasks = [
        (q_index, qx_value, qy_value)
        for q_index, (qx_value, qy_value) in enumerate(zip(qx, qy))
    ]

    # Parallel q-loop
    results = []
    with Pool(
        processes=n_processes,
        initializer=initialize_worker,
        initargs=(kx, ky, energies_k, occupations_k),
    ) as pool:

        for result in tqdm(
            pool.imap_unordered(calculate_one_q, q_tasks, chunksize=1),
            total=num_q_points,
            desc="Computing susceptibility"
        ):
            results.append(result)

    # Sort back into original q-order
    results.sort(key=lambda x: x[0])

    # Print results
    for q_index, qx_scaled, qy_scaled, re_chi, im_chi in results:
        print(q_index, np.round(re_chi, 6), np.round(im_chi, 6))

    # Convert to array for saving
    full_results = np.array([
        [qx_scaled, qy_scaled, re_chi, im_chi]
        for _, qx_scaled, qy_scaled, re_chi, im_chi in results
    ])

    header = (
        "qx_over_a qy_over_b Re_chi Im_chi\n"
        f"num_k={num_k}, "
        f"num_q={num_q}, n_orb={n_orb}, "
        f"kT={kT}, mu={mu}, omega={omega}, eta={eta}, "
        f"n_processes={n_processes}"
    )

    np.savetxt(output_file, full_results, fmt="%.16e", header=header)
    print(f"Data saved to {output_file.resolve()}")