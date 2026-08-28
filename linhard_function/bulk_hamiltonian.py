import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigvalsh




def PauliMatrices():
    p0 = np.eye(2, dtype=complex)
    p1 = np.array([[0, 1], [1, 0]], dtype=complex)
    p2 = np.array([[0, -1j], [1j, 0]], dtype=complex)
    p3 = np.array([[1, 0], [0, -1]], dtype=complex)
    return p0, p1, p2, p3

def nmat(s):
    return np.zeros((s, s), dtype=complex)

def smat(aa, s):
    return np.block([[aa, nmat(s)], [nmat(s), aa]])
p0, p1, p2, p3 = PauliMatrices()

s0 = lambda aa: smat(aa, 2)
s1 = lambda aa: np.block([[nmat(2), aa], [aa, nmat(2)]])
s2 = lambda aa: -1j * np.block([[nmat(2), aa], [-aa, nmat(2)]])
s3 = lambda aa: np.block([[aa, nmat(2)], [nmat(2), -aa]])

s0c = lambda aa: smat(aa, 4)
s1c = lambda aa: np.block([[nmat(4), aa], [aa, nmat(4)]])
s2c = lambda aa: -1j * np.block([[nmat(4), aa], [-aa, nmat(4)]])
s3c = lambda aa: np.block([[aa, nmat(4)], [nmat(4), -aa]])

def ArrayFlatten(matrix):
    return np.block(matrix)

g0 = ArrayFlatten(s0(p0))
g1p = ArrayFlatten((s0(p0) + s0(p3)) / 2)
g1m = ArrayFlatten((s0(p0) - s0(p3)) / 2)
g2p = ArrayFlatten((s1(p0) + 1j * s2(p0) + s1(p3) + 1j * s2(p3)) / 4)
g2m = ArrayFlatten((s1(p0) + 1j * s2(p0) - s1(p3) - 1j * s2(p3)) / 4)
g3 = ArrayFlatten((1j * s1(p2) - s2(p2)) / 2)
g4p = ArrayFlatten((s0(p1) + s3(p1) + 1j * s0(p2) + 1j * s3(p2)) / 4)
g4m = ArrayFlatten((s0(p1) - s3(p1) + 1j * s0(p2) - 1j * s3(p2)) / 4)
g5p = ArrayFlatten((s3(p0) + s3(p3)) / 2)
g5m = ArrayFlatten((s3(p0) - s3(p3)) / 2)
g6 = ArrayFlatten((s1(p1) + 1j * s2(p1)) / 2)


a = 3.477
b = 6.249

rad = np.array([-0.25 * a, 0.32 * b])
rap = np.array([-0.25 * a, -0.07 * b])
rbp = np.array([0.25 * a, 0.07 * b])
rbd = np.array([0.25 * a, -0.32 * b])

tpx = 1.13
tdx = -0.41
tpab = 0.40
tdab = 0.51
t0ab = 0.39
t0abx = 0.29
t0x = 0.14
tpy = 0.13

mup = -1.65
mud = 0.74


l0aby = 0.011
l0y = 0.051
l0z = 0.012
l0py = 0.050
l0pz = 0.012
lpxy = -0.040
lpxz = -0.010
ldxy = -0.031
ldxz = -0.008



def Hamiltonian(kx, ky):

    d1 = np.dot([kx, ky], (rad - rbd))
    d2 = np.dot([kx, ky], (rap - rbp))
    d3 = np.dot([kx, ky], (rad - rbp))
    d4 = np.dot([kx, ky], (rad - rap))


    hsoc = (
        (ldxz * s3c(g5p) + ldxy * s2c(g5p)) * np.sin(a * kx) +
        (lpxz * s3c(g5m) + lpxy * s2c(g5m)) * np.sin(a * kx) -
        1j * l0aby * s2c(g6) * (1 + np.exp(1j * a * kx)) * np.exp(1j * d3) -
        (1j * (l0z * s3c(g4p) + l0y * s2c(g4p)) * np.exp(1j * d4) -
         1j * (l0z * s3c(g4m) + l0y * s2c(g4m)) * np.exp(-1j * d4)) -
        (1j * (l0pz * s3c(g4p) + l0py * s2c(g4p)) * np.exp(-1j * b * ky) * np.exp(1j * d4) -
         1j * (l0pz * s3c(g4m) + l0py * s2c(g4m)) * np.exp(1j * b * ky) * np.exp(-1j * d4))
    )

    hsoc += hsoc.T.conj()

    h0b = (
        s0c((mup / 2 + tpx * np.cos(a * kx) + tpy * np.cos(b * ky)) * g1m +
            (mud / 2 + tdx * np.cos(a * kx)) * g1p +
            tdab * np.exp(-1j * b * ky) * (1.0 + np.exp(1j * a * kx)) * np.exp(1j * d1) * g2p +
            tpab * (1.0 + np.exp(1j * a * kx)) * np.exp(1j * d2) * g2m +
            t0ab * (1.0 - np.exp(1j * a * kx)) * np.exp(1j * d3) * g3 +
            t0abx * (np.exp(-1j * a * kx) - np.exp(2j * a * kx)) * np.exp(1j * d3) * g3 -
            2j * t0x * np.sin(a * kx) * np.exp(1j * d4) * g4p -
            2j * t0x * np.sin(a * kx) * np.exp(-1j * d4) * g4m) +
        s0c((mup / 2 + tpx * np.cos(a * kx) + tpy * np.cos(b * ky)) * g1m +
            (mud / 2 + tdx * np.cos(a * kx)) * g1p +
            tdab * np.exp(-1j * b * ky) * (1.0 + np.exp(1j * a * kx)) * np.exp(1j * d1) * g2p +
            tpab * (1.0 + np.exp(1j * a * kx)) * np.exp(1j * d2) * g2m +
            t0ab * (1.0 - np.exp(1j * a * kx)) * np.exp(1j * d3) * g3 +
            t0abx * (np.exp(-1j * a * kx) - np.exp(2j * a * kx)) * np.exp(1j * d3) * g3 -
            2j * t0x * np.sin(a * kx) * np.exp(1j * d4) * g4p -
            2j * t0x * np.sin(a * kx) * np.exp(-1j * d4) * g4m).T.conj()
    )



    h = hsoc + h0b 
    return h



# =============================================================================
# k-independent 8x8 matrices
# =============================================================================

# SOC matrices
soc_1 = ldxz * s3c(g5p) + ldxy * s2c(g5p)
soc_2 = lpxz * s3c(g5m) + lpxy * s2c(g5m)

soc_3 = s2c(g6)

soc_4p = l0z * s3c(g4p) + l0y * s2c(g4p)
soc_4m = l0z * s3c(g4m) + l0y * s2c(g4m)

soc_5p = l0pz * s3c(g4p) + l0py * s2c(g4p)
soc_5m = l0pz * s3c(g4m) + l0py * s2c(g4m)


# Normal Hamiltonian matrices
hg1m = s0c(g1m)
hg1p = s0c(g1p)
hg2p = s0c(g2p)
hg2m = s0c(g2m)
hg3  = s0c(g3)
hg4p = s0c(g4p)
hg4m = s0c(g4m)


# Coordinate differences
v1 = rad - rbd
v2 = rap - rbp
v3 = rad - rbp
v4 = rad - rap



def Hamiltonian_batch(kx, ky):

    kx = np.asarray(kx)
    ky = np.asarray(ky)

    nk = kx.size

    # ---------------------------------------------------------
    # phases d1, d2, d3, d4 for every k point
    # ---------------------------------------------------------

    d1 = kx * v1[0] + ky * v1[1]
    d2 = kx * v2[0] + ky * v2[1]
    d3 = kx * v3[0] + ky * v3[1]
    d4 = kx * v4[0] + ky * v4[1]


    # ---------------------------------------------------------
    # useful k-dependent quantities
    # ---------------------------------------------------------

    sin_ax = np.sin(a * kx)

    cos_ax = np.cos(a * kx)
    cos_by = np.cos(b * ky)

    exp_ax  = np.exp(1j * a * kx)
    exp_max = np.exp(-1j * a * kx)
    exp_2ax = np.exp(2j * a * kx)

    exp_by  = np.exp(1j * b * ky)
    exp_mby = np.exp(-1j * b * ky)

    exp_d1  = np.exp(1j * d1)
    exp_d2  = np.exp(1j * d2)
    exp_d3  = np.exp(1j * d3)

    exp_d4  = np.exp(1j * d4)
    exp_md4 = np.exp(-1j * d4)


    # ---------------------------------------------------------
    # SOC Hamiltonian
    #
    # shape = (nk, 8, 8)
    # ---------------------------------------------------------

    hsoc = np.zeros((nk, 8, 8), dtype=complex)

    hsoc += sin_ax[:, None, None] * soc_1[None, :, :]

    hsoc += sin_ax[:, None, None] * soc_2[None, :, :]

    hsoc += (
        -1j
        * l0aby
        * ((1.0 + exp_ax) * exp_d3)[:, None, None]
        * soc_3[None, :, :]
    )

    hsoc += (
        -1j
        * exp_d4[:, None, None]
        * soc_4p[None, :, :]
    )

    hsoc += (
        +1j
        * exp_md4[:, None, None]
        * soc_4m[None, :, :]
    )

    hsoc += (
        -1j
        * (exp_mby * exp_d4)[:, None, None]
        * soc_5p[None, :, :]
    )

    hsoc += (
        +1j
        * (exp_by * exp_md4)[:, None, None]
        * soc_5m[None, :, :]
    )


    # Add Hermitian conjugate
    hsoc = hsoc + hsoc.conj().transpose(0, 2, 1)


    # ---------------------------------------------------------
    # Normal Hamiltonian
    # ---------------------------------------------------------

    h0 = np.zeros((nk, 8, 8), dtype=complex)


    # g1m
    coefficient = (
        mup / 2
        + tpx * cos_ax
        + tpy * cos_by
    )

    h0 += coefficient[:, None, None] * hg1m[None, :, :]


    # g1p
    coefficient = (
        mud / 2
        + tdx * cos_ax
    )

    h0 += coefficient[:, None, None] * hg1p[None, :, :]


    # g2p
    coefficient = (
        tdab
        * exp_mby
        * (1.0 + exp_ax)
        * exp_d1
    )

    h0 += coefficient[:, None, None] * hg2p[None, :, :]


    # g2m
    coefficient = (
        tpab
        * (1.0 + exp_ax)
        * exp_d2
    )

    h0 += coefficient[:, None, None] * hg2m[None, :, :]


    # first g3 term
    coefficient = (
        t0ab
        * (1.0 - exp_ax)
        * exp_d3
    )

    h0 += coefficient[:, None, None] * hg3[None, :, :]


    # second g3 term
    coefficient = (
        t0abx
        * (exp_max - exp_2ax)
        * exp_d3
    )

    h0 += coefficient[:, None, None] * hg3[None, :, :]


    # g4p
    coefficient = (
        -2j
        * t0x
        * sin_ax
        * exp_d4
    )

    h0 += coefficient[:, None, None] * hg4p[None, :, :]


    # g4m
    coefficient = (
        -2j
        * t0x
        * sin_ax
        * exp_md4
    )

    h0 += coefficient[:, None, None] * hg4m[None, :, :]


    # Add Hermitian conjugate
    h0 = h0 + h0.conj().transpose(0, 2, 1)


    # ---------------------------------------------------------
    # Complete Hamiltonian
    # ---------------------------------------------------------

    h = hsoc + h0

    return h


# # Set ky to 0 and kx from -pi to pi with 400 divisions
# kx_values = np.linspace(-np.pi, np.pi, 200)
# ky = 0


# mu = 0.06
# def cp(mu):
#     return mu * np.eye(8)
# # Initialize a list to store eigenvalues
# eigenvalues = []

# # Diagonalize Hamiltonian for each kx and store the eigenvalues
# for kx in kx_values:
#     h = Hamiltonian(kx/a, ky/b) - cp(mu)
#     eigvals = eigvalsh(h)*1000  # Use eigvalsh for Hermitian matrices
#     eigenvalues.append(eigvals)

# # Convert eigenvalues list to a numpy array for easier manipulation
# eigenvalues = np.array(eigenvalues)


# plt.rcParams["figure.figsize"] = [5, 4]
# plt.rcParams["figure.autolayout"] = True
# plt.rcParams["font.family"] = "serif"
# # plt.rcParams["font.serif"] = ["Times New Roman"]
# plt.rcParams["mathtext.fontset"] = "dejavuserif"  # or "stix"

# # Plot the band structure
# # plt.figure(figsize=(8, 4))
# for band in range(eigenvalues.shape[1]):
#     plt.plot(kx_values/a, eigenvalues[:, band], color = 'b', lw = 2)
# plt.tick_params(direction='in', length=7, width=1, colors='k', bottom=True,
#                     top=True, left=True, right=True)
# plt.xlabel(r'$k_x$ (1/$\AA$)', fontsize = 18)
# plt.ylabel('Energy (meV)', fontsize = 18)
# plt.ylim(-40,40)
# # plt.xlim(-0.6,0.6)
# # plt.text(-.58, -30, f"$\mu$ = {mu} meV", fontsize = 18)
# # plt.axhline(y = 0, ls = '--', color ='k')
# # plt.xticks([-0.5,0,0.5],fontsize=18)
# # plt.yticks([-40,-20,0,20,40],fontsize=18)
# plt.savefig(f"band.png", dpi=300)
# # plt.savefig(f"band.pdf", dpi=300)
# plt.show()
