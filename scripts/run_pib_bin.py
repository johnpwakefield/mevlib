#!/usr/bin/env python3


from struct import unpack

import numpy as np
from scipy.linalg import solve_triangular
from scipy.interpolate import interp1d
from scipy.integrate import solve_ivp

from mevlib.outfmt.binary import mgn_ints, mgn_diag

import matplotlib.pyplot as plt


#TODO should this script just be scrapped?


plt.rc('font', size=10)
plt.rc('text', usetex=True)
plt.rc('axes', labelsize=10)
plt.rc('legend', fontsize=10)

axsize = (3.0, 2.5)


np.set_printoptions(precision=4)


# constants

intsfile = "example_sph_mevtable_ints.dat"
diagfile = "example_sph_mevtable_diag.dat"

SYMBLEN = 8     # length of species identifier string
voidage = 0.4
VOL = 1.0

y0 = [0.8, 0.0, 0.0, 0.0, 0.0, 0.0]

T = 600.0


# read binary files

with open(intsfile, 'rb') as f:
    mgn = f.read(8)
    if mgn != mgn_ints:
        raise Exception("{} is not a diag file.".format(diagfile))
    M, = unpack('<i', f.read(4))
    f.read(4)
    tmods = np.fromfile(f, dtype=np.float64, count=M, offset=0)
    seffs = np.fromfile(f, dtype=np.float64, count=M, offset=0)
    if len(f.read(1)) != 0:
        raise Exception("Extra bytes at end of file.")

sslookup = interp1d(tmods, seffs)

with open(diagfile, 'rb') as f:
    mgn = f.read(8)
    if mgn != mgn_diag:
        raise Exception("{} is not a diag file.".format(diagfile))
    Ng, = unpack('<B', f.read(1))
    Ns, = unpack('<B', f.read(1))
    M, = unpack('<H', f.read(2))
    f.read(6)
    names = [
        ''.join(
            map(lambda x: x.decode('utf-8'), unpack('<8c', f.read(8)))
        ).strip()
        for i in range(Ng + Ns)
    ]
    molwghts = np.fromfile(f, dtype=np.float64, count=(Ng+Ns), offset=0)
    temps = np.fromfile(f, dtype=np.float64, count=M, offset=0)
    Ds, lams, Rs, BsRs = [], [], [], []
    for i in range(M):
        Ds.append(np.fromfile(f, dtype=np.float64, count=Ng+Ns, offset=0))
        lams.append(np.fromfile(f, dtype=np.float64, count=Ng, offset=0))
        Rs.append(np.fromfile(f, dtype=np.float64, count=Ng*Ng, offset=0))
        BsRs.append(np.fromfile(f, dtype=np.float64, count=Ns*Ng, offset=0))
    Rs = list(map(lambda x: x.reshape((Ng, Ng)).T, Rs))
    BsRs = list(map(lambda x: x.reshape((Ng, Ns)).T, BsRs))
    if len(f.read(1)) != 0:
        raise Exception("Extra bytes at end of file.")


# setup equations

def rhs(t, m):
    i = np.searchsorted(temps, T) + 1
    Zfree = solve_triangular(Rs[i], m[:-1], lower=True)
    if deac:
        yk = m[-1] / cat_mass   # this is not a mass fraction, but is what was used in the CEJ paper
        psi = (
            (1 + YK_coeff * 100 * yk)**(-1.6)
            / (1 + KA * (WA + WR + WASP))
            / (1 + KN * (WN / WCO) * t)
        )
        assert(0 < psi and psi < 1)
    else:
        psi = 1.0
    mult = np.array([1.0 if lam == 0.0 else sslookup(lam) for lam in lams[i]])
    Mdot = -psi * voidage * VOL * Ds[i] * np.dot(
        np.vstack((np.dot(Rs[i], np.diag(lams[i])), BsRs[i])),
        mult * Zfree
    )
    return Mdot


# run pib

soln = solve_ivp(rhs, [0.0, 1200.0], y0, dense_output=True, atol=1e-9)


# make plots

fig, ax = plt.subplots(1, 1, figsize=axsize)
logvals = np.exp(np.linspace(np.log(0.1), np.log(1000.0), 200))
for i, name in enumerate(names):
    #ax.plot(soln.t, soln.y[i, :], label=name)
    solninterp = interp1d(soln.t, soln.y[i, :], kind='cubic')
    ax.semilogx(logvals, [solninterp(t) for t in logvals], label=name)
ax.set_xlabel(r'\( t \)')
ax.grid()
fig.tight_layout()      # call before legend
ax.legend(loc="upper center", ncol=3)


for ext in ["pdf", "svg"]:
    fig.savefig("img/pib_bin.{}".format(ext))


