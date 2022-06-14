#!/usr/bin/env python3


from struct import unpack

import numpy as np
from scipy.linalg import solve_triangular
from scipy.interpolate import interp1d
from scipy.integrate import solve_ivp

from mevlib.outfmt.binary import mgn_ints, mgn_diag

import matplotlib.pyplot as plt


plt.rc('font', size=10)
plt.rc('text', usetex=True)
plt.rc('axes', labelsize=10)
plt.rc('legend', fontsize=8)
np.set_printoptions(precision=4)

axsize = (2.75, 2.75)

def figsize(i, j):
    return (axsize[0] * j, axsize[1] * i)


# constants (all length units in cubic millimeters)

shape = "sph"
intsfile = "example_{}_mevtable_ints.dat".format(shape)
diagfile = "example_{}_mevtable_diag.dat".format(shape)
SYMBLEN = 8         # length of species identifier string
GASCONST = 8.3144                                       # [Pa m^3 / (K mol)]
voidage = 0.319
VOL = (800e-6)**3                                           # [m^3]
# for cyl pVOL = (np.pi * 80.0**2 * 550.0) * 1e-18    # [m^3]
pVOL = np.pi * (400.0e-6)**3 / 6
assert(pVOL < VOL)
epsf = (VOL - pVOL) / VOL
T = 600.0                                                   # [K]
cat_density = 1160e3 * (1 - voidage)                        # [g / m^3]
cat_mass = pVOL * cat_density
tf = 1e4                                                    # [s]
ts = np.linspace(0.0, tf, 2000)                             # [s]
MWs = np.array([444.0, 230.0, 115.0, 52.0, 16.0, 400.0])    # [g / mol]
MW_air = 28.97                                              # [g / mol]
YK_coeff, KA, WA, WR, WASP, KN, WN, WCO = (
    0.49, 0.0022, 35.5, 3.35, 3.35, 2.835, 6.89e-4, 10.0
)
colors = ['blue', 'green', 'red', 'yellow', 'orange', 'brown']
m0 = np.array([0.2, 0.0, 0.0, 0.0, 0.0, 0.0])
P = 101325.0
rhog0 = P / GASCONST / T / (m0[0] / MWs[0] + (1.0 - m0[0]) / MW_air)
mass_air = rhog0 * (VOL - pVOL) * (1.0 - sum(m0))
m0 *= rhog0 * (VOL - pVOL)


# print stuff

print("\n\n")
print("Initial mol S: {}".format(m0[0] * MWs[0]))
print("Initial density: {}".format(rhog0))
print("Particle Volume: {}".format(pVOL))
print("Domain Volume: {}".format(VOL))


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

assert(all(seffs >= 0.0) and all(seffs <= 1.0))
print("minimum eff: {}, maximum eff: {}".format(min(seffs), max(seffs)))
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
    molweights = np.fromfile(f, dtype=np.float64, count=Ng+Ns, offset=0)
    if any(
        abs(mwhard - mwfile) > 1e-8 for mwhard, mwfile in zip(MWs, molweights)
    ):
        print(
            "Molecular weights in binary file do not match hard coded values."
        )
        exit(2)
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


# setup equations (units are mass instead of mass fraction)

# while temperature is constant there's no reason to keep calling the lookup
# function
ref_i = np.searchsorted(temps, T, side='right') - 1
ref_w = T - temps[ref_i]
# if we ever test at a temperature not in the table interpolation will need to
# be added below
assert(abs(ref_w) < 1e-8)
print("conservation check")
print(np.dot(
    Ds[ref_i].reshape((1, -1)),
    np.vstack((
        np.dot(Rs[ref_i], np.diag(lams[ref_i])),
        BsRs[ref_i]
    ))
))

def massfraction(m):
    Y = np.copy(m)
    Y[:-1] /= (sum(m[:-1]) + mass_air)
    Y[-1] = m[-1] / (m[-1] + cat_mass)
    return Y

def gasdensity(m):
    return (sum(m[:-1]) + mass_air) / (VOL - pVOL)

def massrate(t, m, deac):
    i = ref_i
    Zfree = solve_triangular(Rs[i], massfraction(m)[:-1], lower=True)
    if deac:
        # this is not a mass fraction, but is what was used in the CEJ paper
        yk = m[-1] / cat_mass
        psi = (
            (1 + YK_coeff * 100 * yk)**(-1.6)
            / (1 + KA * (WA + WR + WASP))
            / (1 + KN * (WN / WCO) * t)
        )
        assert(0 < psi and psi < 1)
    else:
        psi = 1.0
    EV = np.array([sslookup(psi * lam) for lam in lams[i]])
    rate = - psi * voidage * (sum(m[:-1]) + mass_air) / (VOL - pVOL) * pVOL * (
        Ds[i] * np.dot(
            np.vstack((np.dot(Rs[i], np.diag(lams[i])), BsRs[i])),
            EV * Zfree
        )
    )
    return rate


# run pib

print("initial rate: {}".format(massrate(0.0, m0, deac=False)))

fig, axs = plt.subplots(2, 2, figsize=figsize(2, 2))
fig_cons, axs_cons = plt.subplots(2, 2, figsize=figsize(2, 2))
fig_sepcke, axs_sepcke = plt.subplots(2, 2, figsize=figsize(2, 2))
fig_gasdensity, ax_gasdensity = plt.subplots(1, 1, figsize=figsize(1, 1))
fig_pibpaperwi, ax_pibpaperwi = plt.subplots(1, 1, figsize=figsize(1, 1))
fig_pibpaperno, ax_pibpaperno = plt.subplots(1, 1, figsize=figsize(1, 1))
axs_pibpaper = [ax_pibpaperno, ax_pibpaperwi]
figs_pibpaper = [fig_pibpaperwi, fig_pibpaperno]
for j in range(2):
    soln = solve_ivp(
        lambda t, m: massrate(t, m, deac=j),
        [0.0, tf], np.copy(m0), dense_output=True,
        atol=1e-11, rtol=1e-8, method='LSODA'
    )
    ts, ms = soln.t, soln.y
    deacstring = "With{} Deactivation Due to Coking".format("" if j else "out")
    print("With{} deactivation, initial mass is {} and final is {}.".format(
        "" if j else "out", np.sum(m0), np.sum(ms[:, -1])
    ))
    for i, name in enumerate(names):
        axs[j, 0].semilogx(
            ts, [massfraction(m)[i] for m in ms.T], label=name, color=colors[i]
        )
        axs_pibpaper[j].semilogx(
            ts, [1e9 * m[i] for m in ms.T], label=name, color=colors[i]
        )
        if name != "CK":
            axs[j, 1].loglog(
                ts, [
                    np.abs(massrate(t, m, deac=j)[i])
                    / MWs[i]
                    for t, m in zip(ts, ms.T)
                ],
                label=name, color=colors[i]
            )
        if name != "CK":
            axs_sepcke[j, 0].semilogx(
                ts, [1e9 * m[i] for m in ms.T],
                label=name, color=colors[i]
            )
    assert(names[-1] == "CK")
    axs_sepcke[j, 1].loglog(
        ts, [1e9 * m[-1] for m in ms.T],
        label=name, color=colors[len(names) - 1]
    )
    for k in range(2):
        axs[j, k].set_title(deacstring)
        axs[j, k].set_xlabel(r'\( t \) [s]')
        axs[j, k].grid()
        axs_sepcke[j, k].set_title(deacstring)
        axs_sepcke[j, k].set_xlabel(r'\( t \) [s]')
        axs_sepcke[j, k].grid()
    axs_pibpaper[j].set_xlabel(r'\( t \) [s]')
    axs_pibpaper[j].grid()
    ax_gasdensity.semilogx(
        ts, [gasdensity(m) / rhog0 for m in ms.T], label=deacstring
    )
    ax_gasdensity.set_ylabel("Fraction of Initial Density")
    ax_gasdensity.set_xlabel(r"\( t \) [s]")
    axs[j, 0].set_ylabel("Mass Fraction")
    axs[j, 1].set_ylabel(
        r"Reaction Rate \( \displaystyle \left[\frac{\mathrm{mol}}{\mathrm{s}}"
        r"\right] \)"
    )
    axs[j, 0].set_xlim((1e-3, tf))
    axs[j, 1].set_xlim((1e-2, 1e3))
    axs[j, 0].set_ylim((0.0, 1.0))
    #axs[j, 1].set_ylim((1e-4, 1e2))
    axs_cons[j, 0].plot(ts, 1e9 * np.sum(ms, axis=0))
    axs_cons[j, 1].plot(ts, [
        1e9 * np.sum(massrate(t, m, deac=j)) for t, m in zip(ts, ms.T)
    ])
    axs_sepcke[j, 0].set_ylabel("Mass [ng]")
    axs_sepcke[j, 1].set_ylabel("Mass [ng]")
    axs_pibpaper[j].set_xlim((1e-1, 1e4))
    axs_pibpaper[j].set_ylabel("Mass [ng]")


axs[1, 0].legend()
fig.tight_layout()
axs_sepcke[1, 0].legend()
fig_sepcke.tight_layout()
ax_gasdensity.grid()
ax_gasdensity.legend()
for j in range(2):
    figs_pibpaper[j].tight_layout()      # call before legend
    axs_pibpaper[j].legend(loc="upper center", ncol=3)


for ext in ["pdf", "svg"]:
    fig.savefig("img/pibode_{}.{}".format(shape, ext))
    fig_cons.savefig("img/pibode_{}_masstot_.{}".format(shape, ext))
    fig_sepcke.savefig("img/pibode_{}_sepcke.{}".format(shape, ext))
    fig_gasdensity.savefig("img/pibode_{}_gasdensity.{}".format(shape, ext))
    fig_pibpaperwi.savefig("img/pibode_{}_wideac.{}".format(shape, ext))
    fig_pibpaperno.savefig("img/pibode_{}_nodeac.{}".format(shape, ext))


