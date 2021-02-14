#!/usr/bin/env python3


import numpy as np

from matplotlib import pyplot as plt
import matplotlib.ticker as mtick

from nrel_cylinders import s, s1, s2
from scipy.special import jn_zeros


plt.rc('text', usetex=True)
plt.rc('axes', labelsize=12)


B = - np.array([
    [-1.36211  , -0.0      , -0.0       , -0.0, -0.0, -0.0],
    [+0.259634 , -0.0751174, -0.0       , -0.0, -0.0, -0.0],
    [+0.796909 , +0.030285 , -0.0244072 , -0.0, -0.0, -0.0],
    [+0.213697 , +0.0212921, +0.0119698 , -0.0, -0.0, -0.0],
    [+0.0209471, +0.0054222, +0.00280542, -0.0, -0.0, -0.0],
    [+0.0709262, +0.0181181, +0.00963195, -0.0, -0.0, -0.0]
])

EVs = np.array([
    [+0.840825 , +0.0      , +0.0       , +0.0, +0.0, +0.0],
    [-0.169625 , +0.838677 , +0.0       , +0.0, +0.0, +0.0],
    [-0.497063 , -0.500873 , +0.842312  , +0.0, +0.0, +0.0],
    [-0.124895 , -0.157911 , -0.413088  , +0.0, +0.0, +1.0],
    [-0.0112315, -0.0418322, -0.0968175 , +0.0, +1.0, +0.0],
    [-0.0380112, -0.138062 , -0.332407  , +1.0, +0.0, +0.0]
])

EVsinv = np.linalg.inv(EVs)

Lams = np.diagonal(B)


nexact, kexact = 12, 18
zc = jn_zeros(0, max(nexact, kexact))

bdrycond = np.array([1.0, 0.8, 0.6, 0.4, 0.2, 0.2])


def eff_factor(R, H, ntrunc, ktrunc, Cbdry, zero_cache=None):
    mult = np.array([
        1.0 if lam == 0.0 else
        s(np.sqrt(lam), R, H, ntrunc, ktrunc, zero_cache=zero_cache)
        for lam in Lams
    ])
    res = np.dot(EVs, mult * np.dot(EVsinv, Cbdry)) / Cbdry
    if np.any(np.logical_not(np.isfinite(res))):
        print("Non finite values in result.")
        print(res)
    return res


Npts = 201
ratios = 10**np.linspace(-2.0, 4.0, Npts)
V = np.pi * 2.0**8
Rs = (V / (np.pi * ratios))**(1.0/3)
Hs = ratios * Rs


plotdata_single = np.empty((1, len(Hs)))
for i, (R, H) in enumerate(zip(Rs, Hs)):
    plotdata_single[:,i] = s(
        np.sqrt(Lams[0]),
        R, H, nexact, kexact, zero_cache=zc
    )
fig1, axs1 = plt.subplots(1, 1)
axs1.loglog(Hs / Rs, plotdata_single[0,:], 'k-')
axs1.set_xlabel(r"\( H / R \)")
axs1.yaxis.set_major_formatter(mtick.PercentFormatter(1.0))
axs1.grid()
axs1.set_ylabel("Effectiveness Factor")

plotdata_sepsums_1 = np.empty((1, len(Hs)))
plotdata_sepsums_2 = np.empty((1, len(Hs)))
for i, (R, H) in enumerate(zip(Rs, Hs)):
    plotdata_sepsums_1[:,i] = s1(np.sqrt(Lams[0]), R, H, nexact)
    plotdata_sepsums_2[:,i] = s2(np.sqrt(Lams[0]), R, H, kexact, zero_cache=zc)
fig1, axs1 = plt.subplots(1, 1)
axs1.loglog(Hs / Rs, plotdata_sepsums_1[0,:], 'r-', label="s1")
axs1.loglog(Hs / Rs, plotdata_sepsums_2[0,:], 'b-', label="s2")
axs1.set_xlabel(r"\( H / R \)")
axs1.legend()
axs1.grid()

plotdata_multistage = np.empty((6, len(Hs)))
for i, (R, H) in enumerate(zip(Rs, Hs)):
    plotdata_multistage[:,i] = eff_factor(
        R, H, nexact, kexact, bdrycond, zero_cache=zc
    )
fig2, axs2 = plt.subplots(2, 3)
axs2 = axs2.flatten()
for i, ax in enumerate(axs2):
    component = r"\( \eta_{} \)".format(i+1)
    ax.loglog(
        Hs / Rs, plotdata_multistage[i,:],
        'C{}-'.format(i), label=component
    )
    ax.set_title(component)
    ax.set_xlabel(r"\( H / R \)")
    ax.yaxis.set_major_formatter(mtick.PercentFormatter(1.0))
    ax.set_ylabel("Effectiveness Factor")
    ax.grid()


if True:
    for i, fig in enumerate([fig1, fig2]):
        fig.tight_layout()
        for ext in ['svg', 'pdf']:
            fig.savefig("img/volume-{}.{}".format("fig{}".format(i+1), ext))
else:
    for fig in [fig1, fig2]:
        fig.tight_layout()
    plt.show()


