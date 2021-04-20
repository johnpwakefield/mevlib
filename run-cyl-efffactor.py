#!/usr/bin/env python3


import numpy as np

from matplotlib import pyplot as plt
import matplotlib.ticker as mtick

from nrel_cylinders import s, s1, s2
from scipy.special import jn_zeros


plt.rc('font', size=16)
plt.rc('text', usetex=True)
plt.rc('axes', labelsize=20)


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


nexact, kexact = 128, 128
zc = jn_zeros(0, max(nexact, kexact))

bdrycond = np.array([1.0, 0.7, 0.7, 0.1, 0.1, 0.1])


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

phi2s = [10.0, 1.0, 0.1, 0.0]
plotdata_single = np.empty((len(phi2s), len(Hs)))
for k, phi2 in enumerate(phi2s):
    for i, (R, H) in enumerate(zip(Rs, Hs)):
        plotdata_single[k,i] = s(
            np.sqrt(phi2),
            R, H, nexact, kexact, zero_cache=zc
        )
fig3, axs3 = plt.subplots(1, 1)
for k, phi2 in enumerate(phi2s):
    axs3.semilogx(
        Hs / Rs, plotdata_single[k,:],
        'C{}-'.format(k+1), label=r"\( \phi^2 = {} \)".format(phi2)
    )
axs3.set_xlabel(r"\( H / R \)")
axs3.legend()
axs3.yaxis.set_major_formatter(mtick.PercentFormatter(1.0))
axs3.grid()
axs3.set_ylabel("Effectiveness Factor")

plotdata_check = np.empty((1, len(Hs)))
for i, (R, H) in enumerate(zip(Rs, Hs)):
    plotdata_check[:,i] = s(0.0, R, H, nexact, kexact, zero_cache=zc)
fig5, axs5 = plt.subplots(1, 1)
axs5.semilogx(Hs / Rs, plotdata_check[0,:], 'k-')
axs5.set_xlabel(r"\( H / R \)")
axs5.yaxis.set_major_formatter(mtick.PercentFormatter(1.0))
axs5.grid()
axs5.set_ylabel("Effectiveness Factor")

plotdata_sepsums_1 = np.empty((1, len(Hs)))
plotdata_sepsums_2 = np.empty((1, len(Hs)))
for i, (R, H) in enumerate(zip(Rs, Hs)):
    plotdata_sepsums_1[:,i] = s1(np.sqrt(Lams[0]), R, H, nexact)
    plotdata_sepsums_2[:,i] = s2(np.sqrt(Lams[0]), R, H, kexact, zero_cache=zc)
fig1, axs1 = plt.subplots(1, 1)
axs1.semilogx(Hs / Rs, plotdata_sepsums_1[0,:], 'r-', label=r"\( A_N \)")
axs1.semilogx(Hs / Rs, plotdata_sepsums_2[0,:], 'b-', label=r"\( B_K \)")
axs1.set_xlabel(r"\( H / R \)")
axs1.legend()
axs1.grid()

plotdata_multistage = np.empty((6, len(Hs)))
for i, (R, H) in enumerate(zip(Rs, Hs)):
    plotdata_multistage[:,i] = eff_factor(
        R, H, nexact, kexact, bdrycond, zero_cache=zc
    )
fig2, axs2 = plt.subplots(2, 3, figsize=(8.0, 6.0))
axs2 = axs2.flatten()
for i, ax in enumerate(axs2):
    component = r"\( \eta_{} \)".format(i+1)
    ax.semilogx(
        Hs / Rs, plotdata_multistage[i,:],
        'C{}-'.format(i), label=component
    )
    ax.set_title(component)
    ax.set_xlabel(r"\( H / R \)")
    ax.yaxis.set_major_formatter(mtick.PercentFormatter(1.0))
    ax.set_ylabel("Effectiveness Factor")
    ax.grid()

maxterms = 32
gammas = [0.05, 1.0, 20.0]
Rs = [2.0, 40.0]
phi = 20.0
plotdata_numterms = np.empty((len(Rs),maxterms,len(gammas)))
plotdata_numtermsexact = np.empty((len(gammas),))
for i, gamma in enumerate(gammas):
    for j in range(maxterms):
        for k, R in enumerate(Rs):
            plotdata_numterms[k,j,i] = s(
                phi / R, R, 2.0 * R * gamma, j+1, j+1, zero_cache=zc
            )
    plotdata_numtermsexact[i] = s(
        phi, 1.0, 2.0 * gamma, nexact, kexact, zero_cache=zc
    )
fig4, axs4 = plt.subplots(len(Rs), len(gammas))
for k, R in enumerate(Rs):
    for i, gamma in enumerate(gammas):
        axs4[k,i].semilogy(
            np.arange(maxterms) + 1, plotdata_numterms[k,:,i], 'k.'
        )
        axs4[k,i].set_title(r"\( \gamma = {} \)".format(gamma))
        axs4[k,i].set_xlabel(r"Number of Terms")
        axs4[k,i].set_ylabel(r"Effectiveness Factor")
        axs4[k,i].set_ylim((0.1, 1.0))
fig5idx = 1
fig5, ax5 = plt.subplots(1, 1)
ax5.plot(np.arange(maxterms) + 1, plotdata_numterms[0,:,fig5idx], 'k.')
ax5.plot(
    [1, maxterms], plotdata_numtermsexact[fig5idx] * np.ones((2,)), 'r:'
)
ax5.set_title(r"\( \gamma = {} \)".format(gammas[fig5idx]))
ax5.set_xlabel(r"Number of Terms")
ax5.set_ylabel(r"Effectiveness Factor")
#ax5.set_ylim((0.1, 1.0))
ax5.grid()


if False:
    for i, fig in enumerate([fig1, fig2, fig3, fig4, fig5]):
        fig.tight_layout()
        for ext in ['svg', 'pdf']:
            fig.savefig("img/efffactor-{}.{}".format("fig{}".format(i+1), ext))
else:
    for fig in [fig1, fig2, fig3]:
        fig.tight_layout()
    plt.show()


