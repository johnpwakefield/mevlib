#!/usr/bin/env python3


import numpy as np

from matplotlib import pyplot as plt
import matplotlib.ticker as mtick

from mevlib.scalar import sum_standard
from mevlib.scalar import cyl_intgtd_radial_terms, cyl_intgtd_axial_terms
from scipy.special import jn_zeros

from mevlib.options import showfigs, imgpath


plt.rc('font', size=16)
plt.rc('text', usetex=True)
plt.rc('axes', labelsize=20)


phi2s = [10.0, 1.0, 0.1, 0.0]
phi2nom = 1.36


nexact, kexact = 128, 128
zc = jn_zeros(0, max(nexact, kexact))


def AN(phi2, R, H):
    return sum_standard(cyl_intgtd_radial_terms(phi2, R, H, nexact))

def BK(phi2, R, H):
    return sum_standard(cyl_intgtd_axial_terms(
        phi2, R, H, kexact, zero_cache=zc
    ))

def S(phi2, R, H):
    return AN(phi2, R, H) + BK(phi2, R, H)


Npts = 201
ratios = 10**np.linspace(-2.0, 4.0, Npts)
V = np.pi * 2.0**8
Rs = (V / (np.pi * ratios))**(1.0/3)
Hs = ratios * Rs


plotdata_single = np.empty((len(phi2s), len(Hs)))
for k, phi2 in enumerate(phi2s):
    for i, (R, H) in enumerate(zip(Rs, Hs)):
        plotdata_single[k,i] = S(np.sqrt(phi2), R, H)
fig3, axs3 = plt.subplots(1, 1)
for k, phi2 in enumerate(phi2s):
    axs3.semilogx(
        Hs / Rs, plotdata_single[k,:],
        'C{}-'.format(k+1), label=r"\( \phi^2 = {} \)".format(phi2)
    )
axs3.set_xlabel(r"\( H / (2 R) \)")
axs3.legend()
axs3.yaxis.set_major_formatter(mtick.PercentFormatter(1.0))
axs3.grid()
axs3.set_ylabel("Effectiveness Factor")


plotdata_check = np.empty((3, len(Hs)))
for i, (R, H) in enumerate(zip(Rs, Hs)):
    plotdata_check[0,i] = S(phi2nom, R, H)
    plotdata_check[1,i] = AN(phi2nom, R, H)
    plotdata_check[2,i] = BK(phi2nom, R, H)
fig5, axs5 = plt.subplots(1, 1)
axs5.semilogx(
    0.5 * Hs / Rs, plotdata_check[0,:], 'k-', label=r"\( A_N + B_K \)"
)
axs5.semilogx(0.5 * Hs / Rs, plotdata_check[1,:], 'b-', label=r"\( A_N \)")
axs5.semilogx(0.5 * Hs / Rs, plotdata_check[2,:], 'r-', label=r"\( B_K \)")
axs5.set_xlabel(r"\( H / (2 R) \)")
axs5.yaxis.set_major_formatter(mtick.PercentFormatter(1.0))
axs5.grid()
axs5.set_ylabel("Effectiveness Factor")
axs5.legend()


if showfigs():
    plt.show()
else:
    for ext in ['svg', 'pdf']:
        for i, fig in enumerate([fig3, fig5]):
            fig.tight_layout()
            fig.savefig(imgpath("cyl_singstage_efffactor_{}.{}".format(
                "fig{}".format(i+1), ext
            )))

