#!/usr/bin/env python3


import numpy as np

from matplotlib import pyplot as plt
import matplotlib.ticker as mtick

from mevlib.scalar import sum_standard
from mevlib.scalar import cyl_intgtd_radial_terms, cyl_intgtd_axial_terms
from scipy.special import jn_zeros


plt.rc('font', size=10)
plt.rc('text', usetex=True)
plt.rc('axes', labelsize=10)
plt.rc('legend', fontsize=8)    # maybe should be 10, but 10 is too big

axsize = (2.75, 2.5)

def figsize(i, j):
    return (axsize[0] * j, axsize[1] * i)


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
        plotdata_single[k, i] = S(np.sqrt(phi2), R, H)
fig3, ax3 = plt.subplots(1, 1, figsize=figsize(1, 1))
lines = [
        ax3.semilogx(
            Hs / Rs, plotdata_single[k, :],
            'C{}-'.format(k+1), label=r"\( a^2 = {} \)".format(phi2)
            )
        for k, phi2 in enumerate(phi2s)
        ]
ax3.set_xlabel(r"\( H / (2 R) \)")
#ax3.legend()
annotposns = [(8.0, 0.31), (3.0, 0.44), (2.0, 0.64), (4.0, 0.95)]
for line, xy in zip(lines, annotposns):
    ax3.annotate(
        line[0].get_label(), xy=xy, fontsize=10,
        #color=line[0].get_color()
        color="black"
    )
ax3.yaxis.set_major_formatter(mtick.PercentFormatter(1.0))
ax3.grid()
ax3.set_ylabel("Effectiveness Factor")


plotdata_check = np.empty((3, len(Hs)))
for i, (R, H) in enumerate(zip(Rs, Hs)):
    plotdata_check[0, i] = S(phi2nom, R, H)
    plotdata_check[1, i] = AN(phi2nom, R, H)
    plotdata_check[2, i] = BK(phi2nom, R, H)
fig5, ax5 = plt.subplots(1, 1, figsize=figsize(1, 1))
ax5.semilogx(
    0.5 * Hs / Rs, plotdata_check[0, :], 'k-', label=r"\( A_N + B_K \)"
)
ax5.semilogx(0.5 * Hs / Rs, plotdata_check[1, :], 'b-', label=r"\( A_N \)")
ax5.semilogx(0.5 * Hs / Rs, plotdata_check[2, :], 'r-', label=r"\( B_K \)")
ax5.set_xlabel(r"\( H / (2 R) \)")
ax5.yaxis.set_major_formatter(mtick.PercentFormatter(1.0))
ax5.grid()
ax5.set_ylabel("Effectiveness Factor")
ax5.legend()


for i, fig in enumerate([fig3, fig5]):
    for ext in ['svg', 'pdf']:
        fig.tight_layout()
        fig.savefig("img/cyl_efffactor_singstage_{}.{}".format(
            "fig{}".format(i+1), ext
        ))


