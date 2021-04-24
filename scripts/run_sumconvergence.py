#!/usr/bin/env python3


from itertools import cycle
from matplotlib import pyplot as plt
import matplotlib.ticker as mtick

import numpy as np

from lib_scalar import sum_standard, sum_aitken
from lib_scalar import cyl_intgtd_radial_terms, cyl_intgtd_axial_terms
from lib_scalar import psm_intgtd_xdir


plt.rc('text', usetex=True)
plt.rc('axes', labelsize=12)


markers = cycle(['o', 'v', '>', '<', '^', 's', 'P', '*', 'X', 'd'])


a2 = 0.8**2
R, H = 80.0, 360.0
Lx, Ly, Lz = 120.0, 120.0, 360.0
rttruncs = np.array(list(range(2, 20, 1)))
nexact, kexact, pexact = 400, 800, 800


# cost units based on optimes.jl results
mult_cost = 1.0
divi_cost = 1.3
expn_cost = 2.5
root_cost = 1.5
sine_cost = 8.0
bes0_cost = 30.0
bes1_cost = 30.0


def cyl_ax_std_cost(N):
    term = bes0_cost + bes1_cost + 3 * divi_cost + 11 * mult_cost + root_cost
    return N * term

def cyl_rd_std_cost(K):
    term = 4 * expn_cost + 6 * mult_cost + 4 * divi_cost + root_cost
    return K * term

def cyl_ax_akn_cost(N):
    return cyl_ax_std_cost(N + 2) + (2 * N - 1) * (mult_cost + divi_cost)

def cyl_rd_akn_cost(K):
    return cyl_rd_std_cost(K + 2) + (2 * K - 1) * (mult_cost + divi_cost)

def psm_xd_cost(N):
    return 3 * mult_cost + N * (
        26 * mult_cost + 7 * divi_cost + 4 * expn_cost
        + root_cost + 2 * sine_cost
    )


cyl_radial_exact = sum_standard(cyl_intgtd_radial_terms(a2, R, H, nexact))
cyl_axial_exact = sum_standard(cyl_intgtd_axial_terms(a2, R, H, kexact))
psm_radial_exact = psm_intgtd_xdir(a2, Ly, Lx, Lz, pexact)
psm_axial_exact = psm_intgtd_xdir(a2, Lx, Ly, Lz, pexact)


methods = [     # name, valuefunction, exact, cost function
    (
        "Cylinder - Radial - Standard",
        lambda rtk: sum_standard(cyl_intgtd_radial_terms(a2, R, H, rtk**2)),
        cyl_radial_exact, cyl_rd_std_cost
    ),
    (
        "Cylinder - Axial - Standard",
        lambda rtk: sum_standard(cyl_intgtd_axial_terms(a2, R, H, rtk**2)),
        cyl_axial_exact, cyl_ax_std_cost
    ),
    (
        "Cylinder - Radial - Aitken",
        lambda rtk: sum_aitken(cyl_intgtd_radial_terms(a2, R, H, rtk**2)),
        cyl_radial_exact, cyl_rd_akn_cost
    ),
    (
        "Cylinder - Axial - Aitken",
        lambda rtk: sum_aitken(cyl_intgtd_axial_terms(a2, R, H, rtk**2)),
        cyl_axial_exact, cyl_ax_akn_cost
    ),
    (
        "Prism - Radial",       # one on cross-section
        lambda rtk: 2 * psm_intgtd_xdir(a2, Ly, Lx, Lz, rtk),
        2 * psm_radial_exact, lambda n: 2 * psm_xd_cost(n)
    ),
    (
        "Prism - Axial",        # one on ends
        lambda rtk: psm_intgtd_xdir(a2, Lx, Ly, Lz, rtk),
        psm_axial_exact, psm_xd_cost
    )
]


fig1, ax1 = plt.subplots(1, 1)
fig2, ax2 = plt.subplots(1, 1)
for i, (label, method, exact, cost) in enumerate(methods):
    mkr = next(markers)
    ax1.semilogy(
        rttruncs**2,
        [abs(method(rtk) - exact) / exact for rtk in rttruncs],
        "C{}{}".format(i+1, mkr), label=label
    )
    ax2.semilogy(
        cost(rttruncs**2),
        [abs(method(rtk) - exact) / exact for rtk in rttruncs],
        "C{}{}".format(i+1, mkr), label=label
    )
ax1.set_xlabel("Number of Terms")
ax1.set_ylabel("Relative Error")
ax1.legend()
ax1.yaxis.set_major_formatter(mtick.PercentFormatter(1.0, 4))
ax2.set_xlabel("Estimated Cost")
ax2.set_ylabel("Relative Error")
ax2.legend()
ax2.yaxis.set_major_formatter(mtick.PercentFormatter(1.0, 4))


if False:
    for n, fig in [('nterms', fig1), ('cost', fig2)]:
        fig.tight_layout()
        for ext in ['svg', 'pdf']:
            fig.savefig("img/sumconv-{}.{}", n, ext)
else:
    plt.show()


