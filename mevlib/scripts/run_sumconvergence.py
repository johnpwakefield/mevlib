#!/usr/bin/env python3


from itertools import cycle
from matplotlib import pyplot as plt
import matplotlib.ticker as mtick

import numpy as np

from mevlib.scalar import sum_standard, sum_aitken
from mevlib.scalar import cyl_intgtd_radial_terms, cyl_intgtd_axial_terms
from mevlib.scalar import psm_intgtd_xdir, psm_intgtd_spec
from mevlib.scalar import psm_intgtd_diff
from mevlib.options import imgpath, showfigs


plt.rc('text', usetex=True)
plt.rc('axes', labelsize=12)


markers = cycle(['v', '>', '<', '^', 's', 'd', '*', 'X', '+'])


a2 = 0.8**2
R, H = 80.0, 360.0
Lx, Ly, Lz = 120.0, 120.0, 360.0
nexact, kexact, pexact = 400, 800, 800

cyltrunc = np.array(list(range(3, 24)))
psmseriestrunc = np.array(list(range(3, 8)))
psmspectrunc = np.array(list(range(4, 8)))
psmdifftrunc = np.array(list(range(4, 8)))


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
    return 3 * mult_cost + N**2 * (
        26 * mult_cost + 7 * divi_cost + 4 * expn_cost
        + root_cost + 2 * sine_cost
    )

def psm_spec_cost(N):
    #TODO
    return N

def psm_diff_cost(N):
    #TODO
    return N


cyl_radial_exact = sum_standard(cyl_intgtd_radial_terms(a2, R, H, nexact))
cyl_axial_exact = sum_standard(cyl_intgtd_axial_terms(a2, R, H, kexact))
psm_radial_exact = psm_intgtd_xdir(a2, Ly, Lx, Lz, pexact)
psm_axial_exact = psm_intgtd_xdir(a2, Lx, Ly, Lz, pexact)


methods = [
    # name, valuefunction, exact, termadj, cost function, number of terms
    (
        "Cylinder - Radial - Standard",
        lambda tk: sum_standard(cyl_intgtd_radial_terms(a2, R, H, tk)),
        cyl_radial_exact, lambda tk: tk, cyl_rd_std_cost, cyltrunc
    ),
    (
        "Cylinder - Radial - Aitken",
        lambda tk: sum_aitken(cyl_intgtd_radial_terms(a2, R, H, tk)),
        cyl_radial_exact, lambda tk: tk, cyl_rd_akn_cost, cyltrunc
    ),
    (
        "Cylinder - Axial - Standard",
        lambda tk: sum_standard(cyl_intgtd_axial_terms(a2, R, H, tk)),
        cyl_axial_exact, lambda tk: tk, cyl_ax_std_cost, cyltrunc
    ),
    (
        "Cylinder - Axial - Aitken",
        lambda tk: sum_aitken(cyl_intgtd_axial_terms(a2, R, H, tk)),
        cyl_axial_exact, lambda tk: tk, cyl_ax_akn_cost, cyltrunc
    ),
    (
        "Prism - Radial - Series",      # one on cross-section
        lambda rtk: 2 * psm_intgtd_xdir(a2, Ly, Lx, Lz, rtk),
        2 * psm_radial_exact, lambda rtk: rtk**2, lambda n: 2 * psm_xd_cost(n),
        psmseriestrunc
    ),
    (
        "Prism - Axial - Series",       # one on ends
        lambda rtk: psm_intgtd_xdir(a2, Lx, Ly, Lz, rtk),
        psm_axial_exact, lambda rtk: rtk**2, psm_xd_cost,
        psmseriestrunc
    ),
#   (
#       "Prism - Total - Series",
#       lambda rtk: (
#           2 * psm_intgtd_xdir(a2, Ly, Lx, Lz, rtk)
#           + psm_intgtd_xdir(a2, Lx, Ly, Lz, rtk)
#       ),
#       2 * psm_radial_exact + psm_axial_exact,
#       lambda rtk: rtk**2, psm_xd_cost,
#       psmseriestrunc
#   ),
#   (
#       "Prism - Pseudospectral",
#       lambda tk: psm_intgtd_spec(a2, Lx, Ly, Lz, tk),
#       2 * psm_radial_exact + psm_axial_exact,
#       lambda tk: tk**2, lambda n: psm_spec_cost(n), psmspectrunc
#   ),
#   (
#       "Prism - Finite Difference",
#       lambda tk: psm_intgtd_diff(a2, Lx, Ly, Lz, tk),
#       2 * psm_radial_exact + psm_axial_exact,
#       lambda tk: tk**2, lambda n: psm_diff_cost(n), psmdifftrunc
#   )
]


fig1, ax1 = plt.subplots(1, 1)
fig2, ax2 = plt.subplots(1, 1)
fig3, ax3 = plt.subplots(1, 1)
for i, (label, method, exact, termadj, cost, truncs) in enumerate(methods):
    mkr = next(markers)
    ax1.semilogy(
        termadj(truncs),
        [abs(method(tk) - exact) / exact for tk in truncs],
        "C{}{}".format((i % 9) + 1, mkr), label=label
    )
    ax2.semilogy(
        cost(truncs),
        [abs(method(tk) - exact) / exact for tk in truncs],
        "C{}{}".format((i % 9) + 1, mkr), label=label
    )
    ax3.semilogy(
        termadj(truncs), [method(tk) for tk in truncs],
        "C{}{}".format((i % 9) + 1, mkr), label=label
    )
    print(label)
    print([method(tk) for tk in truncs])
    print(exact)
ax1.grid()
ax1.set_xlabel("Number of Terms")
ax1.set_ylabel("Relative Error")
ax1.legend()
ax1.yaxis.set_major_formatter(mtick.PercentFormatter(1.0, 4))
ax2.grid()
ax2.set_xlabel("Estimated Cost")
ax2.set_ylabel("Relative Error")
ax2.legend()
ax2.yaxis.set_major_formatter(mtick.PercentFormatter(1.0, 4))
ax3.grid()
ax3.set_xlabel("Number of Terms")
ax3.set_ylabel("Effectiveness Factor")
ax3.legend()


if not showfigs():
    for n, fig in [('nterms', fig1), ('cost', fig2), ('efactor', fig3)]:
        fig.tight_layout()
        for ext in ['svg', 'pdf']:
            fig.savefig(imgpath("sumconv-{}.{}".format(n, ext)))
else:
    plt.show()


