#!/usr/bin/env python3


from matplotlib import pyplot as plt
import matplotlib.ticker as mtick

import numpy as np
from scipy.special import i0, i1, i0e, i1e

from nrel_cylinders import s1, s2


#plt.rc('font', size=16)
#plt.rc('text', usetex=True)
#plt.rc('axes', labelsize=20)


plt.rc('text', usetex=True)
plt.rc('axes', labelsize=12)


a = 00.8
R = 04.0
H = 15.0


ntruncs = list(range(2, 24, 1))
ktruncs = list(range(2, 24, 1))
nexact, kexact = 400, 800


s1exact, s2exact = s1(a, R, H, nexact), s2(a, R, H, kexact)


fig1, axs1 = plt.subplots(1, 2)
axs1[0].plot(
    ntruncs,
    [(s1(a, R, H, trunc) - s1exact) / s1exact for trunc in ntruncs],
    'k.'
)
axs1[1].plot(
    ktruncs,
    [(s2(a, R, H, trunc) - s2exact) / s2exact for trunc in ktruncs],
    'k.'
)
axs1[0].set_title(r"\( A_N \)")
axs1[1].set_title(r"\( B_K \)")
axs1[0].set_xlabel(r"\( N \)")
axs1[1].set_xlabel(r"\( K \)")
for i in range(2):
    axs1[i].yaxis.set_major_formatter(mtick.PercentFormatter(1.0))
    axs1[i].set_ylabel("Relative Difference")
    axs1[i].grid()


zs = np.linspace(10.0, 800.0, 201)
fig2, ax2 = plt.subplots(1, 1)
ct, cc = 300.0, 500.0
y0 = i1e(cc) / i0e(cc)
# tangent
m = 1.0 - y0 * (y0 + cc**(-1))
q = 2 * y0 / cc**2 - 2 * y0 - cc**(-1) + 2 * y0**3 + 3 * y0**2 / cc
# TODO secant
ax2.plot(
    zs, i1e(zs) / i0e(zs), 'k-', label=r"\( \frac{I_1(z)}{I_0(z)} \)"
)
ax2.plot(
    zs, np.where(zs < ct, i1e(zs) / i0e(zs), y0 + m * (zs - cc)), 'r--',
    label=r"Linear Approximation"
)
ax2.plot(
    zs, np.where(
        zs < ct, i1e(zs) / i0e(zs), y0 + m * (zs - cc) + q / 2 * (zs - cc)**2
    ),
    'g--', label=r"Quadratic Approximation"
)
ax2.legend()
ax2.set_xlabel(r"\( z \)")
ax2.grid()


if True:
    for i, fig in enumerate([fig1, fig2]):
        fig.tight_layout()
        for ext in ['svg', 'pdf']:
            fig.savefig("img/sumconv-{}.{}".format("fig{}".format(i+1), ext))
else:
    plt.show()

