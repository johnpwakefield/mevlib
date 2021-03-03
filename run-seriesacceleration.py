#!/usr/bin/env python3


from matplotlib import pyplot as plt
import matplotlib.ticker as mtick

import numpy as np
from scipy.special import i0, i1, i0e, i1e, jn_zeros

from nrel_cylinders import s1, s2, gn, ln


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


# cost units based on optimes results
mult_cost = 1.0
divi_cost = 1.3
expn_cost = 2.5
root_cost = 1.5
bes0_cost = 30.0
bes1_cost = 30.0

def s1og_cost(N):
    term = bes0_cost + bes1_cost + 3 * divi_cost + 11 * mult_cost + root_cost
    return N * term

def s2og_cost(K):
    term = 4 * expn_cost + 6 * mult_cost + 4 * divi_cost + root_cost
    return K * term

def s1ak_cost(N):
    return s1og_cost(N + 2) + (2 * N - 1) * (mult_cost + divi_cost)

def s2ak_cost(K):
    return s2og_cost(K + 2) + (2 * K - 1) * (mult_cost + divi_cost)

def s1ct_cost(N, c):
    Ncrit = int(0.5 * (H / (np.pi * R) * np.sqrt(c**2 - a**2 * R**2) - 1.0))
    term = 3 * divi_cost + 13 * mult_cost + root_cost
    return s1og_cost(min(N, Ncrit)) + max(N - Ncrit, 0) * term



s1exact, s2exact = s1(a, R, H, nexact), s2(a, R, H, kexact)

def ANterm(n):
    g = gn(a, H, n)
    return (
        16.0 / (np.pi**2 * R * g * (2 * n + 1)**2)
        * i1e(g * R) / i0e(g * R)
    )

def BKterm(k, alpha):
    l = ln(a, R, alpha)
    return 8.0 * np.tanh(l * H / 2) / (H * l * alpha**2)

def aitkenAN(a, R, H, trunc):
    term = ANterm
    anp1, anp2, anp3 = term(0), term(1), term(2)
    cum = anp1 - anp2**2 / (anp3 - anp2)
    for n in range(trunc-1):
        anp1, anp2, anp3 = anp2, anp3, term(n+3)
        cum += anp1 - anp2**2 / (anp3 - anp2) + anp1**2 / (anp2 - anp1)
    return cum

def aitkenBK(a, R, H, trunc):
    zeros = jn_zeros(0, trunc+3)
    term = lambda k: BKterm(k, zeros[k])
    anp1, anp2, anp3 = term(0), term(1), term(2)
    cum = anp1 - anp2**2 / (anp3 - anp2)
    for n in range(trunc-1):
        anp1, anp2, anp3 = anp2, anp3, term(n+3)
        cum += anp1 - anp2**2 / (anp3 - anp2) + anp1**2 / (anp2 - anp1)
    return cum


fig1, axs1 = plt.subplots(2, 2)
axs1[0,0].semilogy(
    ntruncs,
    [abs((s1(a, R, H, trunc) - s1exact) / s1exact) for trunc in ntruncs],
    'k.', label="Standard"
)
axs1[0,0].semilogy(
    ntruncs,
    [
        abs((aitkenAN(a, R, H, trunc) - s1exact) / s1exact)
        for trunc in ntruncs
    ],
    'gx', label=r"Aitken \( \delta^2 \)"
)
axs1[0,1].semilogy(
    ktruncs,
    [abs((s2(a, R, H, trunc) - s2exact) / s2exact) for trunc in ktruncs],
    'k.'
)
axs1[0,1].semilogy(
    ktruncs,
    [abs((aitkenBK(a, R, H, trunc) - s2exact) / s2exact) for trunc in ktruncs],
    'gx'
)
axs1[0,0].set_title(r"\( A_N \)")
axs1[0,1].set_title(r"\( B_K \)")
axs1[0,0].set_xlabel(r"\( N \)")
axs1[0,1].set_xlabel(r"\( K \)")
axs1[0,0].legend()
for i in range(2):
    axs1[0,i].yaxis.set_major_formatter(mtick.PercentFormatter(1.0))
    axs1[0,i].set_ylabel("Relative Difference")
    axs1[0,i].grid()


axs1[1,0].semilogy(
    [s1og_cost(trunc) for trunc in ntruncs],
    [abs((s1(a, R, H, trunc) - s1exact) / s1exact) for trunc in ntruncs],
    'k.', label="Standard"
)
axs1[1,0].semilogy(
    [s1og_cost(trunc) for trunc in ntruncs],
    [abs((s1(1.0, a * R, a * H, trunc) - s1exact) / s1exact) for trunc in ntruncs],
    'k.', label="Change of Variables"
)
axs1[1,0].semilogy(
    [s1ak_cost(trunc) for trunc in ntruncs],
    [
        abs((aitkenAN(a, R, H, trunc) - s1exact) / s1exact)
        for trunc in ntruncs
    ],
    'gx', label=r"Aitken \( \delta^2 \)"
)
axs1[1,1].semilogy(
    [s2og_cost(trunc) for trunc in ktruncs],
    [abs((s2(a, R, H, trunc) - s2exact) / s2exact) for trunc in ktruncs],
    'k.'
)
axs1[1,1].semilogy(
    [s2ak_cost(trunc) for trunc in ktruncs],
    [abs((aitkenBK(a, R, H, trunc) - s2exact) / s2exact) for trunc in ktruncs],
    'gx'
)
axs1[1,0].set_title(r"\( A_N \)")
axs1[1,1].set_title(r"\( B_K \)")
axs1[1,0].set_xlabel(r"Cost Units")
axs1[1,1].set_xlabel(r"Cost Units")
axs1[1,0].legend()
for i in range(2):
    axs1[1,i].yaxis.set_major_formatter(mtick.PercentFormatter(1.0))
    axs1[1,i].set_ylabel("Relative Difference")
    axs1[1,i].grid()


if True:
    for i, fig in enumerate([fig1]):
        fig.tight_layout()
        for ext in ['svg', 'pdf']:
            fig.savefig(
                "img/seriesaccel-{}.{}".format("fig{}".format(i+1), ext)
            )
else:
    plt.show()


