#!/usr/bin/env python3


import pickle

from math import sqrt

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import colors

from nrel_cylinders import u
from scipy.special import jn_zeros


plt.rc('font', size=16)
plt.rc('text', usetex=True)
plt.rc('axes', labelsize=20)


nexact, kexact = 128, 64
zc = jn_zeros(0, kexact)


D = 5.4423
cases = [
    (240.0, 240.0, 60.0, 1.0),
    (140.0, 840.0, 60.0, 1.0),
    (180.0, 360.0, 60.0, 4.0),
    (180.0, 360.0, 60.0, 0.1),
    (180.0, 360.0, 60.0, 1.0)
]
fns = [
    (
        "../43.08-External_Verification/casestudies_030321/Case{}_A.pickle"
    ).format(i+1)
    for i in range(5)
]
rates = [1.5117e-3, 1.5117e-3, 6.0470e-3, 1.5117e-4, 1.5117e-3]
compphi2s = [rate * L**2 / D for (_, _, L, _), rate in zip(cases, rates)]
cases = [(*case, compphi2) for case, compphi2 in zip(cases, compphi2s)]


for i, ((R, H, L, phi2, cphi2), fn) in enumerate(zip(cases, fns)):
    fig, axs = plt.subplots(1, 3, figsize=(12.0, 6.0))
    cms = [None, None, None]
    def solfunc(r, z):
        return u(
            sqrt(phi2), R / L, H / L, nexact, kexact,
            r / L, z / L,
            zero_cache=zc
        )
    with open(fn, 'rb') as fh:
        refsoln = pickle.load(fh)
    oursoln = solfunc(refsoln['rmesh'], refsoln['zmesh'])
    versoln = refsoln['umesh']
    cms[1] = axs[1].contourf(
        refsoln['rmesh'], refsoln['zmesh'], versoln
    )
    cms[0] = axs[0].contourf(
        refsoln['rmesh'], refsoln['zmesh'], oursoln,
        levels=cms[1].levels
    )
    err = np.abs(oursoln - versoln) / versoln
    print("Case {} error range ({}, {})".format(i+1, err.min(), err.max()))
    lvls = np.power(10.0, np.arange(
        np.floor(np.log10(err.min()) - 1), np.ceil(np.log10(err.max()) + 1)
    ))
    cms[2] = axs[2].contourf(
        refsoln['rmesh'], refsoln['zmesh'], err, lvls, norm=colors.LogNorm()
    )
    axs[0].set_title("Our Solution")
    axs[1].set_title("Reference Solution")
    axs[2].set_title("Relative Error")
    cbs = [plt.colorbar(cm, ax=ax) for ax, cm in zip(axs, cms)]
    for ax in axs:
        ax.set_xlabel("r")
        ax.set_ylabel("z")
#   cbs[2].ax.set_yticklabels([
#       "{:.2%}".format(x)[:-1] + r"\%" for x in cbs[2].get_ticks()
#   ])
    fig.tight_layout()
    for ext in ['svg', 'pdf']:
        fig.savefig("img/comparison_case{}.".format(i+1) + ext)


