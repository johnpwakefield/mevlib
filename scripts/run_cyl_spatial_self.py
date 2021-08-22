#!/usr/bin/env python3


import pkgutil
import pickle

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import colors

from mevlib.shapes import Cylinder
from mevlib.options import imgpath, showfigs


plt.rc('font', size=12)
plt.rc('text', usetex=True)
plt.rc('axes', labelsize=18)


kexact = 200
kspatial = 16
Nr, Nh = 198, 199


D = 5.4423
cases = [
    (240.0, 240.0, 1.0),
    (140.0, 840.0, 1.0),
    (180.0, 360.0, 4.0),
    (180.0, 360.0, 0.1),
    (180.0, 360.0, 1.0)
]


rates = [1.5117e-3, 1.5117e-3, 6.0470e-3, 1.5117e-4, 1.5117e-3]
compphi2s = [rate / D for rate in rates]
cases = [
    (Cylinder(case[0], case[1]), compphi2)
    for case, compphi2 in zip(cases, compphi2s)
]
Ks = np.arange(1, 42, 4)
l2errs = [None for c in cases]  # errors in L2 norm
wterrs = [None for c in cases]  # errors in weighted norm


singfigs, singaxs = zip(*[
    plt.subplots(1, 1, figsize=(4.0, 4.0)) for i in range(len(cases))
])
for i, (shp, phi2) in enumerate(cases):
    def solfunc(r, z, k):
        return np.vectorize(shp.ptwise_axial)(phi2, k, r, z)
    rmesh, zmesh = np.meshgrid(
        np.linspace(0.0, shp.R, Nr), np.linspace(0.0, shp.H, Nh)
    )
    refsoln = solfunc(rmesh, zmesh, kexact)
    assert(not np.any(np.isnan(refsoln)))
    hi1, hi2 = np.array(rmesh.shape) - 1
    l2errs[i] = [
        np.sqrt(np.sum((solfunc(rmesh, zmesh, k) - refsoln).flatten()**2))
        for k in Ks
    ]
    wterrs[i] = [
        np.sqrt(np.sum((
            2 * np.pi * rmesh * (solfunc(rmesh, zmesh, k) - refsoln)**2
        ).flatten()))
        for k in Ks
    ]
    cm = singaxs[i].contourf(rmesh, zmesh, solfunc(rmesh, zmesh, kspatial))
    singaxs[i].set_xlabel(r"\( r \)")
    singaxs[i].set_ylabel(r"\( z \)")
    plt.colorbar(cm, ax=singaxs[i])
for fig in singfigs:
    fig.tight_layout()


l2fig, l2ax = plt.subplots(1, 1, figsize=(6.0, 6.0))
for i, (errs, mkr) in enumerate(zip(l2errs, ['x', 'o', '+', '<', '>'])):
    l2ax.semilogy(
        Ks, errs, 'C{}'.format(i+1) + mkr, label="Case {}".format(i + 1)
    )
wtfig, wtax = plt.subplots(1, 1, figsize=(6.0, 6.0))
for i, (errs, mkr) in enumerate(zip(wterrs, ['x', 'o', '+', '<', '>'])):
    wtax.semilogy(
        Ks, errs, 'C{}'.format(i+1) + mkr, label="Case {}".format(i + 1)
    )
for ax in [l2ax, wtax]:
    ax.grid()
    ax.set_xlabel(r"Number of terms \( K \)")
    ax.set_ylabel(r"\( L^2 \) error")
    ax.legend()


if showfigs():
    plt.show()
else:
    for ext in ['svg', 'pdf']:
        l2fig.savefig(imgpath("spatial_self_l2errs.{}".format(ext)))
        wtfig.savefig(imgpath("spatial_self_wterrs.{}".format(ext)))
        for i, fig in enumerate(singfigs):
            fig.savefig(imgpath("spatial_self_case{}.{}".format(i+1, ext)))


