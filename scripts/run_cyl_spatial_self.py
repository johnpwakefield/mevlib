#!/usr/bin/env python3


# there's a lot of screwy things in this script


import numpy as np
from matplotlib import pyplot as plt
from matplotlib import colors       #noqa W0611

from mevlib.shapes import Cylinder


plt.rc('font', size=10)
plt.rc('text', usetex=True)
plt.rc('axes', labelsize=10)
plt.rc('legend', fontsize=10)

axsize = (2.75, 2.5)

def figsize(i, j):
    return (axsize[0] * j, axsize[1] * i)


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
    plt.subplots(1, 1, figsize=figsize(1, 1)) for i in range(len(cases))
])
for i, (shp, phi2) in enumerate(cases):
    def solfunc(r, z, k):
        return (
            np.vectorize(shp.ptwise_axial)(phi2, k, r, z)
            +
            np.vectorize(shp.ptwise_radial)(phi2, k, r, z)
        )
    rmesh, zmesh = np.meshgrid(
        np.linspace(0.0, shp.R, Nr), np.linspace(0.0, shp.H, Nh)
    )
    refsoln = solfunc(rmesh, zmesh, kexact)
    assert(not np.any(np.isnan(refsoln)))
    hi1, hi2 = np.array(rmesh.shape) - 1
    l2errs[i] = [
        np.sqrt(np.sum((solfunc(rmesh, zmesh, k) - refsoln).flatten()**2))
        * np.sqrt(shp.R / Nr * shp.H / Nh)
        for k in Ks
    ]
    wterrs[i] = [
        np.sqrt(np.sum((
            4 * np.pi * rmesh**2 * shp.R / Nr * shp.H / Nh
            * (solfunc(rmesh, zmesh, k) - refsoln)**2
        ).flatten()))
        for k in Ks
    ]
    cm = singaxs[i].contourf(rmesh, zmesh, solfunc(rmesh, zmesh, kspatial))
    singaxs[i].set_xlabel(r"\( r \)")
    singaxs[i].set_ylabel(r"\( z \)")
    cbar = plt.colorbar(cm, ax=singaxs[i])
    cbar.set_label(r"\( \hat{Y} \)")
for fig in singfigs:
    fig.tight_layout()

l2fig, l2ax = plt.subplots(1, 1, figsize=figsize(2, 2))
for i, (errs, mkr) in enumerate(zip(l2errs, ['x', 'o', '+', '<', '>'])):
    l2ax.semilogy(
        Ks, errs, 'C{}'.format(i+1) + mkr, label="Case {}".format(i + 1)
    )
wtfig, wtax = plt.subplots(1, 1, figsize=figsize(2, 2))
for i, (errs, mkr) in enumerate(zip(wterrs, ['x', 'o', '+', '<', '>'])):
    wtax.semilogy(
        Ks, errs, 'C{}'.format(i+1) + mkr, label="Case {}".format(i + 1)
    )
for ax in [l2ax, wtax]:
    ax.grid()
    ax.set_xlabel(r"Number of terms \( K \)")
    ax.set_ylabel(r"\( L^2 \) error")
    ax.legend()


for ext in ['svg', 'pdf']:
    l2fig.savefig("img/cyl_spatial_self_l2errs.{}".format(ext))
    wtfig.savefig("img/cyl_spatial_self_wterrs.{}".format(ext))
    for i, fig in enumerate(singfigs):
        fig.savefig("img/cyl_spatial_self_case{}.{}".format(i+1, ext))


