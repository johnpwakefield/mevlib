#!/usr/bin/env python3


from math import sqrt

import pkgutil
import pickle

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import colors

from mevlib.shapes import Cylinder


plt.rc('font', size=10)
plt.rc('text', usetex=True)
plt.rc('axes', labelsize=10)
plt.rc('legend', fontsize=10)

axsize = (2.75, 2.5)
arscale = 0.6

def figsize(i, j):
    return (axsize[0] * j, axsize[1] * i)


nexact, kexact = 64, 32
nspatial, kspatial = 32, 16


D = 5.4423
cases = [
    (240.0, 240.0, 1.0),
    (140.0, 840.0, 1.0),
    (180.0, 360.0, 4.0),
    (180.0, 360.0, 0.1),
    (180.0, 360.0, 1.0)
]
fns = [
    (
        "/data/singlestage_cyl_case{}.pickle"
    ).format(i+1)
    for i in range(len(cases))
]
rates = [1.5117e-3, 1.5117e-3, 6.0470e-3, 1.5117e-4, 1.5117e-3]
compphi2s = [rate / D for rate in rates]
pcases = [
    (Cylinder(case[0], case[1]), compphi2)
    for case, compphi2 in zip(cases, compphi2s)
]
Ns = np.arange(1, 42, 4)
l2errs = [None for c in cases]  # errors in L2 norm
wterrs = [None for c in cases]  # errors in weighted norm


multfig, multaxs = plt.subplots(5, 2, figsize=figsize(5, 2))
# note the individual ones are made at roughly the right aspect ratio
singfigs, singaxs = zip(*[
    plt.subplots(
        1, 2, figsize=(arscale * (4 * R * ar + 0.5), arscale * H * ar)
    )
    for R, H, ar in [(R, H, sqrt(16.0 / (R * H))) for R, H, _ in cases]
])
for i, ((shp, phi2), fn) in enumerate(zip(pcases, fns)):
    def solfunc(r, z, n=nspatial, k=kspatial):
        return np.vectorize(shp.ptwise)(
            phi2, {'ntrunc': n, 'ktrunc': k}, r, z
        )
    pkl = pkgutil.get_data('mevlib', fn)
    refsoln = pickle.loads(pkl)
    oursoln = solfunc(refsoln['rmesh'], refsoln['zmesh'])
    versoln = refsoln['umesh']
    hi1, hi2 = np.array(refsoln['rmesh'].shape) - 1
    l2errs[i] = [
        np.sqrt(sum([
            (solfunc(r, z, n=n, k=2*n) - u)**2
            for r, z, u in zip(*map(lambda x: x.flatten(), [
                refsoln['rmesh'][1:hi1, 1:hi2],
                refsoln['zmesh'][1:hi1, 1:hi2],
                refsoln['umesh'][1:hi1, 1:hi2]
            ]))
        ]))
        for n in Ns
    ]
    wterrs[i] = [
        np.sqrt(sum([
            2 * np.pi * r * (solfunc(r, z, n=n, k=2*n) - u)**2
            for r, z, u in zip(*map(lambda x: x.flatten(), [
                refsoln['rmesh'][1:hi1, 1:hi2],
                refsoln['zmesh'][1:hi1, 1:hi2],
                refsoln['umesh'][1:hi1, 1:hi2]
            ]))
        ]))
        for n in Ns
    ]
    meldsoln = np.vstack((np.flip(versoln, axis=0), oursoln))
    meldrmesh = np.vstack((
        refsoln['rmesh'] - shp.R, refsoln['rmesh']
    ))
    meldzmesh = np.vstack((refsoln['zmesh'], refsoln['zmesh']))
    for lax, rax in [
        (multaxs[i, 0], multaxs[i, 1]),
        (singaxs[i][0], singaxs[i][1])
    ]:
        cms = [None, None]
        cms[0] = lax.contourf(meldrmesh, meldzmesh, meldsoln)
        lax.axvline(x=0.0, color='r')
        err = np.abs(oursoln - versoln) / versoln
        print("Case {} error range ({}, {})".format(i+1, err.min(), err.max()))
        lvls = np.power(10.0, np.arange(
            np.floor(np.log10(err.min()) - 1), np.ceil(np.log10(err.max()) + 1)
        ))
        cms[1] = rax.contourf(
            refsoln['rmesh'], refsoln['zmesh'], err,
            lvls, norm=colors.LogNorm(),
            cmap='winter'
        )
        # lax.set_title("Solution Comparison")
        # rax.set_title("Relative Error")
        for ax in [lax, rax]:
            ax.set_xlabel(r"\( r \)")
            ax.set_ylabel(r"\( z \)")
        # these colorbars are the same for each set
        cbs = [plt.colorbar(cm, ax=ax) for ax, cm in zip([lax, rax], cms)]
    for ax in multaxs[i, :]:
        ax.set_xlabel(r"\( r \)")
        ax.set_ylabel(r"\( z \)")
multfig.tight_layout()
for fig in singfigs:
    fig.tight_layout()


l2fig, l2ax = plt.subplots(1, 1, figsize=figsize(2, 2))
for i, (errs, mkr) in enumerate(zip(l2errs, ['x', 'o', '+', '<', '>'])):
    l2ax.semilogy(
        Ns, errs, 'C{}'.format(i+1) + mkr, label="Case {}".format(i + 1)
    )
l2ax.grid()
l2ax.set_xlabel(r"Number of terms \( N \) (\( K = 2 N \))")
l2ax.set_ylabel(r"\( L^2 \) error")
l2ax.legend()

wtfig, wtax = plt.subplots(1, 1, figsize=figsize(2, 2))
for i, (errs, mkr) in enumerate(zip(wterrs, ['x', 'o', '+', '<', '>'])):
    wtax.semilogy(
        Ns, errs, 'C{}'.format(i+1) + mkr, label="Case {}".format(i + 1)
    )
wtax.grid()
wtax.set_xlabel(r"Number of terms \( N \) (\( K = 2 N \))")
wtax.set_ylabel(r"\( L^2 \) error")
wtax.legend()


for ext in ['svg', 'pdf']:
    l2fig.savefig("img/cyl_ssverif_volume_compare_l2errs.{}".format(ext))
    wtfig.savefig("img/cyl_ssverif_volume_compare_wterrs.{}".format(ext))
    multfig.savefig("img/cyl_ssverif_volume_compare_cases.{}".format(ext))
    for i, fig in enumerate(singfigs):
        fig.savefig(
            "img/cyl_ssverif_volume_compare_case{}.{}".format(i+1, ext)
        )


