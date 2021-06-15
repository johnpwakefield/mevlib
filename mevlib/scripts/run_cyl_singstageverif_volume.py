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


nexact, kexact = 64, 32


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
    for i in range(5)
]
rates = [1.5117e-3, 1.5117e-3, 6.0470e-3, 1.5117e-4, 1.5117e-3]
compphi2s = [rate / D for rate in rates]
cases = [
    (Cylinder(case[0], case[1]), compphi2)
    for case, compphi2 in zip(cases, compphi2s)
]
Ns = np.arange(1, 42, 4)
l2errs = [None for c in cases]  # errors in L2 norm
wterrs = [None for c in cases]  # errors in weighted norm


multfig, multaxs = plt.subplots(5, 2, figsize=(8.0, 14.5))
singfigs, singaxs = zip(*[
    plt.subplots(1, 2, figsize=(9.0, 4.0)) for i in range(5)
])
for i, ((shp, phi2), fn) in enumerate(zip(cases, fns)):
    def solfunc(r, z, n=nexact, k=kexact):
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
                refsoln['rmesh'][1:hi1,1:hi2],
                refsoln['zmesh'][1:hi1,1:hi2],
                refsoln['umesh'][1:hi1,1:hi2]
            ]))
        ]))
        for n in Ns
    ]
    wterrs[i] = [
        np.sqrt(sum([
            2 * np.pi * r * (solfunc(r, z, n=n, k=2*n) - u)**2
            for r, z, u in zip(*map(lambda x: x.flatten(), [
                refsoln['rmesh'][1:hi1,1:hi2],
                refsoln['zmesh'][1:hi1,1:hi2],
                refsoln['umesh'][1:hi1,1:hi2]
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
        lax.set_title("Solution Comparison")
        rax.set_title("Relative Error")
        # these colorbars are the same for each set
        cbs = [plt.colorbar(cm, ax=ax) for ax, cm in zip([lax, rax], cms)]
    for ax in multaxs[i, :]:
        ax.set_xlabel(r"\( r \)")
        ax.set_ylabel(r"\( z \)")
multfig.tight_layout()


l2fig, l2ax = plt.subplots(1, 1, figsize=(6.0, 6.0))
for i, (errs, mkr) in enumerate(zip(l2errs, ['x', 'o', '+', '<', '>'])):
    l2ax.semilogy(
        Ns, errs, 'C{}'.format(i+1) + mkr, label="Case {}".format(i + 1)
    )
l2ax.grid()
l2ax.set_xlabel(r"Number of terms \( N \) (\( K = 2 N \))")
l2ax.set_ylabel(r"\( L^2 \) error")
l2ax.legend()

wtfig, wtax = plt.subplots(1, 1, figsize=(6.0, 6.0))
for i, (errs, mkr) in enumerate(zip(wterrs, ['x', 'o', '+', '<', '>'])):
    wtax.semilogy(
        Ns, errs, 'C{}'.format(i+1) + mkr, label="Case {}".format(i + 1)
    )
wtax.grid()
wtax.set_xlabel(r"Number of terms \( N \) (\( K = 2 N \))")
wtax.set_ylabel(r"\( L^2 \) error")
wtax.legend()


if showfigs():
    plt.show()
else:
    for ext in ['svg', 'pdf']:
        l2fig.savefig(imgpath("comparison_l2errs.{}".format(ext)))
        wtfig.savefig(imgpath("comparison_wterrs.{}".format(ext)))
        multfig.savefig(imgpath("comparison_cases.{}".format(ext)))
        for i, fig in enumerate(singfigs):
            fig.savefig(imgpath("comparison_case{}.{}".format(i+1, ext)))


