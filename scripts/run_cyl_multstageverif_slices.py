#!/usr/bin/env python3


import pkgutil
from io import StringIO
import pickle

import numpy as np

from matplotlib import pyplot as plt

from mevlib.mechanisms import Mechanism
from mevlib.shapes import Cylinder
from mevlib.diagonalization import diag_ptwise_setup, diag_ptwise_eval
from mevlib.parsing.auto import parse_attempt


plt.rc('font', size=10)
plt.rc('text', usetex=True)
plt.rc('axes', labelsize=10)
plt.rc('legend', fontsize=10)

axsize = (2.75, 2.5)

def figsize(i, j):
    return (axsize[0] * j, axsize[1] * i)

markersize = 2


cases = [
    ((180.0, 360.0), 600.0, (36.0, 32.0, 28.0, 24.0, 8.0, 0.0)),
    ((180.0, 360.0), 800.0, (36.0, 32.0, 28.0, 24.0, 8.0, 0.0)),
    ((180.0, 360.0), 600.0, (12.0, 16.0, 20.0, 24.0, 16.0, 0.0))
]
mechfiles = [   # this is a hack because of the issue with the verif data
    'data/fcc_multistageverif_36.sbl',
    'data/fcc_multistageverif_36.sbl',
    'data/fcc_multistageverif_12.sbl'
]
drs = ["H", "rad"]
components = ['S', 'D', 'Gas', 'LPG', 'DR']
reffile = "data/multistage_cyl.pickle"


refdata = pickle.loads(pkgutil.get_data('mevlib', reffile))
cases = [(Cylinder(*dims), T, np.array(bdry)) for dims, T, bdry in cases]
figs, axs = zip(*[
    plt.subplots(5, 2, figsize=figsize(5, 2)) for i in range(len(cases))
])
for i, ((cyl, temp, bdry), mechfile) in enumerate(zip(cases, mechfiles)):
    # parse test mechanism
    mechconf = pkgutil.get_data('mevlib', mechfile)
    precision, shape, temperatures, species, reactions = parse_attempt(
        StringIO(mechconf.decode('utf-8')), '.sbl', True, True
    )
    if precision is None:
        precision = {'ntrunc': 128, 'ktrunc': 64}
    mech = Mechanism(species, reactions)
    for k, dr in enumerate(drs):
        # compute our solution, plot slices
        refsoln = refdata["case{}_{}".format(i+1, dr.lower())]
        for j, ref in enumerate([refsoln[k] for k in components]):
            diagdata = diag_ptwise_setup(cyl, mech, bdry, temp, precision)
            xs = ref['x'][1:-1]
            rfys = ref['y'][1:-1]
            if dr == "H":
                our = [
                    diag_ptwise_eval(diagdata, 0.0, 1e6 * z + cyl.H / 2)
                    for z in xs
                ]
            else:
                our = [
                    diag_ptwise_eval(diagdata, 1e6 * r, cyl.H / 2)
                    for r in xs
                ]
            axs[i][j, k].plot(
                1e6 * np.array(xs), rfys, 'bx',
                label="Reference Solution", markersize=markersize
            )
            axs[i][j, k].plot(
                1e6 * np.array(xs), [y[j] for y in our], 'g+',
                label="Our Solution", markersize=markersize
            )
            if dr == "H":
                axs[i][j, k].set_xlabel(r"\( z \)")
            else:
                axs[i][j, k].set_xlabel(r"\( r \)")
            axs[i][j, k].set_ylabel(r"\( C_{} \)".format(j+1))
            axs[i][j, k].grid()
            refrange = max(ref['y']) - min(xs)
            if False:
                axs[i][j, k].set_ylim((
                    min(ref['y']) - 0.1 * refrange,
                    max(ref['y']) + 0.1 * refrange
                ))
    axs[i][-1, -1].legend()


for i, fig in enumerate(figs):
    fig.tight_layout()
    for ext in ['svg', 'pdf']:
        fig.savefig("img/cyl_multistage_case{}.{}".format(i+1, ext))


