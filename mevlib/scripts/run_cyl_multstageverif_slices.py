#!/usr/bin/env python3


import pkgutil
from io import StringIO
import pickle

import numpy as np

from matplotlib import pyplot as plt

from mevlib.mechanisms import Mechanism
from mevlib.shapes import Cylinder
from mevlib.diagonalization import diag_ptwise_setup, diag_ptwise_eval
from mevlib.options import showfigs
from mevlib.parsing.auto import parse_attempt


plt.rc('font', size=12)
plt.rc('text', usetex=True)
plt.rc('axes', labelsize=18)


D = 5.4423


cases = [
#   ((180.0, 360.0), 600.0, (0.9, 0.8, 0.7, 0.6, 0.2, 0.2)),
#   ((180.0, 360.0), 800.0, (0.9, 0.8, 0.7, 0.6, 0.2, 0.2)),
    ((180.0, 360.0), 600.0, (0.3, 0.4, 0.5, 0.6, 0.4, 0.3))
]
drs = ["H", "rad"]
components = ['S', 'D', 'Gas', 'LPG', 'DR']


# build shapes
cases = [(Cylinder(*dims), T, np.array(bdry)) for dims, T, bdry in cases]


# parse test mechanism
mechconf = pkgutil.get_data('mevlib', 'data/verificationcases_mechanism_2.sbl')
precision, shape, temperatures, species, reactions = parse_attempt(
    StringIO(mechconf.decode('utf-8')), '.sbl', True, True
)
mech = Mechanism(species, reactions)
if precision is None:
    precision = {'ntrunc': 128, 'ktrunc': 64}


# compute and make plots
l2errs = [None for c in cases]
figs, axs = zip(*[
    plt.subplots(5, 2, figsize=(8.0, 14.5)) for i in range(len(cases))
])
for i, (shape, temp, bdry) in enumerate(cases):
    for k, dr in enumerate(drs):
        # read reference file
        fn = "data/multistage_cyl_case{}_{}.pickle".format(i+1, dr)
        refsoln = pickle.loads(pkgutil.get_data('mevlib', fn))
        # compute our solution, plot slices
        for j, ref in enumerate([refsoln[k] for k in components]):
            diagdata = diag_ptwise_setup(shape, mech, temp, precision)
            if dr == "H":
                our = [
                    diag_ptwise_eval(
                        diagdata, 40.0 * bdry, 0.0, 1e6 * z + shape.H / 2
                    )
                    for z in ref['x']
                ]
            else:
                our = [
                    diag_ptwise_eval(
                        diagdata, 40.0 * bdry, 1e6 * r, shape.H / 2
                    )
                    for r in ref['x']
                ]
            axs[i][j,k].plot(
                1e6 * np.array(ref['x']), ref['y'], 'bx',
                label="Reference Solution"
            )
            axs[i][j,k].plot(
                1e6 * np.array(ref['x']), [y[j] for y in our], 'g+',
                label="Our Solution"
            )
            refrange = max(ref['y']) - min(ref['x'])
            if False:
                axs[i][j,k].set_ylim((
                    min(ref['y']) - 0.1 * refrange,
                    max(ref['y']) + 0.1 * refrange
                ))


if showfigs():
    plt.show()
else:
    for ext in ['svg', 'pdf']:
        pass
        #TODO
#       l2fig.savefig(imgpath("comparison_l2errs.{}".format(ext)))
#       multfig.savefig(imgpath("comparison_cases.{}".format(ext)))
#       for i, fig in enumerate(singfigs):
#           fig.savefig(imgpath("comparison_case{}.{}".format(i+1, ext)))


