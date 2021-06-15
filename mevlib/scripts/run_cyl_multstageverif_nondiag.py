#!/usr/bin/env python3


import pkgutil
from io import StringIO
import pickle

import numpy as np

from matplotlib import pyplot as plt

from mevlib.mechanisms import Mechanism
from mevlib.diagonalization import diag_ptwise_setup
from mevlib.shapes import Cylinder
from mevlib.options import showfigs
from mevlib.parsing.auto import parse_attempt


plt.rc('font', size=12)
plt.rc('text', usetex=True)
plt.rc('axes', labelsize=18)


cases = [
    ((180.0, 360.0), 600.0, (36.0, 32.0, 28.0, 24.0, 8.0, 0.0)),
    ((180.0, 360.0), 800.0, (36.0, 32.0, 28.0, 24.0, 8.0, 0.0)),
    ((180.0, 360.0), 600.0, (12.0, 16.0, 20.0, 24.0, 16.0, 0.0))
]
mechfiles = [ # this is a hack because of the issue with the verif data
    'data/fcc_multistageverif_36.sbl',
    'data/fcc_multistageverif_36.sbl',
    'data/fcc_multistageverif_12.sbl'
]
drs = ["H", "rad"]
components = ['S', 'D', 'Gas', 'LPG', 'DR']
reffile = "data/multistage_cyl.pickle"


refdata = pickle.loads(pkgutil.get_data('mevlib', reffile))
cases = [(Cylinder(*dims), T, np.array(bdry)) for dims, T, bdry in cases]
l2errs = [None for c in cases]
figs, axs = zip(*[
    plt.subplots(1, 2, figsize=(8.0, 14.5)) for i in range(len(cases))
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
    # compute phi2
    # note this is for fcc_lattanzietal_2020 specifically
    refkijsum = 0.069537
    refDi = 5.138e-9
    refphi2 = refkijsum / refDi
    libDis = mech.getDis(temp)
    libkijs = mech.getkijs(temp)
    np.set_printoptions(formatter={'float': "{0:0.4e}".format})
    print("library computed Dis are {}.".format(libDis))
    print("library computed sum of kijs is {}.".format(
        np.sum(libkijs, axis=0))
    )
    phi2 = sum(libkijs[:, 0] / libDis[0])
    # compare this phi2 to the one provided by diagonalization class
    _, lams, _, _, _ = diag_ptwise_setup(cyl, mech, bdry, temp, precision)
    print(lams) #TODO fix this
    print("library phi2 is {}".format(max(lams)))
    print("manual phi2 is {}.".format(phi2))
    for k, dr in enumerate(drs):
        # compute our solution, plot slices
        ref = refdata["case{}_{}".format(i+1, dr.lower())][components[0]]
        if dr == "H":
            our = [
                bdry * cyl.ptwise(phi2, precision, 0.0, 1e6 * z + cyl.H / 2)
                for z in ref['x']
            ]
        else:
            our = [
                bdry * cyl.ptwise(phi2, precision, 1e6 * r, cyl.H / 2)
                for r in ref['x']
            ]
        axs[i][k].plot(
            1e6 * np.array(ref['x']), ref['y'], 'bx',
            label="Reference Solution"
        )
        axs[i][k].plot(
            1e6 * np.array(ref['x']), [y[0] for y in our], 'g+',
            label="Our Solution"
        )
        if dr == "H":
            axs[i][k].set_xlabel(r"\( z \)")
        else:
            axs[i][k].set_xlabel(r"\( r \)")
        axs[i][k].grid()
        refrange = max(ref['y']) - min(ref['x'])
        if False:
            axs[i][k].set_ylim((
                min(ref['y']) - 0.1 * refrange,
                max(ref['y']) + 0.1 * refrange
            ))
axs[-1][-1].legend()
for fig in figs:
    fig.tight_layout()


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


