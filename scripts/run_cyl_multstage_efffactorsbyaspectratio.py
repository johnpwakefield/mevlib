#!/usr/bin/env python3


import pkgutil
from io import StringIO
from functools import reduce

import numpy as np
import matplotlib.pyplot as plt

from mevlib.mechanisms import Mechanism
from mevlib.shapes import Cylinder
from mevlib.diagonalization import PointDiagonalization, computeintegral
from mevlib.parsing.auto import parse_attempt


plt.rc('font', size=10)
plt.rc('text', usetex=True)
plt.rc('axes', labelsize=10)
plt.rc('legend', fontsize=10)

axsize = (2.75, 2.5)

def figsize(i, j):
    return (axsize[0] * j, axsize[1] * i)


T = 800.0
bdry = np.array([0.3, 0.1, 0.1, 0.0, 0.0, 0.0])
mechfile = 'data/fcc_lattanzietal_2020.sbl'
components = ['S', 'D', 'G', 'LPG', 'DR']

Npts = 201
ratios = 10**np.linspace(-2.0, 4.0, Npts)
V = np.pi * (200.0)**3
Rs = (V / (np.pi * ratios))**(1.0/3)
Hs = ratios * Rs

# parse mechanism
mechconf = pkgutil.get_data('mevlib', mechfile)
precision, _, temperatures, tempspacing, species, reactions = parse_attempt(
    StringIO(mechconf.decode('utf-8')), '.sbl', True, True
)
if precision is None:
    precision = {'ntrunc': 128, 'ktrunc': 64}
mech = Mechanism(species, reactions)
diag = PointDiagonalization(mech, T)

# build data arrays
MEVs = np.empty((Npts, 5))
for i, (R, H) in enumerate(zip(Rs, Hs)):
    cyl = Cylinder(R, H)
    MEVs[i, :] = reduce(np.dot, [
        diag.get_evects(),
        np.diag(computeintegral(cyl, diag.get_evals(), precision)),
        diag.get_evectinv(),
        bdry[:mech.Ng]
    ]) / bdry[:mech.Ng]

# plot
fig, axs = plt.subplots(1, 3, figsize=figsize(1, 3))
for i in range(3):
    axs[i].semilogx(Hs / (2 * Rs), MEVs[:, i], 'C{}-'.format(i))
    axs[i].set_xlabel(r'\( H / (2 R) \)')
    axs[i].set_ylabel(r'\( \eta_{} \)'.format(i+1))
    # axs[i].set_title(r'\( {} \)'.format(components[i]))
    axs[i].grid()
fig.tight_layout()

# save
for ext in ['svg', 'pdf']:
    fig.savefig("img/cyl_multstage_efffactorbyaspectratio.{}".format(ext))


