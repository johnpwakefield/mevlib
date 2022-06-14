#!/usr/bin/env python3


import pkgutil
from io import StringIO

import numpy as np

from mevlib.mechanisms import Mechanism
from mevlib.diagonalization import PointDiagonalization
from mevlib.shapes import Sphere
from mevlib.parsing.auto import parse_attempt


np.set_printoptions(formatter={"float_kind": "{:.4e}".format})


T = 800.0


mechfile = 'data/fcc_lattanzietal_2020.sbl'


# load file

mechconf = pkgutil.get_data('mevlib', mechfile)
precision, shape, temperatures, _, species, reactions = parse_attempt(
    StringIO(mechconf.decode('utf-8')), '.sbl', True, True
)
if precision is None:
    precision = {'ntrunc': 128, 'ktrunc': 64}
mech = Mechanism(species, reactions)


# print reaction rates

print("reaction rates")
print(mech.getkijs(T))


# print diffusion coefficients

print("diffusion coeffs")
print(mech.getDis(T) * 7.0 / 0.319)


# print thiele matrix

print("thiele matrix")
print(mech.getB(T))


# print effectiveness factors

Lam, R = np.linalg.eig(mech.getB(T)[:mech.Ng, :mech.Ng])
print("effectiveness factors")
print([Sphere(200.0).intgtd(a2) for a2 in Lam])


# print conservation check

print("conservation check")
diag = PointDiagonalization(mech, T)
print(np.dot(mech.getDis(T).reshape((1, -1)), np.vstack((
    np.dot(diag.get_evects(), np.diag(diag.get_evals())),
    np.dot(mech.getB(T)[-1, :], diag.get_evects())
))))



