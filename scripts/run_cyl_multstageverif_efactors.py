#!/usr/bin/env python3


import pkgutil
from io import StringIO

import numpy as np

from mevlib.mechanisms import Mechanism
from mevlib.shapes import Cylinder
from mevlib.diagonalization import computetransform
from mevlib.parsing.auto import parse_attempt


cases = [
#   ((180.0, 360.0), 600.0, (36.0, 32.0, 28.0, 24.0, 8.0, 0.0)),
    ((180.0, 360.0), 800.0, (36.0, 32.0, 28.0, 24.0, 8.0, 0.0)),
#   ((180.0, 360.0), 600.0, (12.0, 16.0, 20.0, 24.0, 16.0, 0.0))
]
mechfiles = [   # this is a hack because of the issue with the verif data
#   'data/fcc_multistageverif_36.sbl',
#   'data/fcc_multistageverif_36.sbl',
#   'data/fcc_multistageverif_12.sbl'
    'data/fcc_lattanzietal_2020.sbl'
]
components = ['S', 'D', 'Gas', 'LPG', 'DR']


cases = [(Cylinder(*dims), temp, np.array(bdry)) for dims, temp, bdry in cases]
for i, ((cyl, temp, bdry), mechfile) in enumerate(zip(cases, mechfiles)):
    # parse test mechanism
    mechconf = pkgutil.get_data('mevlib', mechfile)
    precision, shape, temperatures, species, reactions = parse_attempt(
        StringIO(mechconf.decode('utf-8')), '.sbl', True, True
    )
    if precision is None:
        precision = {'ntrunc': 128, 'ktrunc': 64}
    mech = Mechanism(species, reactions)
    mev = np.dot(computetransform(cyl, mech, temp, precision), bdry) / bdry
    print("Case {} MEV: {}".format(i+1, mev))


