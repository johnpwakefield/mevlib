#!/usr/bin/env python3


import pkgutil
from io import StringIO
import pickle

import numpy as np

from mevlib.mechanisms import Mechanism
from mevlib.parsing.auto import parse_attempt


np.set_printoptions(formatter={"float_kind": "{:.4e}".format})


mechfile = 'data/fcc_lattanzietal_2020.sbl'
components = ['S', 'D', 'Gas', 'LPG', 'DR', 'CK']
reffile = "data/multistage_cyl.pickle"
T = 800.0


refdata = pickle.loads(pkgutil.get_data('mevlib', reffile))
mechconf = pkgutil.get_data('mevlib', mechfile)
precision, shape, temperatures, species, reactions = parse_attempt(
    StringIO(mechconf.decode('utf-8')), '.sbl', True, True
)
mech = Mechanism(species, reactions)


print(mech.getmatrix(T))


