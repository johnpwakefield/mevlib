#!/usr/bin/env python3


import pkgutil
from io import StringIO

import numpy as np
from scipy.special import jn_zeros

from mevlib.mechanisms import Mechanism
from mevlib.parsing.auto import parse_attempt


np.set_printoptions(formatter={"float_kind": "{:.4e}".format})


mechfile = 'data/fcc_lattanzietal_2020.sbl'
components = ['S', 'D', 'Gas', 'LPG', 'DR', 'CK']
reffile = "data/multistage_cyl.pickle"
T = 800.0


mechconf = pkgutil.get_data('mevlib', mechfile)
precision, shape, temperatures, species, reactions = parse_attempt(
    StringIO(mechconf.decode('utf-8')), '.sbl', True, True
)
mech = Mechanism(species, reactions)


# TODO this should be rolled into the mevtable entry point as an option


print("Dis at T = {}:".format(T))
print([spc.effective_diffusion(T) for spc in species])


print("kijs at T = {}:".format(T))
print(mech.getkijs(T))


print("Full Matrix at T = {}:".format(T))
print(mech.getmatrix(T))


print("First 300 zeros of J_0L")
print(jn_zeros(0, 300))


