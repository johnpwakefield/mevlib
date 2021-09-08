

import numpy as np
from pickle import dump


# TODO should really be an `all' since the pickle format is pretty flexible

def w_rate_pkl(outfile, diagset, shape, precision, verb=False):
    print("Generating a pickle file for use with 'mevlookup.py'...")
    with open(outfile, 'wb') as f:
        dump({
            'mats':  np.array(diagset.get_transforms(shape, precision)),
            'temps': np.array(diagset.Ts)
        }, f)

