

import numpy as np
from scipy.io import savemat


def mat(outfile, matrices, temperatures, verb=False):
    print("Generating a .mat file for use with 'mevlookup.m'...")
    savemat(
        outfile,
        {'mats': np.array(matrices), 'temps': np.array(temperatures)}
    )


