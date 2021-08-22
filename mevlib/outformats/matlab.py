

import numpy as np
from scipy.io import savemat


# TODO should really just be an `all' since the mat format is pretty flexible

def write_mat_mat(outfile, spcsyms, temperatures, matrices, verb=False):
    print("Generating a .mat file for use with 'mevlookup.m'...")
    savemat(
        outfile,
        {'mats': np.array(matrices), 'temps': np.array(temperatures)}
    )


