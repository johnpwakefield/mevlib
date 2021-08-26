

import numpy as np
from scipy.io import savemat


# TODO should really just be an `all' since the mat format is pretty flexible

def w_mat_mat(outfile, spcsyms, temperatures, matrices, verb=False):
    if verb:
        print("Generating a .mat file for use with 'mevlookup.m'...")
    savemat(
        outfile,
        {'mats': np.array(matrices), 'temps': np.array(temperatures)}
    )

def w_int_mat(outfile, lambdas, ints, verb=False):
    if verb:
        print("Generating a .mat file containing single-stage reaction rates.")
    savemat(
        outfile,
        {'lambdas': np.array(lambdas), 'integrals': np.array(ints)}
    )

def w_dia_mat(outfile, syms, Ts, lambdas, Rs, Rinvs, verb=False):
    if verb:
        print("Generating a .mat file containing diagonalization information.")
    n, _ = Rs[0].shape
    lamat = np.empty((len(Rs), n))
    Rsmat = np.empty((len(Rs), n, n))
    Rimat = np.empty((len(Rs), n, n))
    for m in range(len(Rs)):
        lamat[m, :] = lambdas[m]
        Rsmat[m, :, :], Rimat[m, :, :] = Rs[m], Rinvs[m]
    savemat(outfile, {
        'syms': syms,
        'temps': np.array(Ts),
        'lambdas': lamat,
        'Rs': Rsmat,
        'Rinvs': Rimat
    })


