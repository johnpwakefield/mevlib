

import numpy as np
from scipy.io import savemat


# TODO should really just be an `all' since the mat format is pretty flexible

def w_rate_mat(outfile, diagset, shape, precision, verb=False):
    if verb:
        print("Generating a .mat file for use with 'mevlookup.m'...")
    savemat(outfile, {
        'mats':  np.array(diagset.get_transforms(shape, precision)),
        'temps': np.array(diagset.Ts)
    })

def w_ints_mat(outfile, lams, ints, verb=False):
    if verb:
        print("Generating a .mat file containing single-stage reaction rates.")
    savemat(outfile, {'lambdas': np.array(lams), 'integrals': np.array(ints)})

def w_diag_mat(outfile, diagset, shape, precision, verb=False, inclinvs=False):
    if verb:
        print("Generating a .mat file containing diagonalization information.")
    Ng, Ns, M = diagset.mech.Ng, diagset.mech.Ns, len(diagset.Ts)
    Dimat, lamat, Rsmat, BsRmat = map(np.empty, [
        (M, Ng + Ns), (M, Ng), (M, Ng, Ng), (M, Ns, Ng)
    ])
    for j, (Di, la, R, Bs) in enumerate(zip(
        diagset.get_Dis(), diagset.get_evals(), diagset.get_evects(),
        diagset.Bss
    )):
        Dimat[j, :], lamat[j, :], Rsmat[j, :, :] = Di, la, R
        BsRmat[j, :, :] = np.dot(Bs, R)
    data = {
        'syms': [spc.symb for spc in diagset.mech.spcs],
        'temps': np.array(diagset.Ts),
        'lambdas': lamat,
        'Rs': Rsmat,
    }
    if inclinvs:
        Rimat = np.empty((M, Ng, Ng))
        for j, Ri in enumerate(diagset.get_evectinvs()):
            Rimat[j, :] = Ri
        data.update({'Rinvs': Rimat})
    savemat(outfile, data)

def w_full_mat(outfile, diagset, shape, precision, verb=False):
    w_diag_mat(outfile, diagset, shape, precision, verb=verb, inclinvs=True)

