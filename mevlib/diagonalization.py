"""

Contains functions that yield pointwise and integrated solutions to the scalar
problem for a variety of shapes.

"""

import numpy as np

# the scipy linear algebra package has some advantages, but this can but
# substituted for the numpy version if scipy is not available
from scipy import linalg as la


# matrix transforms for make_table

def computetransform(shpe, mech, T, precision, dropnonreacting=False):
    #TODO interrogate the Mechanism object to choose an ideal method
    lams, R = la.eig(mech.getmatrix(T))
    if np.max(np.abs(lams.imag)) > 1e-8:
        print("Matrix has imaginary eigenvalues (irreversible).")
    lams, R = np.real_if_close(lams), np.real_if_close(R)
    mult = np.array([
        1.0 if lam == 0.0 else shpe.intgtd(lam, precision) for lam in lams
    ])
    A = np.dot(R, np.dot(np.diag(mult), la.inv(R)))
    if dropnonreacting:
        return A[:mech.getnumactive(), :]
    else:
        return A


# pointwise transform (setup and compute functions so we don't recompute at
# each point)

def diag_ptwise_setup(shpe, mech, T, precision):
    #TODO interrogate the Mechanism object to choose an ideal method
    lams, R = la.eig(mech.getmatrix(T))
    if np.max(np.abs(lams.imag)) > 1e-8:
        print("Matrix has imaginary eigenvalues (irreversible).")
    lams, R = np.real_if_close(lams), np.real_if_close(R)
    return shpe, lams, R, precision

def diag_ptwise_eval(setupdata, bdry, *coords):
    shpe, lams, R, precision = setupdata
    mult = np.array([
        1.0 if lam == 0.0 else shpe.ptwise(lam, precision, *coords)
        for lam in lams
    ])
    return np.dot(R, np.dot(np.diag(mult), la.solve(R, bdry.reshape(-1,1))))


