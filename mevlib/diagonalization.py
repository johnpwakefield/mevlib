"""

Contains functions that yield pointwise and integrated solutions to the scalar
problem for a variety of shapes.

"""

import numpy as np
from scipy.linalg import solve
from scipy.sparse.linalg import gmres, spsolve
from scipy.sparse import dok_matrix as spmat
from scipy.interpolate import RegularGridInterpolator as RGI

# the scipy linear algebra package has some advantages, but this can but
# substituted for the numpy version if scipy is not available
from scipy import linalg as la


# diagonalization

def computetransform(shpe, mech, T, precision, dropnonreacting=False):
    #TODO interrogate the Mechanism object to choose an ideal method
    lams, R = la.eig(mech.getmatrix(shpe.charlen, T))
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


