
import numpy as np

# the scipy linear algebra package has some advantages, but this can but
# substituted for the numpy version if scipy is not available
from scipy import linalg as la



#TODO clean these up and make sure standalone scripts still work


# matrix transforms for make_table

def computediagonalization(mech, T):
    lams, R = la.eig(mech.getmatrix(T))
    if np.max(np.abs(lams.imag)) > 1e-8:
        print("Matrix has imaginary eigenvalues (irreversible).")
    lams, R = np.real_if_close(lams), np.real_if_close(R)
    if mech.isreversible():
        pairs = [(lam, R[:,i].copy()) for i, lam in enumerate(lams)]
        pairs = sorted(pairs, key=lambda tpl: np.argwhere(np.abs(tpl[1]) > 1e-6)[0])
        for i, (l, r) in enumerate(pairs):
            lams[i] = l
            R[:,i] = r
    return lams, R

def computeintegral(shpe, lams, precision):
    return np.array([
        1.0 if lam == 0.0 else shpe.intgtd(lam, precision) for lam in lams
        ])

#TODO deprecate and remove this function
def computetransform(shpe, mech, T, precision):
    lams, R = computediagonalization(mech, T)
    mult = computeintegral(shpe, lams, precision)
    A = np.dot(R, np.dot(np.diag(mult), la.inv(R)))
    return A

# pointwise transform (setup and compute functions so we don't recompute at
# each point)

def diag_ptwise_setup(shpe, mech, bdry, T, precision):
    #TODO interrogate the Mechanism object to choose an ideal method
    lams, R = la.eig(mech.getmatrix(T))
    if np.max(np.abs(lams.imag)) > 1e-8:
        print("Matrix has imaginary eigenvalues (irreversible).")
    lams, R = np.real_if_close(lams), np.real_if_close(R)
    ubdry = la.solve(R, bdry.reshape(-1,1))
    return shpe, lams, R, ubdry, precision

def diag_ptwise_eval(setupdata, *coords):
    shpe, lams, R, ubdry, precision = setupdata
    mult = np.array([
        1.0 if lam == 0.0 else shpe.ptwise(lam, precision, *coords)
        for lam in lams
    ])
    return np.dot(R, np.dot(np.diag(mult), ubdry))


