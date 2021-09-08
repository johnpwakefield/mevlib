
from functools import reduce

import numpy as np

# the scipy linear algebra package has some advantages, but this can but
# substituted for the numpy version if scipy is not available
from scipy import linalg as la


# compute diagonalization of the matrix Bg
# TODO at some point in the future this should be replaced by a collection of
# methods that leverage any properties Bg may have

def diagonalize(Bg, reversible=False):
    lams, R = la.eig(Bg)
    if np.max(np.abs(lams.imag)) > 1e-8:
        print("Matrix has imaginary eigenvalues (irreversible).")
    lams, R = np.real_if_close(lams), np.real_if_close(R)
    if reversible:
        pairs = [(lam, R[:, i].copy()) for i, lam in enumerate(lams)]
        pairs = sorted(
            pairs, key=lambda tpl: np.argwhere(np.abs(tpl[1]) > 1e-6)[0]
        )
        for i, (l, r) in enumerate(pairs):
            lams[i] = l
            R[:, i] = r
    return lams, R


# functions for computing integrals and rate transforms

def computeintegral(shpe, lams, precision):
    return np.array([
        1.0 if lam == 0.0 else shpe.intgtd(lam, precision) for lam in lams
        ])

def computetransform(shpe, DboL2, lams, R, Bs, Rinv, precision):
    assert(R.shape[1] == Bs.shape[1])
    assert(R.shape == Rinv.shape)
    Ng, Ns = R.shape[0], Bs.shape[0]
    A = np.empty((Ng + Ns, Ng))
    tail = np.dot(np.diag(computeintegral(shpe, lams, precision)), Rinv)
    A[:Ng, :] = reduce(np.dot, [np.diag(DboL2), R, np.diag(lams), tail])
    if Ns > 0:
        A[-Ns:, :] = reduce(np.dot, [np.diag(DboL2), Bs, R, tail])
    return A


# diagonalization objects

class PointDiagonalization(object):

    Bg, Bs = None, None
    evals, evects, evectinv = None, None, None

    def __init__(self, mech, T):
        self.mech, self.T = mech, T
        B = mech.getB(T)
        self.Bg = B[:mech.Ng, :]
        self.Bs = B[-mech.Ns:, :] if mech.Ns > 0 else np.empty((0, mech.Ng))
        self.evals, self.evects = diagonalize(
            self.Bg, reversible=mech.isreversible()
        )

    def get_evals(self):
        return self.evals

    def get_evects(self):
        return self.evects

    def get_evectinv(self):
        if self.evectinv is None:
            self.evectinv = self.evects.inv()
        return self.evectinv


class DiagonalizationSet(object):

    Bgs, Bss = None, None
    evals, evects, evectinvs = None, None, None

    def __init__(self, mech, Ts):
        self.mech, self.Ts = mech, Ts
        Ng, Ns = self.mech.Ng, self.mech.Ns
        Bs = [self.mech.getB(T) for T in self.Ts]
        assert(all([B.shape == (Ng + Ns, Ng) for B in Bs]))
        self.Bgs = [B[:Ng, :] for B in Bs]
        if Ns > 0:
            self.Bss = [B[-Ns:, :] for B in Bs]
        else:
            self.Bss = [np.empty((0, mech.Ng)) for B in Bs]

    def do_diag(self):
        # TODO this should take advantage of properties of the mechanism
        assert(self.mech.isvalid())
        self.evals, self.evects = zip(*[
            diagonalize(Bg, reversible=self.mech.isreversible)
            for Bg in self.Bgs
        ])

    def do_invs(self):
        if self.evals is None or self.evects is None:
            self.do_diag()
        if self.evectinvs is None:
            self.evectinvs = [la.inv(R) for R in self.evects]

    def get_Dis(self):
        return [self.mech.getDis(T) for T in self.Ts]

    def get_evals(self):
        if self.evals is None:
            self.do_diag()
        return self.evals

    def get_evects(self):
        if self.evects is None:
            self.do_diag()
        return self.evects

    def get_evectinvs(self):
        if self.evectinvs is None:
            self.do_invs()
        return self.evectinvs

    def get_transforms(self, shpe, precision):
        return [
            computetransform(
                shpe, D, evals, evects, Bs, Rinv, precision
            )
            for D, evals, evects, Bs, Rinv in zip(
                self.get_Dis(), self.get_evals(), self.get_evects(),
                self.Bss, self.get_evectinvs()
            )
        ]


# TODO deprecate or clean these up
# pointwise transform (setup and compute functions so we don't recompute at
# each point)

def diag_ptwise_setup(shpe, mech, bdry, T, precision):
    # TODO interrogate the Mechanism object to choose an ideal method
    lams, R = la.eig(mech.getmatrix(T))
    if np.max(np.abs(lams.imag)) > 1e-8:
        print("Matrix has imaginary eigenvalues (irreversible).")
    lams, R = np.real_if_close(lams), np.real_if_close(R)
    ubdry = la.solve(R, bdry.reshape(-1, 1))
    return shpe, lams, R, ubdry, precision

def diag_ptwise_eval(setupdata, *coords):
    shpe, lams, R, ubdry, precision = setupdata
    mult = np.array([
        1.0 if lam == 0.0 else shpe.ptwise(lam, precision, *coords)
        for lam in lams
    ])
    return np.dot(R, np.dot(np.diag(mult), ubdry))


