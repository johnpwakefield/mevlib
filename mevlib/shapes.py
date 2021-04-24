

# TODO it is probably a good idea to keep the functions in 'scalar' separate
# from containing objects for diagnostic purposes, but those functions should
# be used here so we are not duplicating code.


from abc import ABC

import numpy as np
from scipy.special import i0e, i1e, jn_zeros

# the scipy linear algebra package has some advantages, but this can but
# substituted for the numpy version if scipy is not available
from scipy import linalg as la


class Shape(ABC):
    pass

class Cylinder(Shape):

    def __init__(self, R, H, charlen=None):
        if charlen is None:
            charlen = 0.5 * (R * H / (R + H))
        self.R, self.H, self.charlen, self.zero_cache = R, H, charlen, None

    def gn(self, a2, n):
        return np.sqrt(a2 + (np.pi * (2 * n + 1) / self.H)**2)

    def ln(self, a2, alphak):
        return np.sqrt(a2 + (alphak / self.R)**2)

    def s1(self, a2, trunc, bratcutoff=None):
        ns = np.array([n for n in range(trunc-1, -1, -1)])
        gns = self.gn(a2, ns)
        # the ratio I1/I0 tends to 1 (inf/inf) very quickly
        brat = i1e(gns * self.R) / i0e(gns * self.R)
        if bratcutoff is not None:
            brat = np.where(gns * self.R < bratcutoff, brat, 1.0)
        return sum(16.0 / (np.pi**2 * self.R * gns * (2 * ns + 1)**2) * brat)

    def s2(self, a2, trunc):
        if self.zero_cache is None or len(self.zero_cache) < trunc:
            self.zero_cache = jn_zeros(0, trunc)
        lnks = self.ln(a2, self.zero_cache[::-1])
        return 8.0 * sum(
            np.tanh(lnks * self.H / 2)
            / (self.H * lnks * self.zero_cache[::-1]**2)
        )

    def s(self, a2, ntrunc, ktrunc):
        s1val = self.s1(a2, ntrunc)
        s2val = self.s2(a2, ktrunc)
        if np.any(np.logical_not(np.isfinite(s1val))):
            print("Nonfinite values in s1.")
            print(s1val)
        if np.any(np.logical_not(np.isfinite(s2val))):
            print("Nonfinite values in s2.")
            print(s2val)
        return s1val + s2val

    def computetransform(self, mech, T, ntrunc, ktrunc, dropnonreacting=False):
        #TODO interrogate the Mechanism object to choose an ideal method
        lams, R = la.eig(mech.getmatrix(self.charlen, T))
        if np.max(np.abs(lams.imag)) > 1e-8:
            print("Matrix has imaginary eigenvalues (irreversible).")
        lams, R = np.real_if_close(lams), np.real_if_close(R)
        mult = np.array([
            1.0 if lam == 0.0 else self.s(lam, ntrunc, ktrunc) for lam in lams
        ])
        A = np.dot(R, np.dot(np.diag(mult), la.inv(R)))
        if dropnonreacting:
            return A[:mech.getnumactive(), :]
        else:
            return A

