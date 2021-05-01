

from abc import ABC

import numpy as np
import mevlib.scalar as scalar
from scipy.special import jn_zeros

# the scipy linear algebra package has some advantages, but this can but
# substituted for the numpy version if scipy is not available
from scipy import linalg as la


class Shape(ABC):
    pass

class Sphere(Shape):

    def __init__(self, R, charlen=None):
        if charlen is None:
            charlen = R / 3
        self.R, self.charlen = R, charlen

    def s(self, a2):
        return scalar.sph_intgtd(a2, self.R)

    def computetransform(self, mech, T, dropnonreacting=False):
        #TODO extract most of this diagonalization process to a function or the superclass
        #TODO interrogate the Mechanism object to choose an ideal method
        lams, R = la.eig(mech.getmatrix(self.charlen, T))
        if np.max(np.abs(lams.imag)) > 1e-8:
            print("Matrix has imaginary eigenvalues (irreversible).")
        lams, R = np.real_if_close(lams), np.real_if_close(R)
        mult = np.array([
            1.0 if lam == 0.0 else self.s(lam) for lam in lams
        ])
        A = np.dot(R, np.dot(np.diag(mult), la.inv(R)))
        if dropnonreacting:
            return A[:mech.getnumactive(), :]
        else:
            return A

class Cylinder(Shape):

    def __init__(self, R, H, charlen=None):
        if charlen is None:
            charlen = 0.5 * (R * H / (R + H))
        self.R, self.H, self.charlen, self.zero_cache = R, H, charlen, None

    def gn(self, a2, n):
        return np.sqrt(a2 + (np.pi * (2 * n + 1) / self.H)**2)

    def ln(self, a2, alphak):
        return np.sqrt(a2 + (alphak / self.R)**2)

    def s(self, a2, ntrunc, ktrunc):
        if self.zero_cache is None or len(self.zero_cache) < ktrunc:
            self.zero_cache = jn_zeros(0, ktrunc)
        return sum(map(scalar.sum_standard, [
            scalar.cyl_intgtd_radial_terms(a2, self.R, self.H, ntrunc),
            scalar.cyl_intgtd_axial_terms(
                a2, self.R, self.H, ktrunc, zero_cache=self.zero_cache
            )
        ]))

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

