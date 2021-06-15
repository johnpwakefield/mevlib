

from abc import ABC

#import numpy as np
from scipy.special import jn_zeros
from mevlib.scalar import sum_standard
from mevlib.scalar import sph_ptwise, sph_intgtd
from mevlib.scalar import cyl_ptwise_radial_terms
from mevlib.scalar import cyl_ptwise_axial_terms
from mevlib.scalar import cyl_intgtd_radial_terms, cyl_intgtd_axial_terms
from mevlib.scalar import psm_ptwise_series, psm_intgtd_series


# shapes

class Shape(ABC):
    pass

class Sphere(Shape):

    def __init__(self, R):
        self.R = R

    def ptwise(self, a2, r):
        return sph_ptwise(a2, self.R, r)

    def intgtd(self, a2, precision=None):
        # precision included for consistency; not needed here
        return sph_intgtd(a2, self.R)

class Cylinder(Shape):

    def __init__(self, R, H):
        self.R, self.H, self.zero_cache = R, H, None

    def ptwise_radial_terms(self, a2, trunc, r, z):
        return cyl_ptwise_radial_terms(a2, self.R, self.H, trunc, r, z)

    def ptwise_axial_terms(self, a2, trunc, r, z):
        if self.zero_cache is None or len(self.zero_cache) < trunc:
            self.zero_cache = jn_zeros(0, trunc)
        return cyl_ptwise_axial_terms(
            a2, self.R, self.H, trunc, r, z, zero_cache=self.zero_cache
        )

    def ptwise_radial(self, a2, ntrunc, r, z):
        return sum_standard(self.ptwise_radial_terms(a2, ntrunc, r, z))

    def ptwise_axial(self, a2, ktrunc, r, z):
        return sum_standard(self.ptwise_axial_terms(a2, ktrunc, r, z))

    def ptwise(self, a2, precision, r, z):
        return sum(map(sum_standard, [
            self.ptwise_radial_terms(a2, precision['ntrunc'], r, z),
            self.ptwise_axial_terms(a2, precision['ktrunc'], r, z)
        ]))

    def intgtd_radial_terms(self, a2, trunc, bratcutoff=None):
        return cyl_intgtd_radial_terms(
            a2, self.R, self.H, trunc, bratcutoff=bratcutoff
        )

    def intgtd_axial_terms(self, a2, trunc):
        if self.zero_cache is None or len(self.zero_cache) < trunc:
            self.zero_cache = jn_zeros(0, trunc)
        return cyl_intgtd_axial_terms(
            a2, self.R, self.H, trunc, zero_cache=self.zero_cache
        )

    def intgtd(self, a2, precision):
        ntrunc, ktrunc = precision['ntrunc'], precision['ktrunc']
        return sum(map(sum_standard, [
            self.intgtd_radial_terms(a2, ntrunc),
            self.intgtd_axial_terms(a2, ktrunc)
        ]))

class Prism(Shape):

    def __init__(self, Lx, Ly, Lz):
        self.Lx, self.Ly, self.Lz = Lx, Ly, Lz

    @staticmethod
    def split_prec(precision):
        if 'trunc' in precision and not any([
            k in precision.keys() for k in ['xtrunc', 'ytrunc', 'ztrunc']
        ]):
            return precision['trunc'], precision['trunc'], precision['trunc']
        else:
            return (precision[k] for k in ['xtrunc', 'ytrunc', 'ztrunc'])

    def ptwise(self, a2, precision, x, y, z):
        return psm_ptwise_series(
            a2, self.Lx, self.Ly, self.Lz, self.split_prec(precision), x, y, z
        )

    def intgtd(self, a2, precision):
        return psm_intgtd_series(
            a2, self.Lx, self.Ly, self.Lz, *self.split_prec(precision)
        )


