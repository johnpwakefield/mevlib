"""

Contains functions that yield pointwise and integrated solutions to the scalar
problem for a variety of shapes.

"""

import numpy as np
from scipy.special import j0, j1, i0e, i1e, jn_zeros


# general functions

def sum_standard(terms):
    return np.sum(terms[::-1])

def sum_aitken(terms):
    anp1, anp2, anp3 = terms[0], terms[1], terms[2]
    cml = anp1 - anp2**2 / (anp3 - anp2)
    for n in range(len(terms)-3):
        anp1, anp2, anp3 = anp2, anp3, terms[n+3]
        cml += anp1 - anp2**2 / (anp3 - anp2) + anp1**2 / (anp2 - anp1)
    return cml


# cylinder functions

def cyl_gn(a2, H, n):
    return np.sqrt(a2 + (np.pi * (2 * n + 1) / H)**2)

def cyl_ln(a2, R, alphak):
    return np.sqrt(a2 + (alphak / R)**2)

def cyl_ptwise_radial_terms(a2, R, H, trunc, r, z):
    terms = np.empty((trunc,))
    for n in range(trunc-1, -1, -1):
        gn = cyl_gn(a2, H, n)
        terms[n] = (
            4.0 / (np.pi * (2 * n + 1))
            * np.sin(np.pi * (2 * n + 1) * z / H)
            * np.exp(gn * (r - R)) * i0e(gn * r) / i0e(gn * R)
        )
    return terms

def cyl_ptwise_axial_terms(a2, R, H, trunc, r, z, zero_cache=None):
    if zero_cache is not None and len(zero_cache) >= trunc:
        alpha = zero_cache
    else:
        alpha = jn_zeros(0, trunc)
    terms = np.empty((trunc,))
    for k in range(trunc, 0, -1):
        lk = cyl_ln(a2, R, alpha[k-1])
        terms[k-1] = (
            2.0 * np.cosh(lk * (z - H / 2)) / np.cosh(lk * H / 2)
            * j0(alpha[k-1] / R * r) / (j1(alpha[k-1]) * alpha[k-1])
        )
    return terms

def cyl_ptwise_radial(a2, R, H, ntrunc, r, z):
    return sum_standard(cyl_ptwise_radial_terms(a2, R, H, ntrunc, r, z))

def cyl_ptwise_axial(a2, R, H, ktrunc, r, z, zero_cache=None):
    return sum_standard(cyl_ptwise_axial_terms(
        a2, R, H, ktrunc, r, z, zero_cache=zero_cache
    ))

def cyl_ptwise(a2, R, H, ntrunc, ktrunc, r, z, zero_cache=None):
    return sum(map(sum_standard, [
        cyl_ptwise_radial_terms(a2, R, H, ntrunc, r, z),
        cyl_ptwise_axial_terms(a2, R, H, ktrunc, r, z, zero_cache)
    ]))

cyl_ptwise_radial_vect = np.vectorize(cyl_ptwise_radial)
cyl_ptwise_axial_vect = np.vectorize(cyl_ptwise_axial)
cyl_ptwise_vect = np.vectorize(cyl_ptwise)

def cyl_intgtd_radial_terms(a2, R, H, trunc, bratcutoff=None):
    ns = np.arange(trunc)
    gns = cyl_gn(a2, H, ns)
    # the ratio I1/I0 tends to 1 (inf/inf) very quickly
    #TODO instead of computing then checking, just compute what is needed
    brat = i1e(gns * R) / i0e(gns * R)
    if bratcutoff is not None:
        brat = np.where(gns * R < bratcutoff, brat, 1.0)
    return 16.0 / (np.pi**2 * R * gns * (2 * ns + 1)**2) * brat

def cyl_intgtd_axial_terms(a2, R, H, trunc, zero_cache=None):
    if zero_cache is not None and len(zero_cache) >= trunc:
        alpha = zero_cache
    else:
        alpha = jn_zeros(0, trunc)
    lnks = cyl_ln(a2, R, alpha)
    return 8.0 * np.tanh(lnks * H / 2) / (H * lnks * alpha**2)

def cyl_intgtd(a2, R, H, ntrunc, ktrunc, bratcutoff=None, zero_cache=None):
    return sum(map(sum_standard, [
        cyl_intgtd_radial_terms(a2, R, H, ntrunc, bratcutoff=bratcutoff),
        cyl_intgtd_axial_terms(a2, R, H, ktrunc, zero_cache=zero_cache)
    ]))


# prism functions

def psm_beta(a2, L1, L2, L3, ms, ns):
    return np.sqrt(
        (np.pi * (2 * ms + 1) * (L1 / L2))**2 +
        (np.pi * (2 * ns + 1) * (L1 / L3))**2 +
        L1**2 * a2
    )

def psm_ptwise_xdir(a2, Lx, Ly, Lz, trunc, x, y, z):
    ms, ns = np.meshgrid(np.arange(trunc)[::-1], np.arange(trunc)[::-1])
    betas = psm_beta(a2, Lx, Ly, Lz, ms, ns)
    return 16.0 / np.pi**2 * np.sum(
        np.sin(np.pi * (2 * ms + 1) * y / Ly) *
        np.sin(np.pi * (2 * ns + 1) * z / Lz) *
        (
            #TODO add truncation for large values
            np.cosh(betas * (x / Lx - 0.5)) / np.cosh(betas * 0.5)
        ) / ((2 * ms + 1) * (2 * ns + 1))
    )

def psm_ptwise(a2, Lx, Ly, Lz, truncs, x, y, z):
    if np.isscalar(truncs):
        truncs = truncs, truncs, truncs
    else:
        assert(len(truncs) == 3)
    return sum([
        psm_ptwise_xdir(a2, L1, L2, L3, trunc, c1, c2, c3)
        for c1, c2, c3, L1, L2, L3, trunc in [
            (x, y, z, Lx, Ly, Lz, truncs[0]),
            (y, z, x, Ly, Lz, Lx, truncs[1]),
            (z, x, y, Lz, Lx, Ly, truncs[2])
        ]
    ])

def psm_intgtd_xdir(a2, Lx, Ly, Lz, trunc):
    ms, ns = np.meshgrid(np.arange(trunc)[::-1], np.arange(trunc)[::-1])
    betas = psm_beta(a2, Lx, Ly, Lx, ms, ns)
    return 32.0 / np.pi**4 * np.sum(
        np.tanh(betas * 0.5) / (betas * (2 * ms + 1)**2 * (2 * ns + 1)**2)
    )

def psm_intgtd(a2, Lx, Ly, Lz, truncs):
    if np.isscalar(truncs):
        truncs = truncs, truncs, truncs
    else:
        assert(len(truncs) == 3)
    return sum([
        psm_intgtd_xdir(a2, L1, L2, L3, trunc)
        for L1, L2, L3, trunc in [
            (Lx, Ly, Lz, truncs[0]),
            (Ly, Lz, Lx, truncs[1]),
            (Lz, Lx, Ly, truncs[2])
        ]
    ])


# sphere functions

def sph_ptwise(a2, R, r):
    return (R / r) * np.sinh(np.sqrt(a2) * r) / np.sinh(np.sqrt(a2) * R)

def sph_intgtd(a2, R):
    a = np.sqrt(a2)
    return 3.0 / (R * a) * (np.coth(R * a) - (R * a)**(-1))


