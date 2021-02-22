

import numpy as np
from scipy.special import j0, j1, i0, i1, jn_zeros


# solution pointwise

def u1(a, R, H, trunc, r, z):
    return sum([
        4.0 / (np.pi * (2 * n + 1))
        * np.sin(np.pi * (2 * n + 1) * z / H)
        * i0(np.sqrt(a**2 + (np.pi / H * (2 * n + 1))**2) * r)
        / i0(np.sqrt(a**2 + (np.pi / H * (2 * n + 1))**2) * R)
        for n in range(trunc-1, -1, -1)
    ])

def u2(a, R, H, trunc, r, z, zero_cache=None):
    if zero_cache is not None and len(zero_cache) >= trunc:
        alpha = zero_cache
    else:
        alpha = jn_zeros(0, trunc)
    return 2.0 * sum([
        np.cosh(np.sqrt(a**2 + (alpha[k-1] / R)**2) * (z - H / 2))
        / np.cosh(np.sqrt(a**2 + (alpha[k-1] / R)**2) * H / 2)
        * j0(alpha[k-1] / R * r) / (alpha[k-1] * j1(alpha[k-1]))
        for k in range(trunc, 0, -1)
    ])

def u(a, R, H, ntrunc, ktrunc, r, z, zero_cache=None):
    return (
        u1(a, R, H, ntrunc, r, z)
        +
        u2(a, R, H, ktrunc, r, z, zero_cache=zero_cache)
    )


# residual pointwise

def u1laplacian(a, R, H, trunc, r, z):
    def gn(n):
        return np.sqrt(a**2 + (np.pi / H * (2 * n + 1))**2)
    return sum([
        4.0 * np.pi * (2 * n + 1) * gn(n)**2 / H**2
        * np.sin(np.pi * (2 * n + 1) * z / H)
        * i0(gn(n) * r) / i0(gn(n) * R)
        for n in range(trunc-1, -1, -1)
    ])

def u2laplacian(a, R, H, trunc, r, z):
    alpha = jn_zeros(0, trunc)
    def lk(k):
        return np.sqrt(a**2 + (alpha[k-1] / R)**2)
    return - 2.0 / R**2 * sum([
        lk(k)**2 * alpha[k-1]
        * np.cosh(lk(k) * (z - H / 2)) / np.cosh(lk(k) * H / 2)
        * j0(alpha[k-1] / R * r) / j1(alpha[k-1])
        for k in range(trunc, 0, -1)
    ])

def uresidual(a, R, H, ntrunc, ktrunc, r, z):
    return (
        + u1laplacian(a, R, H, ntrunc, r, z)
        - a**2 * u1(a, R, H, ntrunc, r, z)
        + u2laplacian(a, R, H, ntrunc, r, z)
        - a**2 * u2(a, R, H, ktrunc, r, z)
    )


# solution integrated over domain

def gn(a, H, n):
    return np.sqrt(a**2 + (np.pi * (2 * n + 1) / H)**2)

def ln(a, R, alphak):
    return np.sqrt(a**2 + (alphak / R)**2)

def s1(a, R, H, trunc, bratcutoff=400.0):
    ns = np.array([n for n in range(trunc-1, -1, -1)])
    gns = gn(a, H, ns)
    # the ratio I1/I0 tends to 1 (inf/inf) very quickly
    brat = np.where(gns * R < bratcutoff, i1(gns * R) / i0(gns * R), 1.0)
    return sum(16.0 / (np.pi**2 * R * gns * (2 * ns + 1)**2) * brat)

def s2(a, R, H, trunc, zero_cache=None):
    if zero_cache is not None and len(zero_cache) >= trunc:
        alpha = zero_cache
    else:
        alpha = jn_zeros(0, trunc)
    lnks = ln(a, R, alpha)
    return 8.0 * sum(np.tanh(lnks * H / 2) / (H * lnks * alpha**2))

def s(a, R, H, ntrunc, ktrunc, zero_cache=None):
    s1val = s1(a, R, H, ntrunc)
    s2val = s2(a, R, H, ktrunc, zero_cache=zero_cache)
    if np.any(np.logical_not(np.isfinite(s1val))):
        print("Nonfinite values in s1.")
        print(s1val)
    if np.any(np.logical_not(np.isfinite(s2val))):
        print("Nonfinite values in s2.")
        print(s2val)
    return s1val + s2val


