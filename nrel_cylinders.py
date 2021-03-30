

import numpy as np
from scipy.special import j0, j1, i0e, i1e, jn_zeros


# solution pointwise

def u1(a, R, H, trunc, r, z):
    cum = 0.0
    for n in range(trunc-1, -1, -1):
        gn = np.sqrt(a**2 + (np.pi / H * (2 * n + 1))**2)
        cum += (
            4.0 / (np.pi * (2 * n + 1))
            * np.sin(np.pi * (2 * n + 1) * z / H)
            * np.exp(gn * (r - R)) * i0e(gn * r) / i0e(gn * R)
        )
    return cum

def u2(a, R, H, trunc, r, z, zero_cache=None):
    if zero_cache is not None and len(zero_cache) >= trunc:
        alpha = zero_cache
    else:
        alpha = jn_zeros(0, trunc)
    cum = 0.0
    for k in range(trunc, 0, -1):
        lk = np.sqrt(a**2 + (alpha[k-1] / R)**2)
        cum += (
            np.cosh(lk * (z - H / 2)) / np.cosh(lk * H / 2)
            * j0(alpha[k-1] / R * r) / (j1(alpha[k-1]) * alpha[k-1])
        )
    return 2.0 * cum

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
        * np.exp(gn(n) * (r - R)) * i0e(gn(n) * r) / i0e(gn(n) * R)
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

def u1residual(a, R, H, ntrunc, r, z):
    return (
        + u1laplacian(a, R, H, ntrunc, r, z)
        - a**2 * u1(a, R, H, ntrunc, r, z)
    )

def u2residual(a, R, H, ktrunc, r, z):
    return (
        + u2laplacian(a, R, H, ktrunc, r, z)
        - a**2 * u2(a, R, H, ktrunc, r, z)
    )

def uresidual(a, R, H, ntrunc, ktrunc, r, z):
    return (
        + u1residual(a, R, H, ntrunc, r, z)
        + u2residual(a, R, H, ktrunc, r, z)
    )


# solution integrated over domain

def gn(a, H, n):
    return np.sqrt(a**2 + (np.pi * (2 * n + 1) / H)**2)

def ln(a, R, alphak):
    return np.sqrt(a**2 + (alphak / R)**2)

def s1(a, R, H, trunc, bratcutoff=None):
    ns = np.array([n for n in range(trunc-1, -1, -1)])
    gns = gn(a, H, ns)
    # the ratio I1/I0 tends to 1 (inf/inf) very quickly
    #TODO instead of computing then checking, just compute what is needed
    brat = i1e(gns * R) / i0e(gns * R)
    if bratcutoff is not None:
        brat = np.where(gns * R < bratcutoff, brat, 1.0)
    return sum(16.0 / (np.pi**2 * R * gns * (2 * ns + 1)**2) * brat)

def s2(a, R, H, trunc, zero_cache=None):
    if zero_cache is not None and len(zero_cache) >= trunc:
        alpha = zero_cache
    else:
        alpha = jn_zeros(0, trunc)
    lnks = ln(a, R, alpha[::-1])
    return 8.0 * sum(np.tanh(lnks * H / 2) / (H * lnks * alpha[::-1]**2))

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


# reaction mechanisms

class Mechanism(object):

    gasconst = 8.31446261815324e-3  # kJ per mol Kelvin

    def Di(self, T):
        return (
            48.5 * self.Dpore * np.sqrt(T / self.Ms) * self.epspore / self.tau
        )

    def kij(self, T):
        return (
            self.k0s
            * np.exp(-self.Eas / self.gasconst * (1 / T - 1 / self.T0))
        )

    def makematrix(self, T, L):
        kijs, Dis = self.kij(T), self.Di(T)
        phiij2s = kijs * L**2 / (np.tile(Dis.reshape((-1,1)),6))
        B = np.diag(np.sum(phiij2s, axis=0)) - phiij2s
        return B


class CEJMechanism(Mechanism):

    k0s = np.array([            # 1 / s
        [0.000, 0.000, 0.000, 0.000, 0.000, 0.000],
        [1.413, 0.000, 0.000, 0.000, 0.000, 0.000],
        [4.337, 0.229, 0.000, 0.000, 0.000, 0.000],
        [1.163, 0.161, 0.128, 0.000, 0.000, 0.000],
        [0.114, 0.041, 0.030, 0.000, 0.000, 0.000],
        [0.386, 0.137, 0.103, 0.000, 0.000, 0.000]
    ])

    Eas = np.array([            # kJ / mol
        [ 0.0,  0.0,  0.0, 0.0, 0.0, 0.0],      # noqa E201
        [47.6,  0.0,  0.0, 0.0, 0.0, 0.0],
        [43.4, 54.1,  0.0, 0.0, 0.0, 0.0],
        [38.5, 62.9, 80.5, 0.0, 0.0, 0.0],
        [30.2, 66.7, 85.2, 0.0, 0.0, 0.0],
        [30.0, 65.0, 77.3, 0.0, 0.0, 0.0]
    ])

    Ms = np.array([             # g / mol
                444.0, 230.0, 115.0, 52.0, 16.0, 400.0
                ])

    T0 = 773                        # Kelvin
    Dpore = 2.0                     # nm
    epspore = 0.319                 # nondim
    tau = 7.0                       # nondim


class FakeReversibleMechanism(Mechanism):

    Ms = CEJMechanism.Ms
    T0 = CEJMechanism.T0
    Dpore = CEJMechanism.Dpore
    epspore = CEJMechanism.epspore
    tau = CEJMechanism.tau

    def __init__(self):
        self.rand = np.random.rand(*CEJMechanism.k0s.shape)
        self.k0s = CEJMechanism.k0s + 0.5 * self.rand * CEJMechanism.k0s.T
        self.Eas = CEJMechanism.Eas + 2.0 * self.rand * CEJMechanism.Eas.T


