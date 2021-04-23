

import numpy as np


# helper function

def beta(ms, ns, a2, L1, L2, L3):
    return np.sqrt(
        (np.pi * (2 * ms + 1) * (L1 / L2))**2 +
        (np.pi * (2 * ns + 1) * (L1 / L3))**2 +
        L1**2 * a2
    )


# solution pointwise

def upartscalar(a2, Lx, Ly, Lz, trunc, x, y, z):
    ms, ns = np.meshgrid(np.arange(trunc)[::-1], np.arange(trunc)[::-1])
    return 16.0 / np.pi**2 * np.sum(
        np.sin(np.pi * (2 * ms + 1) * y / Ly) *
        np.sin(np.pi * (2 * ns + 1) * z / Lz) *
        np.cosh(beta(ms, ns, a2, Lx, Ly, Lz) * (x / Lx - 0.5)) /
        np.cosh(beta(ms, ns, a2, Lx, Ly, Lz) * 0.5) /
        ((2 * ms + 1) * (2 * ns + 1))
    )
upart = np.vectorize(upartscalar)
def uscalar(a2, Lx, Ly, Lz, truncs, x, y, z):
    if np.isscalar(truncs):
        truncs = truncs, truncs, truncs
    else:
        assert(len(truncs) == 3)
    return sum([
        upartscalar(a2, L1, L2, L3, trunc, c1, c2, c3)
        for c1, c2, c3, L1, L2, L3, trunc in [
            (x, y, z, Lx, Ly, Lz, truncs[0]),
            (y, z, x, Ly, Lz, Lx, truncs[1]),
            (z, x, y, Lz, Lx, Ly, truncs[2])
        ]
    ])
u = np.vectorize(uscalar)


# solution integrated over domain

def s(a2, Lx, Ly, Lz, trunc):
    ms, ns = np.meshgrid(np.arange(trunc), np.arange(trunc))
    return 32.0 / np.pi**4 * sum([
        np.sum(
            np.tanh(beta(ms, ns, a2, L1, L2, L3) * 0.5) /
            (beta(ms, ns, a2, L1, L2, L3) * (2 * ms + 1)**2 * (2 * ns + 1)**2)
        )
        for L1, L2, L3 in [(Lx, Ly, Lz), (Ly, Lz, Lx), (Lz, Lx, Ly)]
    ])


