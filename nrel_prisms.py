

import numpy as np


# helper function

def beta(ms, ns, a, L2, L3):
    return np.sqrt(
        (np.pi * (2 * ms + 1) / L2)**2 +
        (np.pi * (2 * ns + 1) / L3)**2 +
        a**2
    )


# solution pointwise

def u(a, Lx, Ly, Lz, trunc, x, y, z):
    ms, ns = np.meshgrid(np.arange(trunc), np.arange(trunc))
    return 16.0 / np.pi**2 * sum([
        np.sum(
            ((2 * ms + 1) * (2 * ns + 1))**(-1) *
            np.sin(np.pi * (2 * ms + 1) * c2 / L2) *
            np.sin(np.pi * (2 * ns + 1) * c3 / L3) *
            np.cosh(beta(ms, ns, a, L2, L3) * (c1 - L1 / 2)) /
            np.cosh(beta(ms, ns, a, L2, L3) * L1 / 2)
        )
        for c1, c2, c3, L1, L2, L3 in [
            (x, y, z, Lx, Ly, Lz), (y, z, x, Ly, Lz, Lx), (z, x, y, Lz, Lx, Ly)
        ]
    ])


# solution integrated over domain

def s(a, Lx, Ly, Lz, trunc):
    ms, ns = np.meshgrid(np.arange(trunc), np.arange(trunc))
    return 128.0 / np.pi**4 * sum([
        np.sum(
            np.tanh(beta(ms, ns, a, L2, L3) * L1 / 2) /
            (L1 * beta(ms, ns, a, L2, L3) * (2 * ms + 1)**2 * (2 * ns + 1)**2)
        )
        for L1, L2, L3 in [(Lx, Ly, Lz), (Ly, Lz, Lx), (Lz, Lx, Ly)]
    ])


