#!/usr/bin/env python3


import numpy as np

from nrel_prisms import s


# raw data

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
gasconst = 8.31446261815324e-3  # kJ per mol Kelvin
Dpore = 2.0                     # nm
epspore = 0.319                 # nondim
tau = 7.0                       # nondim
L = 60.0                        # micrometers
trunc = 64                      # truncation

cases = [       # Lx, Ly, Lz, T, C1, C2, C3, C4, C5, C6
    (120.0, 120.0, 230.0, 600.0, (0.9, 0.8, 0.7, 0.6, 0.2, 0.2)),
    (120.0, 140.0, 330.0, 800.0, (0.9, 0.8, 0.7, 0.6, 0.2, 0.2)),
    (120.0,  80.0, 280.0, 600.0, (0.3, 0.4, 0.5, 0.6, 0.4, 0.3))
]


# functions

def Di(T):
    return 48.5 * Dpore * np.sqrt(T / Ms) * epspore / tau

def kij(k0, Ea, T):
    return k0 * np.exp(-Ea / gasconst * (1 / T - 1 / T0))

def makematrix(T):
    kijs, Dis = kij(k0s, Eas, T), Di(T)
    phiij2s = kijs * L**2 / (np.tile(Dis.reshape((-1,1)),6))
    B = np.diag(np.sum(phiij2s, axis=0)) - phiij2s
    return B

def factor(B):
    EVs = np.eye(6)
    for k in range(6):
        for l in range(k+1,6):
            if B[k,k] != B[l,l]:
                EVs[l,k] = np.dot(EVs[k:l,k], B[l,k:l]) / (B[k,k] - B[l,l])
    return EVs

def backsub(EVs, x):
    y = np.empty((6,))
    y[0] = x[0] / EVs[0,0]
    for k in range(1, 6):
        y[k] = (x[k] - np.dot(EVs[k,:k], y[:k])) / EVs[k,k]
    return y

def eff_factor_general(Lx, Ly, Lz, T, C):
    lams, EVs = np.linalg.eig(makematrix(T))
    u = np.linalg.solve(EVs, np.array(C).reshape((-1,1))).reshape((-1,))
    mult = np.array([
        1.0 if lam == 0.0 else
        s(lam, Lx, Ly, Lz, trunc)
        for lam in lams
    ])
    res = np.dot(EVs, (mult * u).reshape((-1,1))).reshape((-1,)) / np.array(C)
    if np.any(np.logical_not(np.isfinite(res))):
        print("Non finite values in result.")
        print(res)
    return res

def eff_factor_fast(Lx, Ly, Lz, T, C):
    B = makematrix(T)
    EVs = factor(B)
    ufree = backsub(EVs, np.array(C)).reshape((-1,))
    ubar = ufree * np.array([
        1.0 if a2 == 0.0 else
        s(a2, Lx, Ly, Lz, trunc)
        for a2 in np.diag(B)
    ])
    etas = (
        np.dot(EVs, ubar.reshape((-1,1))).reshape((-1,))
        / np.array(C).reshape((-1,))
    )
    return etas


# compute for varied cases

def stringifyetas(etas):
    return ["{:.4e}".format(eta) for eta in etas]

for efactor, name in [
    (eff_factor_general, "General:"),
    (eff_factor_fast, "Fast:")
]:
    print(name)
    print(" & ".join(list(map(
        lambda s: r"\( {} \)".format(s),
        [
            "Lx", "Ly", "Lz", "T",
            *[r"C_{}".format(i+1) for i in range(6)],
            *[r"\eta_{}".format(i+1) for i in range(6)]
        ]
    ))) + r" \\\hline\hline")
    for Lx, Ly, Lz, T, C in cases:
        print(
            " & ".join(list(map(str,
                                [
                                    Lx, Ly, Lz, T, *C,
                                    *stringifyetas(efactor(Lx, Ly, Lz, T, C))
                                ]
                                )))
            + r" \\\hline"
        )


