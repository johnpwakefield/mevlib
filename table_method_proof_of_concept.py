#!/usr/bin/env python3


# we need to
#   1.  measure error in the reconstructed matrix
#   2.  determine how many iterations are required with a nearby inverse used
#           as a preconditioner


import numpy as np

import matplotlib.pyplot as plt

from nrel_cylinders import CEJMechanism, FakeReversibleMechanism


mechs = [
    (CEJMechanism(), "CEJ"),
    (FakeReversibleMechanism(), "Fake Reversible"),
]
T0, T1 = 700.0, 710.0
Ts = np.linspace(T0, T1, 12)
dTs = [T - T0 for T in Ts]
sortol = 1e-6
L = 60.0
Cs = [np.random.rand(6,1) for i in range(4)]


def backsub(A, x):
    y = np.empty_like(x)
    y[0] = x[0] / A[0,0]
    for k in range(1, 6):
        y[k] = (x[k] - np.dot(A[k,:k], y[:k])) / A[k,k]
    return y

def SORiters(A, b, x, omega, tol, Nmax=10000):
    b, x = b.reshape((-1,1)), x.reshape((-1,1))
    D = np.diag(np.diag(A))
    L, U = np.tril(A) - D, np.triu(A) - D
    for n in range(Nmax):
        if np.linalg.norm(np.dot(A, x) - b) < tol:
            return n
        x = backsub(
            D + omega * L,
            omega * b - np.dot(omega * U + (omega - 1) * D, x)
        )
    return Nmax


fig, axs = plt.subplots(2, 3, figsize=(6.0, 6.0))

for i, (mech, mechname) in enumerate(mechs):

    A0, A1 = mech.makematrix(T0, L), mech.makematrix(T1, L)
    Am = (A1 - A0) / (T1 - T0)
    L0, R0 = np.linalg.eig(A0)
    L1, R1 = np.linalg.eig(mech.makematrix(T1, L))
    L0, L1 = L0[::-1], L1[::-1]
    R0, R1 = R0[:,::-1], R1[:,::-1]
    Lm, Rm = (L1 - L0) / (T1 - T0), (R1 - R0) / (T1 - T0)
    R0inv = np.linalg.inv(R0)
    jaciter = np.eye(6) - np.dot(np.diag(np.diag(R0)**(-1)), R0)
    jacspec = np.linalg.eig(jaciter)[0]
    specrad = np.max(np.abs(jacspec))       # 0, which seems obvious now
#   print(specrad)
#   omega = 1.0 + (specrad / (1.0 + np.sqrt(1.0 - specrad**2)))**2
    omega = 1.25

    # error in reconstructed matrix
    eigerrs = list(map(np.linalg.norm, [
        np.dot(
            np.dot(
                R0 + (T - T0) * Rm,
                np.diag(L0 + (T - T0) * Lm)
            ),
            np.linalg.inv(R0 + (T - T0) * Rm)
        )
        - mech.makematrix(T, L)
        for T in Ts
    ]
    ))
    linerrs = list(map(np.linalg.norm, [
        A0 + (T - T0) * Am - mech.makematrix(T, L) for T in Ts
    ]
    ))
    axs[i,0].plot(dTs, linerrs, 'k.')
    axs[i,0].plot(dTs, eigerrs, 'b.')
    axs[i,0].set_xlabel(r"\( \Delta t \)")
    axs[i,0].set_ylabel(r"Error")
    axs[i,0].set_title("Matrix Error for {}".format(mechname))

    # condition number of preconditioned matrix
    axs[i,1].plot(dTs, [np.linalg.cond(R0 + (T - T0) * Rm) for T in Ts], 'k.')
    axs[i,1].plot(
        dTs, [np.linalg.cond(np.dot(R0inv, R0 + (T - T0) * Rm)) for T in Ts],
        'b.'
    )
    axs[i,1].set_xlabel(r"\( \Delta T \)")
    axs[i,1].set_ylabel("Condition Number")
    axs[i,1].set_title("Condition Number of Preconditioned Matrix")

    # number of iterations required
    for j, C in enumerate(Cs):
        axs[i,2].plot(
            dTs, [
                SORiters(R0 + (T - T0) * Rm, C, C, omega, sortol)
                for T in Ts
            ], 'r.', label="Unconditioned"
        )
        axs[i,2].plot(
            dTs, [
                SORiters(
                    np.dot(R0inv,R0 + (T - T0) * Rm), np.dot(R0inv,C),
                    np.dot(R0inv, C), omega, sortol
                )
                for T in Ts
            ], 'b.', label="Preconditioned"
        )
    axs[i,2].set_xlabel(r"\( \Delta T \)")
    axs[i,2].set_ylabel("SOR Iterations")


plt.show()


