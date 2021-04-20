#!/usr/bin/env python3


from random import random

import numpy as np
from matplotlib import pyplot as plt

from nrel_cylinders import u, uresidual, u1residual, u2residual


plt.rc('text', usetex=True)
plt.rc('axes', labelsize=12)


def checkpoint(a, R, H, ntrunc, ktrunc, dr, dz, r, z):
    return (
        - (
            (r + 0.5 * dr) * (
                u(a, R, H, ntrunc, ktrunc, r + dr, z)
                -
                u(a, R, H, ntrunc, ktrunc, r, z)
            )
            -
            (r - 0.5 * dr) * (
                u(a, R, H, ntrunc, ktrunc, r, z)
                -
                u(a, R, H, ntrunc, ktrunc, r - dr, z)
            )
        ) / (2.0 * dr**2)
        - (
            + u(a, R, H, ntrunc, ktrunc, r, z + dz)
            - 2 * u(a, R, H, ntrunc, ktrunc, r, z)
            + u(a, R, H, ntrunc, ktrunc, r, z - dz)
        ) / dz**2
        +
        a**2 * u(a, R, H, ntrunc, ktrunc, r, z)
    )


a = 00.8
R = 04.0
H = 15.0


numcheckpts = 6000
chkntrunc, chkktrunc = 24, 24
chkdr, chkdz = 1e-4 * R, 1e-4 * H


fig6, ax6 = plt.subplots(1, 1, figsize=(3.0,3.0))
chkrs = [R * random() for i in range(numcheckpts)]
chkzs = [H * random() for i in range(numcheckpts)]
vals = [
    abs(checkpoint(a, R, H, chkntrunc, chkktrunc, chkdr, chkdz, r, z))
    for (r, z) in zip(chkrs, chkzs)
]
sc = ax6.scatter(chkrs, chkzs, c=vals, s=3.6)
fig6.colorbar(sc, ax=ax6)


truncs = [6, 12, 24, 48]

rs, zs = np.linspace(0.0, 1.05 * R, 120), np.linspace(-0.05 * H, 1.05 * H, 320)
rmesh, zmesh = np.meshgrid(rs, zs)
fig1, axs1 = plt.subplots(len(truncs), len(truncs))
fig2, axs2 = plt.subplots(len(truncs), len(truncs))
for fig, axs, xlims, ylims, zlim, in [
    (fig1, axs1, None, None, None),
    (fig2, axs2, (0.0, R), (0.0, H), None)
     ]:
    for i, trunc1 in enumerate(truncs):
        for j, trunc2 in enumerate(truncs):
            vals = uresidual(a, R, H, trunc1, trunc2, rmesh, zmesh)
            if zlim is not None:
                minv, maxv = 0.5 - zlim, 0.5 + zlim
            elif xlims is not None and ylims is not None:
                idxs = np.logical_and(
                    np.logical_and(xlims[0] < rmesh, rmesh < xlims[1]),
                    np.logical_and(ylims[0] < zmesh, zmesh < ylims[1])
                )
                minv, maxv = np.min(vals[idxs]), np.max(vals[idxs])
            else:
                minv, maxv = -np.inf, np.inf
            vals = np.minimum(np.maximum(vals, minv), maxv)
            cf = axs[i, j].contourf(rmesh, zmesh, vals)
            axs[i, j].set_xlim(xlims)
            axs[i, j].set_ylim(ylims)
            fig.colorbar(cf, ax=axs[i, j])
            axs[i, j].set_title("{}, {}".format(trunc1, trunc2))
            axs[i, j].set_xticklabels([])
            axs[i, j].set_yticklabels([])
#           axs[i, j].set_xlabel("r")
#           axs[i, j].set_ylabel("z")


maxterms, meshsize = 128, 400
rs, zs = np.linspace(0.0, R, meshsize+1), np.linspace(0.0, H, meshsize+1)
rs = 0.5 * (rs[1:] + rs[:-1])
zs = 0.5 * (zs[1:] + zs[:-1])
rmesh, zmesh = np.meshgrid(rs, zs)
fig3, axs3 = plt.subplots(1, 2)
for i, f in enumerate([u1residual, u2residual]):
    resids = np.array([
        np.sum(np.abs(f(a, R, H, n+1, rmesh, zmesh)))
        for n in range(maxterms)
    ]) / meshsize**2
    axs3[i].plot(np.arange(maxterms) + 1, resids, 'k.')
    axs3[i].set_xlabel(
        r"Number of Terms (\( {} \))".format("N" if i == 0 else "K")
    )
    axs3[i].set_ylabel("Mean Residual")
    axs3[i].grid()
    axs3[i].set_title(
        "Integrated Error in {} Sum".format("First" if i == 0 else "Second")
    )


if False:
    for i, fig in enumerate([fig1, fig2, fig3, fig6]):
        fig.tight_layout()
        for ext in ['svg', 'pdf']:
            fig.savefig("img/ptcheck-fig{}.{}".format(i+1,ext))
else:
    plt.show()


