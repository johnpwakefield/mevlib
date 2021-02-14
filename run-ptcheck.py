#!/usr/bin/env python3


from random import random

from matplotlib import pyplot as plt

from nrel_cylinders import u


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


numcheckpts = 2000
chkntrunc, chkktrunc = 24, 24
chkdr, chkdz = 0.01 * R, 0.01 * H


fig6, ax6 = plt.subplots(1, 1, figsize=(3.0,3.0))
chkrs = [R * random() for i in range(numcheckpts)]
chkzs = [H * random() for i in range(numcheckpts)]
vals = [
    abs(checkpoint(a, R, H, chkntrunc, chkktrunc, chkdr, chkdz, r, z))
    for (r, z) in zip(chkrs, chkzs)
]
sc = ax6.scatter(chkrs, chkzs, c=vals, s=3.6)
fig6.colorbar(sc, ax=ax6)


if True:
    fig6.tight_layout()
    for ext in ['svg', 'pdf']:
        fig6.savefig("img/run-ptcheck.{}".format(ext))
else:
    plt.show()


