#!/usr/bin/env python3


import numpy as np

from matplotlib import pyplot as plt

from lib_scalar import psm_ptwise, psm_ptwise_xdir


plt.rc('text', usetex=True)
plt.rc('axes', labelsize=12)


a2 = 1.3
Lx, Ly, Lz = 0.8, 1.8, 1.2
trunc = 70
xs = np.linspace(0.0, Lx, 80)
ys = np.linspace(0.0, Ly, 81)
z = Lz / 3
xmesh, ymesh = np.meshgrid(xs, ys)


fig4, axs4 = plt.subplots(1, 3, figsize=(9.0, 3.0))
for i, (c1, c2, c3, L1, L2, L3) in enumerate([
    (xmesh, ymesh, z, Lx, Ly, Lz),
    (ymesh, z, xmesh, Ly, Lz, Lx),
    (z, xmesh, ymesh, Lz, Lx, Ly)
]):
    axs4[i].contourf(
        xmesh, ymesh,
        np.vectorize(psm_ptwise_xdir)(a2, L1, L2, L3, trunc, c1, c2, c3)
    )

fig5, ax5 = plt.subplots(1, 1, figsize=(3.0,3.0))
vals = np.vectorize(psm_ptwise)(a2, Lx, Ly, Lz, trunc, xmesh, ymesh, z)
cf = ax5.contourf(xmesh, ymesh, np.minimum(np.maximum(vals, -0.2), 1.2))
fig5.colorbar(cf, ax=ax5)
ax5.set_xlabel(r"\( x \)")
ax5.set_ylabel(r"\( y \)")


plt.show()


