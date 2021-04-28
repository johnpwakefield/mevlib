#!/usr/bin/env python3


import numpy as np

from matplotlib import pyplot as plt

from mevlib.scalar import psm_ptwise_xdir, psm_ptwise_series
from mevlib.scalar import psm_spec_coeffs, psm_spec_eval
from mevlib.scalar import psm_diff_coeffs, psm_diff_eval
from mevlib.scalar import psm_diff_axes
from mevlib.scripts import imgpath, showfigs


plt.rc('text', usetex=True)
plt.rc('axes', labelsize=12)


a2 = 1.2**2
Lx, Ly, Lz = 0.3, 0.4, 1.7
#series_trunc, spec_trunc, diff_trunc = 32, 16, 2
series_trunc, spec_trunc, diff_trunc = 8, 4, 32
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
        np.vectorize(psm_ptwise_xdir)(a2, L1, L2, L3, series_trunc, c1, c2, c3)
    )

fig5, ax5 = plt.subplots(1, 3, figsize=(3.0,3.0))
vals = np.vectorize(psm_ptwise_series)(
    a2, Lx, Ly, Lz, series_trunc, xmesh, ymesh, z
)
vals[np.logical_and(
    np.minimum(xmesh, Lx - xmesh) < 0.05 * Lx,
    np.minimum(ymesh, Ly - ymesh) < 0.05 * Ly
)] = np.nan
cf = ax5[0].contourf(xmesh, ymesh, np.minimum(np.maximum(vals, -0.2), 1.2))
fig5.colorbar(cf, ax=ax5[0])
ax5[0].set_xlabel(r"\( x \)")
ax5[0].set_ylabel(r"\( y \)")
ax5[0].set_title("Series with {} Terms".format(3 * series_trunc**2))
coeffs = psm_spec_coeffs(a2, Lx, Ly, Lz, spec_trunc)
vals = np.empty(xmesh.shape)
for i in range(xmesh.shape[0]):
    for j in range(ymesh.shape[1]):
        vals[i, j] = psm_spec_eval(coeffs, xmesh[i, j], ymesh[i, j], z)
cf = ax5[1].contourf(xmesh, ymesh, vals)
fig5.colorbar(cf, ax=ax5[1])
ax5[1].set_xlabel(r"\( x \)")
ax5[1].set_ylabel(r"\( y \)")
ax5[1].set_title("Pseudospectral with {} Terms".format(spec_trunc**3))
coeffs = psm_diff_coeffs(a2, Lx, Ly, Lz, diff_trunc)
#als = np.empty(xmesh.shape)
#or i in range(xmesh.shape[0]):
#   for j in range(ymesh.shape[1]):
#       vals[i, j] = psm_diff_eval(coeffs, xmesh[i, j], ymesh[i, j], z)
#f = ax5[2].contourf(xmesh, ymesh, vals)
#TODO DEBUG
xs, ys, zs = psm_diff_axes(Lx, Ly, Lz, diff_trunc, diff_trunc, diff_trunc)
cf = ax5[2].contourf(*np.meshgrid(xs, ys), coeffs[-1][:, :, 3])
#TODO DEBUG
fig5.colorbar(cf, ax=ax5[2])
ax5[2].set_xlabel(r"\( x \)")
ax5[2].set_ylabel(r"\( y \)")
ax5[2].set_title("Finite Differences with {} Terms".format(spec_trunc**3))


if showfigs():
    plt.show()
else:
    for i, fig in enumerate([fig4, fig5]):
        for ext in ['pdf', 'svg']:
            fig.savefig(imgpath("psm_solvis-{}.{}".format(i+1, ext)))


