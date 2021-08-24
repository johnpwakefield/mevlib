#!/usr/bin/env python3


import numpy as np

from matplotlib import pyplot as plt

from mevlib.scalar import psm_ptwise_xdir, psm_ptwise_series
from mevlib.scalar import psm_spec_coeffs, psm_spec_eval
from mevlib.scalar import psm_diff_coeffs, psm_diff_axes
from mevlib.scalar import psm_ptwise_diff2
from mevlib.options import imgpath, showfigs


plt.rc('text', usetex=True)
plt.rc('axes', labelsize=12)


def clip(arr):
    return np.minimum(1.0, np.maximum(0.0, arr))


a2 = 0.15**2
Lx, Ly, Lz = 0.3 * 200, 0.4 * 200, 1.7 * 200
series_trunc, spec_trunc, diff_trunc = 12, 4, 16
xs = np.linspace(0.0, Lx, 80)
ys = np.linspace(0.0, Ly, 81)
zs = np.linspace(0.0, Lz, 160)
z = Lz / 3


xmesh, zmesh = np.meshgrid(xs, zs)
slicefig, sliceax = plt.subplots(1, 1, figsize=(3.0, 4.0))
cf = sliceax.contourf(
    xmesh, zmesh, np.vectorize(psm_ptwise_series)(
        a2,
        Lx, Ly, Lz, series_trunc, series_trunc, series_trunc,
        xmesh, Ly / 3, zmesh
    )
)
sliceax.set_xlabel(r"\( x \)")
sliceax.set_ylabel(r"\( z \)")
slicefig.colorbar(cf, ax=sliceax)
slicefig.tight_layout()


xmesh, ymesh = np.meshgrid(xs, ys)
pcfig, pcaxs = plt.subplots(1, 3, figsize=(9.0, 3.0))
for i, (c1, c2, c3, L1, L2, L3) in enumerate([
    (xmesh, ymesh, z, Lx, Ly, Lz),
    (ymesh, z, xmesh, Ly, Lz, Lx),
    (z, xmesh, ymesh, Lz, Lx, Ly)
]):
    pcaxs[i].contourf(
        xmesh, ymesh, np.vectorize(psm_ptwise_xdir)(
            a2, L1, L2, L3, series_trunc, c1, c2, c3
        )
    )
    pcaxs[i].set_title(r"\( u^{}(x, y, z^*) \)".format(["x", "y", "z"][i]))
    pcaxs[i].set_xlabel(r"\( x \)")
    pcaxs[i].set_ylabel(r"\( y \)")
pcfig.tight_layout()

levels = np.linspace(0.0, 1.0, 20, endpoint=True)
errfig, erraxs = plt.subplots(1, 3, figsize=(9.0, 3.0))
vals = np.vectorize(psm_ptwise_series)(
    a2, Lx, Ly, Lz, series_trunc, series_trunc, series_trunc, xmesh, ymesh, z
)
cf = erraxs[0].contourf(
    xmesh, ymesh, vals, levels=levels
)
erraxs[0].set_xlabel(r"\( x \)")
erraxs[0].set_ylabel(r"\( y \)")
erraxs[0].set_title("Series with {} Terms".format(3 * series_trunc**2))
coeffs = psm_spec_coeffs(a2, Lx, Ly, Lz, spec_trunc)
vals = np.empty(xmesh.shape)
for i in range(xmesh.shape[0]):
    for j in range(ymesh.shape[1]):
        vals[i, j] = psm_spec_eval(coeffs, xmesh[i, j], ymesh[i, j], z)
cf = erraxs[1].contourf(xmesh, ymesh, vals, levels=levels)
erraxs[1].set_xlabel(r"\( x \)")
erraxs[1].set_ylabel(r"\( y \)")
erraxs[1].set_title("Pseudospectral ({} Points)".format(spec_trunc**3))
coeffs = psm_ptwise_diff2(a2, Lx, Ly, Lz, diff_trunc, diff_trunc, diff_trunc)
xs, ys, zs = psm_diff_axes(Lx, Ly, Lz, diff_trunc, diff_trunc, diff_trunc)
cf = erraxs[2].contourf(
    *np.meshgrid(xs, ys), coeffs[:, :, 3], levels=levels
)
errfig.colorbar(cf, ax=erraxs[2])
erraxs[2].set_xlabel(r"\( x \)")
erraxs[2].set_ylabel(r"\( y \)")
erraxs[2].set_title("Finite Differences ({} Points)".format(spec_trunc**3))


if showfigs():
    plt.show()
else:
    for i, fig in enumerate([slicefig, pcfig, errfig]):
        for ext in ['pdf', 'svg']:
            fig.savefig(imgpath("psm_solvis-{}.{}".format(i+1, ext)))


