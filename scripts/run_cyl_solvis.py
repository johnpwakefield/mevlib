#!/usr/bin/env python3


import numpy as np

from matplotlib import pyplot as plt

from mevlib.shapes import Cylinder


plt.rc('font', size=10)
plt.rc('text', usetex=True)
plt.rc('axes', labelsize=10)
plt.rc('legend', fontsize=10)

axsize = (2.75, 2.5)

def figsize(i, j):
    return (axsize[0] * j, axsize[1] * i)


a2 = 0.8**2
R, H = 04.0, 15.0
cyl = Cylinder(R, H)
rs, zs = np.linspace(0.0, 1.05 * R, 120), np.linspace(-0.05 * H, 1.05 * H, 320)
rmesh, zmesh = np.meshgrid(rs, zs)


truncs = [2, 6, 12, 24]
slicetrunc = 80


fig1, axs1 = plt.subplots(
    len(truncs), len(truncs), figsize=figsize(len(truncs), len(truncs))
)
fig2, axs2 = plt.subplots(
    len(truncs), len(truncs), figsize=figsize(len(truncs), len(truncs))
)

for fig, axs, xlims, ylims, zlim, in [
    (fig1, axs1, None, None, 2.0),
    (fig2, axs2, (0.0, R), (0.0 * H, 1.0 * H), None)
     ]:
    for i, trunc1 in enumerate(truncs):
        for j, trunc2 in enumerate(truncs):
            vals = (
                np.vectorize(cyl.ptwise_radial)(a2, trunc1, rmesh, zmesh)
                +
                np.vectorize(cyl.ptwise_axial)(a2, trunc2, rmesh, zmesh)
            )
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


fig3, axs3 = plt.subplots(2, len(truncs), figsize=figsize(2, len(truncs)))
for i, trunc in enumerate(truncs):
    for j, f in enumerate(
        map(np.vectorize, [cyl.ptwise_radial, cyl.ptwise_axial])
    ):
        cf = axs3[j, i].contourf(rmesh, zmesh, f(a2, trunc, rmesh, zmesh))
        axs3[j, i].set_xlim((0.0, R))
        axs3[j, i].set_ylim((0.0, H))
        fig3.colorbar(cf, ax=axs3[j, i])
        axs3[j, i].set_xlabel("r")
        axs3[j, i].set_ylabel("z")
    axs3[0, i].set_title("{}".format(trunc))


fig4, axs4 = plt.subplots(2, 2, figsize=figsize(2, 2))
fig4_notitle, axs4_notitle = plt.subplots(2, 2, figsize=figsize(2, 2))
for n, (fig, axs) in enumerate([(fig4, axs4), (fig4_notitle, axs4_notitle)]):
    for i, f in enumerate(map(
        np.vectorize, [cyl.ptwise_radial, cyl.ptwise_axial]
    )):
        nd = "First" if i == 0 else "Second"
        axs[i, 0]
        axs[i, 0].plot(rs, f(a2, slicetrunc, rs, H/3))
        axs[i, 0].set_ylim((-0.5, 1.5))
        axs[i, 0].set_xlabel(r"\( r \)")
        axs[i, 0].set_ylabel(r"\( \hat{Y} \)")
        if n == 0:
            axs[i, 0].set_title(
                r"{} sum, slice at \( z = H / 3 \)".format(nd)
            )
        axs[i, 1].plot(zs, f(a2, slicetrunc, R/3, zs))
        axs[i, 1].set_ylim((-0.5, 1.5))
        axs[i, 1].set_xlabel(r"\( z \)")
        axs[i, 1].set_ylabel(r"\( \hat{Y} \)")
        if n == 0:
            axs[i, 1].set_title(
                r"{} sum, slice at \( r = R / 3 \)".format(nd)
            )
        axs[i, 0].grid()
        axs[i, 1].grid()


fig5, ax5 = plt.subplots(1, 1, figsize=figsize(1, 1))
vals = (
    np.vectorize(cyl.ptwise_radial)(a2, slicetrunc, rmesh, zmesh) +
    np.vectorize(cyl.ptwise_axial)(a2, slicetrunc, rmesh, zmesh)
)
vals = np.maximum(np.minimum(vals, 1.5), -0.5)
cf = ax5.contourf(rmesh, zmesh, vals)
ax5.set_xlim(xlims)
ax5.set_ylim(ylims)
cbar = fig5.colorbar(cf, ax=ax5)
cbar.set_label(r"\( \hat{Y} \)")
ax5.set_xlabel(r"\( r \)")
ax5.set_ylabel(r"\( z \)")


for fig, figname in [
    (fig1, "fig1"),
    (fig2, "fig2"),
    (fig3, "fig3"),
    (fig4, "fig4"), (fig4_notitle, "fig4_notitle"),
    (fig5, "fig5")
]:
    fig.tight_layout()
    for ext in ['svg', 'pdf']:
        fig.savefig("img/cyl_solvis_{}.{}".format(figname, ext))


