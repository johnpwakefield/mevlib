#!/usr/bin/env python3


import numpy as np

from matplotlib import pyplot as plt

from nrel_cylinders import u1, u2


plt.rc('text', usetex=True)
plt.rc('axes', labelsize=12)


a = 00.8
R = 04.0
H = 15.0
rs, zs = np.linspace(0.0, 1.05 * R, 120), np.linspace(-0.05 * H, 1.05 * H, 320)
rmesh, zmesh = np.meshgrid(rs, zs)


truncs = [2, 6, 12, 24]
slicetrunc = 80


fig1, axs1 = plt.subplots(len(truncs), len(truncs))
fig2, axs2 = plt.subplots(len(truncs), len(truncs))

for fig, axs, xlims, ylims, zlim, in [
    (fig1, axs1, None, None, 2.0),
    (fig2, axs2, (0.0, R), (0.0 * H, 1.0 * H), None)
     ]:
    for i, trunc1 in enumerate(truncs):
        for j, trunc2 in enumerate(truncs):
            vals = (
                u1(a, R, H, trunc1, rmesh, zmesh) +
                u2(a, R, H, trunc2, rmesh, zmesh)
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


fig3, axs3 = plt.subplots(2, len(truncs))
for i, trunc in enumerate(truncs):
    for j, f in enumerate([u1, u2]):
        cf = axs3[j, i].contourf(rmesh, zmesh, f(a, R, H, trunc, rmesh, zmesh))
        axs3[j,i].set_xlim((0.0, R))
        axs3[j,i].set_ylim((0.0, H))
        fig3.colorbar(cf, ax=axs3[j,i])
        axs3[j,i].set_xlabel("r")
        axs3[j,i].set_ylabel("z")
    axs3[0, i].set_title("{}".format(trunc))


fig4, axs4 = plt.subplots(2, 2)
for i, f in enumerate([u1, u2]):
    nd = "First" if i == 0 else "Second"
    axs4[i,0].plot(rs, f(a, R, H, slicetrunc, rs, H/3))
    axs4[i,0].set_xlabel(r"\( r \)")
    axs4[i,0].set_ylim((-0.5, 1.5))
    axs4[i,0].set_title(r"{} sum, slice at \( z = H / 3 \)".format(nd))
    axs4[i,1].plot(zs, f(a, R, H, slicetrunc, R/3, zs))
    axs4[i,1].set_xlabel(r"\( z \)")
    axs4[i,1].set_ylim((-0.5, 1.5))
    axs4[i,1].set_title(r"{} sum, slice at \( r = R / 3 \)".format(nd))
    axs4[i,0].grid()
    axs4[i,1].grid()


fig5, ax5 = plt.subplots(1, 1, figsize=(3.0,3.0))
vals = (
    u1(a, R, H, slicetrunc, rmesh, zmesh) +
    u2(a, R, H, slicetrunc, rmesh, zmesh)
)
vals = np.maximum(np.minimum(vals, 1.5), -0.5)
cf = ax5.contourf(rmesh, zmesh, vals)
ax5.set_xlim(xlims)
ax5.set_ylim(ylims)
fig5.colorbar(cf, ax=ax5)
ax5.set_xlabel(r"\( r \)")
ax5.set_ylabel(r"\( z \)")


if False:
    for i, fig in enumerate([fig1, fig2, fig3, fig4, fig5]):
        fig.tight_layout()
        for ext in ['svg', 'pdf']:
            fig.savefig("img/graphic-{}.{}".format("fig{}".format(i+1), ext))
else:
    plt.show()


