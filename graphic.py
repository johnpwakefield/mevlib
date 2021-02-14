#!/usr/bin/env python3


from random import random

import numpy as np

from matplotlib import pyplot as plt
import plotly.graph_objects as go

from nrel_cylinders import u1, u2, u


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
rs, zs = np.linspace(0.0, 1.05 * R, 120), np.linspace(-0.05 * H, 1.05 * H, 320)
rmesh, zmesh = np.meshgrid(rs, zs)


truncs = [2, 6, 12, 24]
slicetrunc = 80
numcheckpts = 2000
chkntrunc, chkktrunc = 24, 24
chkdr, chkdz = 0.01 * R, 0.01 * H


fig1, axs1 = plt.subplots(len(truncs), len(truncs))
fig2, axs2 = plt.subplots(len(truncs), len(truncs))

for fig, axs, xlims, ylims, zlim, in [
    (fig1, axs1, None, None, 1.0),
    (fig2, axs2, (0.0, 0.95 * R), (0.05 * H, 0.95 * H), 1.0)
     ]:
    for i, trunc1 in enumerate(truncs):
        for j, trunc2 in enumerate(truncs):
            vals = (
                u1(a, R, H, trunc1, rmesh, zmesh) +
                u2(a, R, H, trunc2, rmesh, zmesh)
            )
            if zlim is not None:
                vals = np.minimum(np.maximum(vals, 0.5 - zlim), 0.5 + zlim)
            cf = axs[i, j].contourf(rmesh, zmesh, vals)
            axs[i, j].set_xlim(xlims)
            axs[i, j].set_ylim(ylims)
            fig.colorbar(cf, ax=axs[i, j])
            axs[i, j].set_title("{}, {}".format(trunc1, trunc2))
            axs[i, j].set_xlabel("r")
            axs[i, j].set_ylabel("z")


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
vals = np.maximum(np.minimum(
    u1(a, R, H, slicetrunc, rmesh, zmesh) +
    u2(a, R, H, slicetrunc, rmesh, zmesh),
    1.5), -0.5)
cf = ax5.contourf(rmesh, zmesh, vals)
ax5.set_xlim(xlims)
ax5.set_ylim(ylims)
fig5.colorbar(cf, ax=ax5)
ax5.set_xlabel(r"\( r \)")
ax5.set_ylabel(r"\( z \)")


fig6, ax6 = plt.subplots(1, 1, figsize=(3.0,3.0))
chkrs = [R * random() for i in range(numcheckpts)]
chkzs = [H * random() for i in range(numcheckpts)]
vals = [
    abs(checkpoint(a, R, H, chkntrunc, chkktrunc, chkdr, chkdz, r, z))
    for (r, z) in zip(chkrs, chkzs)
]
sc = ax6.scatter(chkrs, chkzs, c=vals, s=3.6)
fig.colorbar(sc, ax=ax6)


N = 40

X, Y, Z = map(lambda l: l.flatten(), np.mgrid[-R:R:N*1j,-R:0:N*1j,0:H:N*1j])
vals = np.array([
    abs(u(a, R, H, slicetrunc, slicetrunc, np.sqrt(x**2 + y**2), z))
    for (x, y, z) in zip(X, Y, Z)
])
fig7 = go.Figure(data=go.Volume(
    x=X.flatten(),
    y=Y.flatten(),
    z=Z.flatten(),
    value=vals.flatten(),
    isomin=-0.1,
    isomax=1.0,
    opacity=0.2,        # needs to be small to see through all surfaces
    surface_count=21,   # needs to be a large number for good volume rendering
))

yc = R / 8
thetac, zc = np.mgrid[0:2*np.pi:N*1j,0:H:N*1j]
xc = R * np.cos(thetac)
yc = np.minimum(R * np.sin(thetac), yc)
cc = np.empty_like(thetac)
thetatb, rtb = np.mgrid[0:2*np.pi:N*1j,0:R:N*1j]
xt, yt = rtb * np.cos(thetatb), np.minimum(rtb * np.sin(thetatb), yc)
zt, zb = np.zeros_like(thetatb), H * np.ones_like(thetatb)
ct, cb = np.empty_like(thetatb), np.empty_like(thetatb)
for i in range(N):
    for j in range(N):
        cc[i,j] = u(
            a, R, H, slicetrunc, slicetrunc,
            np.sqrt(xc[i,j]**2 + yc[i,j]**2), zc[i,j]
        )
        ct[i,j] = u(
            a, R, H, slicetrunc, slicetrunc,
            np.sqrt(xt[i,j]**2 + yt[i,j]**2), zt[i,j]
        )
fig8 = go.Figure(data=[
    go.Surface(x=xc, y=yc, z=zc, surfacecolor=cc),
    go.Surface(x=xt, y=yt, z=zt, surfacecolor=ct),
    go.Surface(x=xt, y=yt, z=zb, surfacecolor=ct)
])


if True:
    for i, fig in enumerate([fig1, fig2, fig3, fig4, fig5, fig6]):
        fig.tight_layout()
        for ext in ['svg', 'pdf']:
            fig.savefig("img/graphic-{}.{}".format("fig{}".format(i+1), ext))
    for i, fig in enumerate([fig7, fig8]):
        fig.write_html("img/graphic-fig{}.html".format(7+i))
else:
    for fig in [fig7, fig8]:
        fig.show()
    plt.show()


