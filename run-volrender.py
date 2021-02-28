#!/usr/bin/env python3


import numpy as np

import plotly.graph_objects as go

from nrel_cylinders import u


a = 00.8
R = 04.0
H = 15.0
rs, zs = np.linspace(0.0, 1.05 * R, 120), np.linspace(-0.05 * H, 1.05 * H, 320)
rmesh, zmesh = np.meshgrid(rs, zs)


slicetrunc = 80
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
    for i, fig in enumerate([fig7, fig8]):
        fig.write_html("img/volrender-fig{}.html".format(7+i))
else:
    for fig in [fig7, fig8]:
        fig.show()


