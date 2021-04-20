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
N = 240
yc = R / 8
zc = 0.4*H


fullscale = [
    [0.000, "green"],
    [0.001, "green"],
    [0.002, "blue"],
    [1.000, "orange"]
]
green = [
    [0.0, "green"],
    [1.0, "green"]
]

cyl1theta, cyl1z = np.mgrid[0:2*np.pi:N*1j,0:zc:N*1j]
cyl1u = np.zeros_like(cyl1theta)

flat1theta, flat1r = np.mgrid[0:2*np.pi:N*1j,0:R:N*1j]
flat1z = zc * np.ones_like(flat1theta)
flat1u = u(a, R, H, slicetrunc, slicetrunc, flat1r, zc)
flat1u[flat1r*np.sin(flat1theta) < -yc] = 0.0

cyl2theta, cyl2z = np.mgrid[0:2*np.pi:N*1j,zc:H:N*1j]
cyl2x = R * np.cos(cyl2theta)
cyl2y = np.minimum(R * np.sin(cyl2theta), yc)
cyl2u = u(a, R, H, slicetrunc, slicetrunc, np.sqrt(cyl2x**2 + cyl2y**2), cyl2z)
cyl2u[cyl2y != yc] = 1e-7
print(np.max(cyl2u), np.min(cyl2u))

flat2theta, flat2r = np.mgrid[0:2*np.pi:N*1j,0:R:N*1j]
flat2x = R * np.cos(flat2theta)
flat2y = np.minimum(R * np.sin(flat2theta), yc)
flat2z = H * np.ones_like(flat2theta)
flat2u = np.zeros_like(flat2theta)
flat2u[flat2r * np.sin(flat2theta) > yc] = np.nan


fig8 = go.Figure(data=[
    go.Surface(
        x=R*np.cos(cyl1theta), y=R*np.sin(cyl1theta), z=cyl1z,
        surfacecolor=cyl1u, colorscale=green
    ),
    go.Surface(
        x=cyl2x, y=cyl2y, z=cyl2z,
        surfacecolor=cyl2u, colorscale=fullscale
    ),
    go.Surface(
        x=flat1r*np.cos(flat1theta), y=flat1r*np.sin(flat1theta), z=flat1z,
        surfacecolor=flat1u, colorscale=fullscale
    ),
    go.Surface(
        x=flat2x, y=flat2y, z=flat2z,
        surfacecolor=flat2u, colorscale=green
    )
])
fig8.update_layout(
    scene={
        'xaxis': {'showticklabels': False, 'showgrid': False, 'title': ''},
        'yaxis': {'showticklabels': False, 'showgrid': False, 'title': ''},
        'zaxis': {'showticklabels': False, 'showgrid': False, 'title': ''}
    }
)


if True:
    for i, fig in enumerate([fig8]):
        fig.write_html("img/volrender-fig{}.html".format(8+i))
else:
    for fig in [fig8]:
        fig.show()


