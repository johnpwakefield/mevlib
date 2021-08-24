#!/usr/bin/env python3


import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

from mevlib.scalar import rct_ptwise_diff


dims = (4.0, 2.0)
M, N = 321, 177
a = 2.7


fig, ax = plt.subplots(1, 1, figsize=dims)

xs = np.linspace(0.0, dims[0], M+2, endpoint=True)
ys = np.linspace(0.0, dims[1], N+2, endpoint=True)
xs, ys = xs[1:-1], ys[1:-1]
xmesh, ymesh = np.meshgrid(xs, ys, indexing='xy')
zmesh1 = rct_ptwise_diff(a**2, dims[0], dims[1], M, N)
zmesh2 = (
    np.exp(-a * xmesh) + np.exp(-a * (dims[0] - xmesh)) +
    np.exp(-a * ymesh) + np.exp(-a * (dims[1] - ymesh))
)
zmesh1[np.logical_and(ymesh > dims[1] / 2, xmesh < dims[0] / 2)] = 0.0
zmesh1[np.logical_and(ymesh < dims[1] / 2, xmesh > dims[0] / 2)] = 0.0
zmesh2[np.logical_and(ymesh <= dims[1] / 2, xmesh <= dims[0] / 2)] = 0.0
zmesh2[np.logical_and(ymesh >= dims[1] / 2, xmesh >= dims[0] / 2)] = 0.0
zmesh = zmesh1 + zmesh2
cl = ax.contourf(
    xmesh, ymesh, zmesh,
    cmap=cm.Blues, levels=42,
    antialiased=True
)
for c in cl.collections:
    c.set_edgecolor("face")
ax.axis('off')

plt.savefig('img/logobg.pdf')
plt.savefig('img/logobg.png')




