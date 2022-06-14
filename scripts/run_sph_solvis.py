#!/usr/bin/env python3


import numpy as np
import matplotlib.pyplot as plt
from mevlib.scalar import sph_ptwise


plt.rc('font', size=10)
plt.rc('text', usetex=True)
plt.rc('axes', labelsize=10)
plt.rc('legend', fontsize=10)

axsize = (2.75, 2.5)

def figsize(i, j):
    return (axsize[0] * j, axsize[1] * i)


R = 120.0
a2 = 0.04**2
rs = np.linspace(0.0, R, 300)


fig1, ax1 = plt.subplots(1, 1, figsize=figsize(1, 1))
ax1.plot(rs, sph_ptwise(a2, R, rs), 'b-')
ax1.grid()
ax1.set_xlabel(r'\( r \)')
ax1.set_ylabel(r'\( \hat{Y} \)')


fig2, ax2 = plt.subplots(1, 1, figsize=figsize(1, 1))
xs = np.linspace(-R, R, 543)
z = 0.05 * R
xmesh, ymesh = np.meshgrid(xs, xs)
umesh = sph_ptwise(a2, R, np.sqrt(xmesh**2 + ymesh**2 + z**2))
umesh[R**2 < xmesh**2 + ymesh**2 + z**2] = np.nan
clp = ax2.contourf(xmesh, ymesh, umesh)
ax2.set_xlabel(r'\( x \)')
ax2.set_ylabel(r'\( y \)')
cbar = fig2.colorbar(clp, fraction=0.046, pad=0.04)
cbar.set_label(r"\( \hat{Y} \)")
ax2.set_aspect('equal')


for i, fig in enumerate([fig1, fig2]):
    fig.tight_layout()
    for ext in ['pdf', 'svg']:
        fig.savefig("img/sph_solvis-{}.{}".format(i+1, ext))


