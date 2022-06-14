#!/usr/bin/env python3


from functools import reduce

import pkgutil
from io import StringIO

import numpy as np
import scipy.linalg as la

from matplotlib import pyplot as plt
import matplotlib.ticker as mtick

from mevlib.mechanisms import Mechanism
from mevlib.shapes import Cylinder
from mevlib.parsing.auto import parse_attempt
from mevlib.diagonalization import PointDiagonalization, computeintegral


plt.rc('font', size=10)
plt.rc('text', usetex=True)
plt.rc('axes', labelsize=10)
plt.rc('legend', fontsize=8)
plt.rc('xtick', labelsize=8)
plt.rc('ytick', labelsize=8)


rad, hgt = 90.0, 360.0
Ts = np.linspace(500.0, 900.0, 300)
Cbdry = np.array([36.0, 32.0, 28.0, 24.0, 8.0])
mechfile = 'data/fcc_lattanzietal_2020.sbl'
components = ['S', 'D', 'Gas', 'LPG', 'DR']


shpe = Cylinder(rad, hgt)


chi = 2
beta = 0.0
def I(a2, R):           #noqa E743
    a = np.sqrt(a2)
    th = np.tanh(a * R)
    return (R / a - a**(-2) * th) / (beta * (a * R - th) + R * th)
def surfappx(shpe, mech, T, precision):
    lams, R = la.eig(mech.getB(T)[:mech.Ng, :])
    if np.max(np.abs(lams.imag)) > 1e-8:
        print("Matrix has imaginary eigenvalues (irreversible).")
    lams, R = np.real_if_close(lams), np.real_if_close(R)
    SAoV = 2 * (rad + hgt) / (rad * hgt)
    mult = np.array([
        SAoV * I(a2, np.sqrt(0.5 * rad**2 + 0.5 * rad * hgt)) for a2 in lams
    ])
    mult = np.minimum(np.where(np.isnan(mult), 1.0, mult), 1.0)
    A = np.dot(R, np.dot(np.diag(mult), la.inv(R)))
    return A


figfull, axsfull = plt.subplots(5, 2, figsize=(5.0, 9.0))
figtrnc, axstrnc = plt.subplots(3, 2, figsize=(5.0, 5.0))


# parse test mechanism
mechconf = pkgutil.get_data('mevlib', mechfile)
precision, shape, temperatures, _, species, reactions = parse_attempt(
    StringIO(mechconf.decode('utf-8')), '.sbl', True, True
)
if precision is None:
    precision = {'ntrunc': 128, 'ktrunc': 64}
mech = Mechanism(species, reactions)


# compute our solution, plot effectiveness factors
libmevs = np.empty((len(Ts), len(components)))
for i, T in enumerate(Ts):
    diag = PointDiagonalization(mech, T)
    libmevs[i, :] = reduce(np.dot, [
        diag.get_evects(),
        np.diag(computeintegral(shpe, diag.get_evals(), precision)),
        diag.get_evectinv(),
        Cbdry[:mech.Ng]
    ]) / Cbdry[:mech.Ng]
apxmevs = np.vstack([
    np.dot(surfappx(shpe, mech, T, precision), Cbdry) / Cbdry[0]
    for T in Ts
])

for i in range(len(components)):
    for axs in [axsfull] + ([axstrnc] if i < 3 else []):
        axs[i, 0].plot(Ts, libmevs[:, i], label="Actual")
        axs[i, 0].plot(Ts, apxmevs[:, i], label="Approximate")
        axs[i, 1].plot(
            Ts, np.abs(libmevs[:, i] - apxmevs[:, i]) / libmevs[:, i], "k-"
        )
        axs[i, 1].yaxis.set_major_formatter(mtick.PercentFormatter(1.0))
        axs[i, 0].set_ylabel(r"\( \eta_{} \)".format(i+1))
        axs[i, 1].set_ylabel("Relative Error")
        for j in range(2):
            axs[i, j].set_xlabel("Temperature")
            axs[i, j].grid()
        axs[i, 0].set_ylim((0.3, 1.2))
        axs[i, 1].set_ylim((0.0, 0.3))


axsfull[0, 0].legend()
axstrnc[0, 0].legend()


for fig, n in [(figfull, ""), (figtrnc, "_trunc")]:

    fig.tight_layout()

    for ext in ['svg', 'pdf']:
        fig.savefig("img/SA_appx_multistage_cyl{}.{}".format(n, ext))


