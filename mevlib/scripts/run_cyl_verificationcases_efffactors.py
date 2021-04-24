#!/usr/bin/env python3


from matplotlib import pyplot as plt

from mevlib.scalar import cyl_intgtd
from scipy.special import jn_zeros


plt.rc('font', size=16)
plt.rc('text', usetex=True)
plt.rc('axes', labelsize=20)


nexact, kexact = 128, 128
zc = jn_zeros(0, kexact)


D = 5.4423
cases = [
    (240.0, 240.0, 60.0, 1.0),
    (140.0, 840.0, 60.0, 1.0),
    (180.0, 360.0, 60.0, 4.0),
    (180.0, 360.0, 60.0, 0.1),
    (180.0, 360.0, 60.0, 1.0)
]
rates = [1.5117e-3, 1.5117e-3, 6.0470e-3, 1.5117e-4, 1.5117e-3]
compphi2s = [rate * L**2 / D for (_, _, L, _), rate in zip(cases, rates)]
cases = [(*case, compphi2) for case, compphi2 in zip(cases, compphi2s)]


etas = [
    cyl_intgtd(phi2, R / L, H / L, nexact, kexact, zero_cache=zc)
    for R, H, L, phi2, cphi2 in cases
]


print("Case\tR\tH\tphi2\teta")
for i, ((R, H, L, phi2, cphi2), eta) in enumerate(zip(cases, etas)):
    print("\t".join(list(map(str, [i+1, R, H, phi2, eta]))))


