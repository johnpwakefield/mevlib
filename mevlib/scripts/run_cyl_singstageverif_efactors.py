#!/usr/bin/env python3


from mevlib.shapes import Cylinder


nexact, kexact = 128, 128


D = 5.4423
cases = [
    (240.0, 240.0, 1.0),
    (140.0, 840.0, 1.0),
    (180.0, 360.0, 4.0),
    (180.0, 360.0, 0.1),
    (180.0, 360.0, 1.0)
]
cyls = map(lambda tpl: Cylinder(*tpl[:2]), cases)
rates = [1.5117e-3, 1.5117e-3, 6.0470e-3, 1.5117e-4, 1.5117e-3]
phi2s = [rate / D for rate in rates]


etas = [
    cyl.intgtd(phi2, {'ntrunc': nexact, 'ktrunc': kexact})
    for cyl, phi2 in zip(cyls, phi2s)
]


print("Case\tR\tH\tphi2\teta")
for i, ((R, H, phi2), eta) in enumerate(zip(cases, etas)):
    print("\t".join(list(map(str, [i+1, R, H, phi2, eta]))))


