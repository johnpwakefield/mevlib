#!/usr/bin/env python3


import numpy as np
from mevlib.maketable import compute_matrices
from getvals import getvals


fname = "../examples/example_file.sbl"

conc = np.array([
    [0.1, 0.2, 0.3, 0.4, 0.5, 0.6],
    [0.0, 0.1, 0.3, 0.4, 0.2, 0.0]
]).T


pointmats, pointtemps = compute_matrices(fname, True, False)


w = 0.666
intertemps = [
    (1 - w) * pointtemps[i] + w * pointtemps[i+1]
    for i in range(len(pointtemps) - 1)
]
intermats = [
    (1 - w) * pointmats[i] + w * pointmats[i+1]
    for i in range(len(pointtemps) - 1)
]

pointtbl = getvals(
    np.tile(conc, len(pointtemps)),
    np.repeat(np.array(pointtemps), conc.shape[1])
)
intertbl = getvals(
    np.tile(conc, len(intertemps)),
    np.repeat(np.array(intertemps), conc.shape[1])
)
pointlib = np.array([
    np.dot(mat, conc[:,j]) for mat in pointmats for j in range(conc.shape[1])
]).T
interlib = np.array([
    np.dot(mat, conc[:,j]) for mat in intermats for j in range(conc.shape[1])
]).T

pointerr = max(np.linalg.norm(pointtbl - pointlib, axis=0))
intererr = max(np.linalg.norm(intertbl - interlib, axis=0))

print(pointerr, intererr)



