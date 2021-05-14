#!/usr/bin/env python3


import subprocess


cmd = "./compile.sh"
process = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE)
out, err = process.communicate()


import numpy as np
from getvals import getvals


concs = np.array([
    [0.1, 0.2, 0.3, 0.4, 0.5, 0.6],
    [0.0, 0.1, 0.3, 0.4, 0.2, 0.0]
]).T
temps = np.array([678.9, 777.7])


print(getvals(concs, temps))


