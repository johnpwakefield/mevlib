

import numpy as np


# solution pointwise

def u(a, R, r):
    return (R / r) * np.sinh(a * r) / np.sinh(a * R)


# solution integrated over domain

def s(a, R):
    return 3.0 / (R * a) * (np.coth(R * a) - (R * a)**(-1))


