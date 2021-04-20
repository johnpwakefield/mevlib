

import numpy as np


# solution pointwise

def u(a, R, r):
    return np.cosh(a * r) / np.cosh(a * R)


# solution integrated over domain

def s(a, R):
    return (
        (3.0 / a + 18.0 / a**3 / R**2) * np.tanh(a * R)
        - (9.0 / R / a**2 + 18.0 / R**3 / a**4)
    )


