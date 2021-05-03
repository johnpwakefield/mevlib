"""

Contains functions that yield pointwise and integrated solutions to the scalar
problem for a variety of shapes.

"""

import numpy as np
from scipy.special import j0, j1, i0e, i1e, jn_zeros
from scipy.linalg import solve
from scipy.sparse.linalg import gmres, spsolve
from scipy.sparse import dok_matrix as spmat
from scipy.interpolate import RegularGridInterpolator as RGI


# general functions

def sum_standard(terms):
    return np.sum(terms[::-1])

def sum_aitken(terms):
    anp1, anp2, anp3 = terms[0], terms[1], terms[2]
    cml = anp1 - anp2**2 / (anp3 - anp2)
    for n in range(len(terms)-3):
        anp1, anp2, anp3 = anp2, anp3, terms[n+3]
        cml += anp1 - anp2**2 / (anp3 - anp2) + anp1**2 / (anp2 - anp1)
    return cml


# cylinder functions

def cyl_gn(a2, H, n):
    return np.sqrt(a2 + (np.pi * (2 * n + 1) / H)**2)

def cyl_ln(a2, R, alphak):
    return np.sqrt(a2 + (alphak / R)**2)

def cyl_ptwise_radial_terms(a2, R, H, trunc, r, z):
    terms = np.empty((trunc,))
    for n in range(trunc-1, -1, -1):
        gn = cyl_gn(a2, H, n)
        terms[n] = (
            4.0 / (np.pi * (2 * n + 1))
            * np.sin(np.pi * (2 * n + 1) * z / H)
            * np.exp(gn * (r - R)) * i0e(gn * r) / i0e(gn * R)
        )
    return terms

def cyl_ptwise_axial_terms(a2, R, H, trunc, r, z, zero_cache=None):
    if zero_cache is not None and len(zero_cache) >= trunc:
        alpha = zero_cache
    else:
        alpha = jn_zeros(0, trunc)
    terms = np.empty((trunc,))
    for k in range(trunc, 0, -1):
        lk = cyl_ln(a2, R, alpha[k-1])
        terms[k-1] = (
            2.0 * np.cosh(lk * (z - H / 2)) / np.cosh(lk * H / 2)
            * j0(alpha[k-1] / R * r) / (j1(alpha[k-1]) * alpha[k-1])
        )
    return terms

def cyl_ptwise_radial(a2, R, H, ntrunc, r, z):
    return sum_standard(cyl_ptwise_radial_terms(a2, R, H, ntrunc, r, z))

def cyl_ptwise_axial(a2, R, H, ktrunc, r, z, zero_cache=None):
    return sum_standard(cyl_ptwise_axial_terms(
        a2, R, H, ktrunc, r, z, zero_cache=zero_cache
    ))

def cyl_ptwise(a2, R, H, ntrunc, ktrunc, r, z, zero_cache=None):
    return sum(map(sum_standard, [
        cyl_ptwise_radial_terms(a2, R, H, ntrunc, r, z),
        cyl_ptwise_axial_terms(a2, R, H, ktrunc, r, z, zero_cache)
    ]))

cyl_ptwise_radial_vect = np.vectorize(cyl_ptwise_radial)
cyl_ptwise_axial_vect = np.vectorize(cyl_ptwise_axial)
cyl_ptwise_vect = np.vectorize(cyl_ptwise)

def cyl_intgtd_radial_terms(a2, R, H, trunc, bratcutoff=None):
    ns = np.arange(trunc)
    gns = cyl_gn(a2, H, ns)
    # the ratio I1/I0 tends to 1 (inf/inf) very quickly
    #TODO instead of computing then checking, just compute what is needed
    brat = i1e(gns * R) / i0e(gns * R)
    if bratcutoff is not None:
        brat = np.where(gns * R < bratcutoff, brat, 1.0)
    return 16.0 / (np.pi**2 * R * gns * (2 * ns + 1)**2) * brat

def cyl_intgtd_axial_terms(a2, R, H, trunc, zero_cache=None):
    if zero_cache is not None and len(zero_cache) >= trunc:
        alpha = zero_cache
    else:
        alpha = jn_zeros(0, trunc)
    lnks = cyl_ln(a2, R, alpha)
    return 8.0 * np.tanh(lnks * H / 2) / (H * lnks * alpha**2)

#TODO functions like this might be confusing since they aren't used directly
# in shapes.py
def cyl_intgtd(a2, R, H, ntrunc, ktrunc, bratcutoff=None, zero_cache=None):
    return sum(map(sum_standard, [
        cyl_intgtd_radial_terms(a2, R, H, ntrunc, bratcutoff=bratcutoff),
        cyl_intgtd_axial_terms(a2, R, H, ktrunc, zero_cache=zero_cache)
    ]))


# prism functions

def psm_beta(a2, L1, L2, L3, ms, ns):
    return np.sqrt(
        (np.pi * (2 * ms + 1) * (L1 / L2))**2 +
        (np.pi * (2 * ns + 1) * (L1 / L3))**2 +
        L1**2 * a2
    )

def psm_ptwise_xdir(a2, Lx, Ly, Lz, trunc, x, y, z):
    ms, ns = np.meshgrid(np.arange(trunc)[::-1], np.arange(trunc)[::-1])
    betas = psm_beta(a2, Lx, Ly, Lz, ms, ns)
    return 16.0 / np.pi**2 * np.sum(
        np.sin(np.pi * (2 * ms + 1) * y / Ly) *
        np.sin(np.pi * (2 * ns + 1) * z / Lz) *
        (
            #TODO add truncation for large values
            np.cosh(betas * (x / Lx - 0.5)) / np.cosh(betas * 0.5)
        ) / ((2 * ms + 1) * (2 * ns + 1))
    )

def psm_ptwise_series(a2, Lx, Ly, Lz, xtrunc, ytrunc, ztrunc, x, y, z):
    return sum([
        psm_ptwise_xdir(a2, L1, L2, L3, trunc, c1, c2, c3)
        for c1, c2, c3, L1, L2, L3, trunc in [
            (x, y, z, Lx, Ly, Lz, xtrunc),
            (y, z, x, Ly, Lz, Lx, ytrunc),
            (z, x, y, Lz, Lx, Ly, ztrunc)
        ]
    ])

def psm_intgtd_xdir(a2, Lx, Ly, Lz, trunc):
    ms, ns = np.meshgrid(np.arange(trunc)[::-1], np.arange(trunc)[::-1])
    betas = psm_beta(a2, Lx, Ly, Lx, ms, ns)
    cml = 0.0
    for dg in range(2 * trunc - 1):
        for k in range(trunc - abs(trunc - 1 - dg)):
            if dg < trunc:
                i, j = dg - k, k
            else:
                i, j = trunc - 1 - k, dg - trunc + 1 + k
            cml += (
                np.tanh(betas[i,j] * 0.5) /
                (betas[i,j] * (2 * ms[i,j] + 1)**2 * (2 * ns[i,j] + 1)**2)
            )
    return 32.0 / np.pi**4 * cml
#   return 32.0 / np.pi**4 * np.sum(
#       np.tanh(betas * 0.5) / (betas * (2 * ms + 1)**2 * (2 * ns + 1)**2)
#   )

def psm_intgtd_series(a2, Lx, Ly, Lz, xtrunc, ytrunc, ztrunc):
    return sum([
        psm_intgtd_xdir(a2, L1, L2, L3, trunc)
        for L1, L2, L3, trunc in [
            (Lx, Ly, Lz, xtrunc), (Ly, Lz, Lx, ytrunc), (Lz, Lx, Ly, ztrunc)
        ]
    ])

def Tk(k, x):
    return np.cos(k * np.arccos(x))
def Tkxx(k, x):
    return k * (
        k * Tk(k, x) -
        x * np.sin((k + 1) * np.arccos(x)) / np.sqrt(1.0 - x**2)
    ) / (x**2 - 1.0)
def phi(n, x):
    return Tk(2 * (n + 1), x) - Tk(2 * n, x)
def phixx(n, x):
    return Tkxx(2 * (n + 1), x) - Tkxx(2 * n, x)

def psm_spec_intmesh(Ns):
    return np.meshgrid(*map(np.arange, Ns))

def psm_spec_coeffs(a2, Lx, Ly, Lz, truncs):
    if np.isscalar(truncs):
        truncs = truncs, truncs, truncs
    else:
        assert(len(truncs) == 3)
    Nx, Ny, Nz = truncs
    xs, ys, zs = map(
        lambda N: -np.cos(0.5 * np.pi * (1.0 + np.arange(N)) / N), [Nx, Ny, Nz]
    )
    xmesh, ymesh, zmesh = np.meshgrid(xs, ys, zs)
    xflat, yflat, zflat = map(lambda arr: arr.flatten(), [xmesh, ymesh, zmesh])
    del xmesh, ymesh, zmesh
    mmesh, nmesh, pmesh = psm_spec_intmesh((Nx, Ny, Nz))
    mflat, nflat, pflat = map(lambda arr: arr.flatten(), [mmesh, nmesh, pmesh])
    del mmesh, nmesh, pmesh
    mgrid, xgrid = np.meshgrid(mflat, xflat)
    del mflat, xflat
    ngrid, ygrid = np.meshgrid(nflat, yflat)
    del nflat, yflat
    pgrid, zgrid = np.meshgrid(pflat, zflat)
    del pflat, zflat
    xc, yc, zc = map(lambda L: 4.0 / L**2, [Lx, Ly, Lz])
    H = (
        a2 * phi(mgrid, xgrid) * phi(ngrid, ygrid) * phi(pgrid, zgrid)
        - xc * phixx(mgrid, xgrid) * phi(ngrid, ygrid) * phi(pgrid, zgrid)
        - yc * phi(mgrid, xgrid) * phixx(ngrid, ygrid) * phi(pgrid, zgrid)
        - zc * phi(mgrid, xgrid) * phi(ngrid, ygrid) * phixx(pgrid, zgrid)
    )
    sol = solve(H, -a2 * np.ones((Nx * Ny * Nz, 1)))
    sol, info = gmres(H, -a2 * np.ones((Nx * Ny * Nz, 1)), atol=1e-11, x0=sol)
    return (Lx, Ly, Lz, sol.reshape((Nx, Ny, Nz)))

def psm_spec_eval(coeffs, x, y, z):
    xh, yh, zh = map(
        lambda cL: (cL[0] - cL[1] / 2) / (cL[1] / 2),
        zip((x, y, z), coeffs[:3])
    )
    mmesh, nmesh, pmesh = psm_spec_intmesh(coeffs[-1].shape)
    return 1.0 + np.sum(
        coeffs[-1] * phi(mmesh, xh) * phi(nmesh, yh) * phi(pmesh, zh)
    )

def psm_intgtd_spec_eval(coeffs):
    Lx, Ly, Lz, coeffs = coeffs
    mmesh, nmesh, pmesh = map(
        lambda arr: arr[::-1, ::-1, ::-1], psm_spec_intmesh(coeffs.shape)
    )
    return 1.0 - 64.0 * np.sum(
        coeffs[::-1, ::-1, ::-1] * (
            (1.0 - 4 * mmesh**2) * (3.0 + 2.0 * mmesh) *
            (1.0 - 4 * nmesh**2) * (3.0 + 2.0 * nmesh) *
            (1.0 - 4 * pmesh**2) * (3.0 + 2.0 * pmesh)
        )**(-1)
    )

def psm_intgtd_spec(a2, Lx, Ly, Lz, truncs, coeffs=None):
    if coeffs is None:
        coeffs = psm_spec_coeffs(a2, Lx, Ly, Lz, truncs)
    return psm_intgtd_spec_eval(coeffs)

def psm_diff_axes(Lx, Ly, Lz, Nx, Ny, Nz):
    xs, ys, zs = map(
        lambda LN: LN[0] * (0.5 + np.arange(LN[1])) / LN[1],
        [(Lx, Nx), (Ly, Ny), (Lz, Nz)]
    )
    return xs, ys, zs

def psm_diff_coeffs(a2, Lx, Ly, Lz, truncs):
    if np.isscalar(truncs):
        truncs = truncs, truncs, truncs
    else:
        assert(len(truncs) == 3)
    Nx, Ny, Nz = truncs
    xs, ys, zs = psm_diff_axes(Lx, Ly, Lz, Nx, Ny, Nz)
    hx, hy, hz = xs[1] - xs[0], ys[1] - ys[0], zs[1] - zs[0]
    A = spmat((Nx * Ny * Nz, Nx * Ny * Nz))
    rhs = np.zeros((Nx, Ny, Nz))
    for m in range(Nx):
        for n in range(Ny):
            for p in range(Nz):
                row = np.zeros((Nx, Ny, Nz))
                row[m,n,p] = -2.0 * (hx**-2 + hy**-2 + hz**-2 + a2 / 2)
                if m + 1 != Nx:
                    row[m+1,n,p] = hx**-2
                else:
                    rhs[m,n,p] = -hx**-2
                if m != 0:
                    row[m-1,n,p] = hx**-2
                else:
                    rhs[m,n,p] = -hx**-2
                if n + 1 != Ny:
                    row[m,n+1,p] = hy**-2
                else:
                    rhs[m,n,p] = -hy**-2
                if n != 0:
                    row[m,n-1,p] = hy**-2
                else:
                    rhs[m,n,p] = -hy**-2
                if p + 1 != Nz:
                    row[m,n,p+1] = hz**-2
                else:
                    rhs[m,n,p] = -hz**-2
                if p != 0:
                    row[m,n,p-1] = hz**-2
                else:
                    rhs[m,n,p] = -hz**-2
                A[np.ravel_multi_index((m, n, p), (Nx, Ny, Nz)), :] = (
                    row.flatten()
                )
    rhs = rhs.flatten()
    return (Lx, Ly, Lz, spsolve(A, rhs).reshape((Nx, Ny, Nz)))

def psm_diff_eval(coeffs, x, y, z):
    Lx, Ly, Lz = coeffs[:3]
    Nx, Ny, Nz = coeffs[-1].shape
    xs, ys, zs = psm_diff_axes(Lx, Ly, Lz, Nx, Ny, Nz)
    cfs = coeffs[-1].flatten().reshape((Nx, Ny, Nz), order='F')
    interp = RGI((xs, ys, zs), cfs, bounds_error=False, fill_value=-1)
    return interp([x, y, z])

def psm_intgtd_diff_eval(coeffs):
    Lx, Ly, Lz = coeffs[:3]
    Nx, Ny, Nz = coeffs[-1].shape
    return np.sum(coeffs[-1]) / (Nx * Ny * Nz)

def psm_intgtd_diff(a2, Lx, Ly, Lz, truncs, coeffs=None):
    if coeffs is None:
        coeffs = psm_diff_coeffs(a2, Lx, Ly, Lz, truncs)
    return psm_intgtd_diff_eval(coeffs)

psm_ptwise = psm_ptwise_series
psm_intgtd = psm_intgtd_series


# sphere functions

def sph_ptwise(a2, R, r):
    return (R / r) * np.sinh(np.sqrt(a2) * r) / np.sinh(np.sqrt(a2) * R)

def sph_intgtd(a2, R):
    a = np.sqrt(a2)
    return 3.0 / (R * a) * (np.coth(R * a) - (R * a)**(-1))


