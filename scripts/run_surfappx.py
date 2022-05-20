#!/usr/bin/env python3


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick

from mevlib.scalar import sph_intgtd, cyl_intgtd


plt.rc('font', size=10)
plt.rc('text', usetex=True)
plt.rc('axes', labelsize=10)
plt.rc('legend', fontsize=8)
plt.rc('xtick', labelsize=8)
plt.rc('ytick', labelsize=8)


# the single fig has a huge legend; make it the same size
sidebysidefigsize = (5.5, 2.5)


V = np.pi * 140.0**2 * 840.0
beta = 0.0
lams = np.linspace(3.4e-5, 1.7e-3, 400)


def surf_rate(a2, depth_corr=None):
    if depth_corr is None:
        dc = 1.0
    else:
        dc = 1.0 - np.exp(- np.sqrt(a2) * depth_corr)
    return (beta * a2 + np.sqrt(a2))**(-1) * dc

chi = 2
def I(a2, R):           #noqa E743
    a = np.sqrt(a2)
    th = np.tanh(a * R)
    return (R / a - a**(-2) * th) / (beta * (a * R - th) + R * th)
def cylappx(cr, ch):
    SA, V = 2 * np.pi * cr * (cr + ch), np.pi * cr**2 * ch
    SAoV = 2 * (cr + ch) / (cr * ch)
    cyl_int = np.array([cyl_intgtd(a2, cr, ch, 128, 128) for a2 in lams])
    cyl_srt = np.minimum(np.array([surf_rate(a2) * SAoV for a2 in lams]), 1.0)
    cyl_dct = np.minimum(np.array([
        surf_rate(a2, depth_corr=cr) * SAoV for a2 in lams
    ]), 1.0)
    cyl_crv = np.minimum(np.array([
        SAoV * I(a2, np.sqrt(0.5 * cr**2 + 0.5 * cr * ch)) for a2 in lams
    ]), 1.0)
    cyl_crv1 = np.minimum(np.array([
        (
            4 * np.pi * cr * I(a2, np.sqrt(cr)) + surf_rate(a2) * SA
        ) / V for a2 in lams
    ]), 1.0)
    return cyl_int, cyl_srt, cyl_dct, cyl_crv, cyl_crv1


sr = (3 * V / 4 / np.pi)**(1/3)

# cylinder 1 is tall, cylinder 2 is short
a1r = 3.0
c1r = (V / (2 * a1r * np.pi))**(1/3)
c1h = 2 * a1r * c1r
a2r = 1.0
c2r = (V / (2 * a2r * np.pi))**(1/3)
c2h = 2 * a2r * c2r


sph_int = np.array([sph_intgtd(a2, sr) for a2 in lams])
sph_srt = np.minimum(np.array([surf_rate(a2) * 3.0 / sr for a2 in lams]), 1.0)
sph_dct = np.minimum(np.array([
    surf_rate(a2, depth_corr=sr) * 3.0 / sr for a2 in lams
]), 1.0)
cyl1_int, cyl1_srt, cyl1_dct, cyl1_crv, cyl1_crv1 = cylappx(c1r, c1h)
cyl2_int, cyl2_srt, cyl2_dct, cyl2_crv, cyl2_crv1 = cylappx(c2r, c2h)


fig1, ax1 = plt.subplots(1, 2, figsize=sidebysidefigsize)
ax1[0].plot(lams, sph_int, 'b-', label="Sphere")
ax1[0].plot(lams, sph_srt, 'k-', label="SA Approximation")
ax1[0].plot(lams, sph_dct, 'g-', label="Depth Corrected SA Appx")
ax1[0].set_xlabel(r'\( a^2 \)')
ax1[0].set_ylabel(r'\( \eta \)')
ax1[0].set_ylim((0.0, 1.0))
ax1[0].legend()
ax1[0].grid()
ax1[1].semilogy(lams, sph_srt - sph_int, 'k-', label="SA Approximation")
ax1[1].semilogy(lams, sph_dct - sph_int, 'g-', label="Depth Corrected SA Appx")
ax1[1].set_xlabel(r'\( a^2 \)')
ax1[1].set_ylabel('Difference')
ax1[1].yaxis.set_major_formatter(mtick.PercentFormatter(1.0))
ax1[1].grid()


fig2, ax2 = plt.subplots(1, 2, figsize=sidebysidefigsize)
ax2[0].plot(lams, cyl1_int, 'C1-', label=r"Cylinder (\( \gamma = 3.0 \))")
ax2[0].plot(lams, cyl1_srt, 'C2-', label="SA Approximation")
ax2[0].plot(lams, cyl1_dct, 'C3-', label="Depth Corrected SA Appx")
ax2[0].plot(lams, cyl1_crv, 'C4-', label="Curvature Approximation")
ax2[0].plot(lams, cyl1_crv1, 'C5-', label="Cylinder Curvature")
ax2[0].set_xlabel(r'\( a^2 \)')
ax2[0].set_ylabel(r'\( \eta \)')
ax2[0].set_ylim((0.0, 1.0))
ax2[0].legend()
ax2[0].grid()
ax2[1].semilogy(
    lams, cyl1_srt - cyl1_int, 'C2-', label="SA Approximation"
)
ax2[1].semilogy(
    lams, cyl1_dct - cyl1_int, 'C3-', label="Depth Corrected SA Appx"
)
ax2[1].semilogy(
    lams, cyl1_crv - cyl1_int, 'C4-', label="Curvature Approximation"
)
ax2[1].semilogy(
    lams, np.abs(cyl1_crv1 - cyl1_int), 'C5-', label="Cylinder Curvature"
)
ax2[1].set_xlabel(r'\( a^2 \)')
ax2[1].set_ylabel('Difference')
ax2[1].yaxis.set_major_formatter(mtick.PercentFormatter(1.0))
ax2[1].grid()


fig3, ax3 = plt.subplots(1, 2, figsize=sidebysidefigsize)
ax3[0].plot(lams, cyl2_int, 'b-', label=r"Cylinder (\( \gamma = 1.0 \))")
ax3[0].plot(lams, cyl2_srt, 'k-', label="SA Approximation")
ax3[0].plot(lams, cyl2_dct, 'g-', label="Depth Corrected SA Appx")
ax3[0].plot(lams, cyl2_crv, 'r-', label="Curvature Approximation")
ax3[0].set_xlabel(r'\( a^2 \)')
ax3[0].set_ylabel(r'\( \eta \)')
ax3[0].set_ylim((0.0, 1.0))
ax3[0].legend()
ax3[0].grid()
ax3[1].semilogy(
    lams, cyl2_srt - cyl2_int, 'k-', label="SA Approximation"
)
ax3[1].semilogy(
    lams, cyl2_dct - cyl2_int, 'g-', label="Depth Corrected SA Appx"
)
ax3[1].semilogy(
    lams, cyl2_crv - cyl2_int, 'r-', label="Curvature Approximation"
)
ax3[1].set_xlabel(r'\( a^2 \)')
ax3[1].set_ylabel('Difference')
ax3[1].yaxis.set_major_formatter(mtick.PercentFormatter(1.0))
ax3[1].grid()


all_fig, all_ax = plt.subplots(1, 1, figsize=sidebysidefigsize)
all_ax.loglog(lams, sph_srt - sph_int, 'b:', label="Sphere SA")
all_ax.loglog(lams, sph_dct - sph_int, 'g:', label="Sphere Depth Corrected")
all_ax.loglog(lams, cyl1_srt - cyl1_int, 'b-', label="Tall Cylinder SA")
all_ax.loglog(
    lams, cyl1_dct - cyl1_int, 'g-', label="Tall Cylinder Depth Corrected"
)
all_ax.loglog(
    lams, cyl1_crv - cyl1_int, 'r-', label="Tall Cylinder Curvature"
)
all_ax.loglog(lams, cyl2_srt - cyl2_int, 'b--', label="Squat Cylinder SA")
all_ax.loglog(
    lams, cyl2_dct - cyl2_int, 'g--', label="Squat Cylinder Depth Corrected"
)
all_ax.loglog(
    lams, cyl2_crv - cyl2_int, 'r--', label="Squat Cylinder Curvature"
)
all_ax.yaxis.set_major_formatter(mtick.PercentFormatter(1.0))
all_ax.grid()
all_ax.legend(bbox_to_anchor=(1.04, 1), loc="upper left", borderaxespad=0)
all_ax.set_xlabel(r'\( a^2 \)')
all_ax.set_ylabel("Error")


print((cyl2_crv - cyl2_int)[-4:])


for fig, fn in [
    (fig1, 'SA_appx_sph'),
    (fig2, 'SA_appx_cyl_tall'),
    (fig3, 'SA_appx_cyl_short'),
    (all_fig, 'SA_appx_all')
]:
    for ext in ['pdf', 'svg']:
        fig.tight_layout()
        fig.savefig("img/surfappx_{}.{}".format(fn, ext))


