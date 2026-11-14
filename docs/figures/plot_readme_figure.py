#!/usr/bin/env python3
# Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science and Technology, Abu Dhabi, UAE)
"""Figure of the README: the effective potentials of the massive Ellis-Bronnikov wormhole,

    r0^2 V_S(x) = e^(2f) [ (l(l + 1) + 1)/(1 + x^2) - (2x - c)^2/(4(1 + x^2)^2) ]      (scalar)
    r0^2 V_A(x) = e^(2f) [ (l(l + 1) - 3)/(1 + x^2) + 3(2x - c)^2/(4(1 + x^2)^2) ]    (axial)
    f(x) = c (arctan x - pi/2),     x = r/r0,     c = 2M/r0,

equations (4), (6), (34) and (36) of the article, against x on both sides of the throat. From this
folder,
    python3 plot_readme_figure.py
writes ellisbronnikov-potentials.png (Python 3 with NumPy and Matplotlib). The style is the stamped
copy manuscript_style.py beside it.
"""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import manuscript_style as hs  # noqa: E402

C_VALUES = (0.0, 0.25, 0.5, 1.0)


def redshift(x, c):
    return c * (np.arctan(x) - np.pi / 2)


def scalar(x, L, c):
    return np.exp(2 * redshift(x, c)) * ((L * (L + 1) + 1) / (1 + x**2) - (2 * x - c) ** 2 / (4 * (1 + x**2) ** 2))


def axial(x, L, c):
    return np.exp(2 * redshift(x, c)) * ((L * (L + 1) - 3) / (1 + x**2) + 3 * (2 * x - c) ** 2 / (4 * (1 + x**2) ** 2))


hs.configure(family="libertine")
plt.rcParams["savefig.dpi"] = 260
plt.rcParams["figure.facecolor"] = hs.PAPER
plt.rcParams["savefig.facecolor"] = hs.PAPER

fig, (left, right) = plt.subplots(1, 2, figsize=(hs.WIDTH, 2.7))
x = np.linspace(-6.0, 8.0, 2801)
for axis, potential, L, title in ((left, scalar, 0, r"(a) scalar, $\ell = 0$"),
                                  (right, axial, 2, r"(b) axial, $\ell = 2$")):
    for k, c in enumerate(C_VALUES):
        style = hs.series(k)
        style["marker"] = None
        values = potential(x, L, c)
        axis.plot(x, values, lw=1.1, label=rf"$\mathfrak{{c}} = {c:g}$", **style)
        print(f"{title}: c = {c:g}: min {values.min():.5f} at x = {x[values.argmin()]:.3f}, "
              f"max {values.max():.5f} at x = {x[values.argmax()]:.3f}")
    axis.axhline(0.0, color=hs.SLATE, lw=0.6, linestyle=(0, (1.6, 1.6)), zorder=0)
    axis.axvline(0.0, color=hs.MIST, lw=0.6, zorder=0)
    axis.set_xlim(-6.0, 8.0)
    axis.set_xlabel(r"$x = r/r_0$")
    axis.set_ylabel(r"$r_0^2\,V(x)$")
    hs.panel(axis, title)
left.set_ylim(-0.03, 1.05)
right.set_ylim(-0.1, 3.2)
hs.legend_outside(fig, where="lower", ncol=4)

out = HERE / "ellisbronnikov-potentials.png"
hs.save(fig, out)
# The background is opaque: without its alpha channel, and compressed a second time without
# loss, the file is smaller.
Image.open(out).convert("RGB").save(out, optimize=True)
print("written", out)
