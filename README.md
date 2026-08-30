<!--
  Supplementary materials for
  "Instability analysis of massive static phantom wormholes via the spectral method"

  Authors:
    Dr. Davide Batic (Department of Mathematics, Khalifa University of Science
      and Technology, Abu Dhabi, UAE)
    Dr. Denys Dutykh (Department of Mathematics, Khalifa University of Science
      and Technology, Abu Dhabi, UAE)
-->

# Ellis–Bronnikov Wormholes: Quasi-Normal Modes by a Spectral Method

[![Eur. Phys. J. C](https://img.shields.io/badge/Eur.%20Phys.%20J.%20C-85%20%282025%29%20144-005a9c)](https://link.springer.com/article/10.1140/epjc/s10052-025-13867-x)
[![DOI](https://img.shields.io/badge/DOI-10.1140%2Fepjc%2Fs10052--025--13867--x-1a7f37)](https://doi.org/10.1140/epjc/s10052-025-13867-x)
[![arXiv](https://img.shields.io/badge/arXiv-2502.05486-b31b1b)](https://arxiv.org/abs/2502.05486)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-4c1)](LICENSE)
[![Maple](https://img.shields.io/badge/Maple-18%2B-ff6f00)](https://www.maplesoft.com/products/Maple/)
[![MATLAB](https://img.shields.io/badge/MATLAB-Advanpix%20MCT-0076a8)](https://www.advanpix.com/)
[![QNMs Hall of Fame](https://img.shields.io/badge/QNMs-Hall%20of%20Fame-8a2be2)](https://qnms.denys-dutykh.com/)
[![GitHub](https://img.shields.io/badge/GitHub-dutykh%2FEllisBronnikov-181717?logo=github)](https://github.com/dutykh/EllisBronnikov)
[![Last commit](https://img.shields.io/github/last-commit/dutykh/EllisBronnikov)](https://github.com/dutykh/EllisBronnikov/commits/main)

Symbolic derivation and arbitrary-precision numerical computation of the scalar
and axial quasi-normal mode (QNM) spectra of massive static phantom
(Ellis–Bronnikov) wormholes, by a Chebyshev collocation method combined with a
polynomial (quadratic) eigenvalue solver.

This repository holds the complete supplementary material for the article
[*Instability analysis of massive static phantom wormholes via the spectral
method*](https://doi.org/10.1140/epjc/s10052-025-13867-x), Eur. Phys. J. C
**85**, 144 (2025).

---

## Table of contents

- [Repository at a glance](#repository-at-a-glance)
- [Scientific context](#scientific-context)
- [The computational pipeline](#the-computational-pipeline)
- [The discrete eigenvalue problem](#the-discrete-eigenvalue-problem)
- [Coefficient functions and endpoint regularity](#coefficient-functions-and-endpoint-regularity)
- [Scalar versus axial perturbations](#scalar-versus-axial-perturbations)
- [Requirements](#requirements)
- [Quick start](#quick-start)
- [Parameter reference](#parameter-reference)
- [Interpreting and validating the output](#interpreting-and-validating-the-output)
- [Practical notes and troubleshooting](#practical-notes-and-troubleshooting)
- [QNMs Hall of Fame](#qnms-hall-of-fame)
- [Citation](#citation)
- [License](#license)
- [Authors](#authors)

---

## Repository at a glance

| File | Language | Role |
| --- | --- | --- |
| `SCALAR-PHANTOM-USING-ODE2.mw` | Maple worksheet | Full symbolic derivation for **scalar** perturbations: nondimensionalization, asymptotic factorization, compactification, extraction of the nine ODE coefficients, endpoint regularity checks. |
| `AXIAL-PHANTOM-USING-ODE2.mw` | Maple worksheet | The same derivation for **axial** (gravitational) perturbations. |
| `matrixassembler-scalar.mpl` | Maple source | Standalone procedure `MatrixAssembler` that collocates the scalar equation and exports the three matrices `M0`, `M1`, `M2`. |
| `matrixassembler-axial.mpl` | Maple source | The same for the axial case. |
| `chaseeigs.m` | MATLAB script | Reads the exported matrices in multiprecision, solves the quadratic eigenvalue problem, and overlays the spectra obtained at several resolutions. |
| `LICENSE` | Text | GNU General Public License, version 3. |

The two worksheets are the *derivation*; the two `.mpl` files are the
*production code* distilled from them, with the resulting coefficient functions
hard-coded so that no worksheet needs to be re-executed to reproduce the
spectra.

---

## Scientific context

An Ellis–Bronnikov wormhole is a static, spherically symmetric solution of the
Einstein equations sourced by a phantom scalar field, that is, a scalar with the
wrong sign of the kinetic term. The massive (Bronnikov) branch carries a
Schwarzschild-like mass in addition to the throat radius, and the single
dimensionless combination that controls its dynamics is the ratio of the
Schwarzschild radius to the throat radius, denoted `c` throughout the code.

Linearizing the field equations about this background and separating variables
reduces each perturbation sector to a single second-order ordinary differential
equation on the whole line, one asymptotically flat region on each side of the
throat. Quasi-normal modes are the frequencies `ω` for which that equation
admits a solution that is purely outgoing at both mouths. Because the resulting
boundary-value problem is not self-adjoint, the admissible frequencies are
complex, and the sign of `Im ω` decides stability: a mode with `Im ω > 0` grows
in time and signals an instability.

The published analysis finds purely imaginary quasi-normal modes that earlier
studies had missed. For scalar perturbations, instabilities appear once `c ⩾ 1`;
for axial perturbations, the onset occurs at smaller values of `c`, reflecting
the extra sensitivity of the geometry to gravitational-wave content. Wormholes
whose throat is much larger than the Schwarzschild radius remain stable.

Everything in this repository is written in dimensionless variables. Lengths are
measured in units of the throat radius `r0`, and frequencies in the reciprocal of
that unit, so `c` and `ω` are both pure numbers and the exponents appearing in
the asymptotic factorization are dimensionless by construction.

---

## The computational pipeline

The two Maple worksheets carry out the same five-stage derivation, and the
`.mpl` assemblers implement its outcome.

### Stage 1: nondimensionalization

The master radial equation `ODE2` for the perturbation amplitude `Z(r)` is
rewritten with `C = r0·c` and `r = r0·x` through `PDEtools[dchange]`, producing
`ODE2x` in the dimensionless radial coordinate `x = r/r0`, with `x ∈ (-∞, +∞)`
covering both asymptotic regions and `x = 0` at the throat.

### Stage 2: asymptotic factorization

The coefficients of `ODE2x` are expanded with `series(·, x = ±∞)` under the
assumption `c > 0`, which exposes the outgoing-wave behaviour at each mouth. That
behaviour is then factored out analytically,

```text
Z(x)  =  A(x; ω, c) · Φ(x)
```

where `A` collects the oscillatory and algebraic asymptotics (built from
`arctan x` and an exponential in `ω` and `c`). This is the decisive step: the
quasi-normal boundary conditions, which are a radiation condition rather than a
decay condition, become *properties of the ansatz* instead of constraints that
would have to be imposed numerically. The remaining amplitude `Φ` is bounded on
the whole line, and the equation it satisfies is called `ODEx`.

### Stage 3: compactification

The algebraic map

```text
x  =  tan(π y / 2),        y ∈ (-1, 1)
```

sends the infinite line onto the open unit interval, with the two mouths at
`y = ±1` and the throat at `y = 0`. Applying it to `ODEx` gives `ODE3`, an
equation posed on a bounded interval and therefore ready for a Chebyshev
expansion. The map is analytic and its inverse is `y = (2/π)·arctan x`.

### Stage 4: extraction and rescaling of the coefficients

Substituting the auxiliary symbols `Φ ↦ Z`, `Φ′ ↦ Y`, `Φ″ ↦ X` turns `ODE3` into a
linear form whose coefficients are read off with `coeff`, giving `SS0`, `SS1`,
`SS2`. Series expansions about `y = ±1` reveal the order of the zeros that these
coefficients develop at the endpoints; dividing by the appropriate power of
`(1 - y²)` yields the rescaled, endpoint-regular coefficients `M0`, `M1`, `M2`
and the equation

```text
M2(y) Φ″(y)  +  M1(y) Φ′(y)  +  M0(y) Φ(y)  =  0
```

Collecting this in the frequency splits each coefficient into a quadratic
polynomial, and the nine functions that the assemblers use are precisely the
resulting blocks:

```text
L0j  :  frequency-independent part      (multiplies Φ, Φ′, Φ″ for j = 0, 1, 2)
L1j  :  part linear in the frequency
L2j  :  part quadratic in the frequency
```

### Stage 5: endpoint regularity check

The worksheets close by evaluating the one-sided limits of all nine
coefficients at `y = ±1` with `limit(·, y = ∓1, right/left)`. All nine limits are
finite, which is what makes the scheme of the next section well posed. The
values are tabulated [below](#coefficient-functions-and-endpoint-regularity).

---

## The discrete eigenvalue problem

`MatrixAssembler` expands the unknown amplitude in the Chebyshev basis of the
first kind,

```text
                n-1
Φ(y)  ≈  P(y) =  Σ   a_j T_j (y)
                j=0
```

and enforces the differential equation at the `n` roots of `T_n`,

```text
y_i  =  cos( (2i - 1) π / (2n) ),        i = 1, …, n
```

These nodes are strictly interior, so no equation is ever evaluated at the
endpoints. Each collocation condition contributes one row, and the three
matrices are assembled entry by entry as the coefficients of the unknowns `a_j`:

```text
M0[i, j+1]  =  L00(y_i) T_j (y_i)  +  L01(y_i) T′_j (y_i)  +  L02(y_i) T″_j (y_i)
M1[i, j+1]  =  L10(y_i) T_j (y_i)  +  L11(y_i) T′_j (y_i)
M2[i, j+1]  =  L20(y_i) T_j (y_i)
```

for `i = 1, …, n` and `j = 0, …, n-1`. The terms in `L12`, `L21`, `L22` are
absent because those coefficients vanish identically.

The frequencies are then the eigenvalues of the quadratic (polynomial)
eigenvalue problem

```text
( M0  +  i ω M1  +  ω² M2 ) a  =  0
```

which `chaseeigs.m` solves with `polyeig(M0, mp('1i')*M1, M2)` in multiprecision
arithmetic. A problem of size `n × n` returns `2n` eigenvalues.

**No boundary rows are substituted into the matrices.** This is a deliberate
feature of the formulation rather than an omission. At `y = ±1` the coefficients
of `Φ″` and `Φ′` in the frequency-independent block, and the coefficient of `Φ″`
in the linear block, all vanish, while `L00`, `L11` and `L20` stay finite. The
equation therefore degenerates at each endpoint into an algebraic relation
between `Φ(±1)` and `Φ′(±1)`, a Robin-type condition that the expansion inherits
automatically. Combined with the factorization of Stage 2, this means the
quasi-normal boundary conditions are built into the operator itself.

---

## Coefficient functions and endpoint regularity

With `λ = L(L+1)` the separation constant of the angular problem, the one-sided
limits of the nine coefficients at the two mouths are as follows. They were
obtained symbolically and hold for both perturbation sectors.

| Coefficient | Multiplies | Limit as `y → -1⁺` | Limit as `y → +1⁻` |
| --- | --- | --- | --- |
| `L00` | `Φ` | `-π²λ/16` | `-π²λ/16` |
| `L01` | `Φ′` | `0` | `0` |
| `L02` | `Φ″` | `0` | `0` |
| `L10` | `Φ` | `0` | `0` |
| `L11` | `Φ′` | `-(π/4)·e^(πc)` | `π/4` |
| `L12` | `Φ″` | `0` (identically) | `0` (identically) |
| `L20` | `Φ` | `π c (π c e^(πc) + 2 e^(πc) + 1) e^(πc) / 8` | `π c (π c - e^(πc) - 2) / 8` |
| `L21` | `Φ′` | `0` (identically) | `0` (identically) |
| `L22` | `Φ″` | `0` (identically) | `0` (identically) |

Two consequences are worth recording. First, every limit is finite, so the
collocation rows stay bounded as the nodes crowd towards the endpoints and the
conditioning of the pencil degrades only through the usual `O(n⁴)` Chebyshev
differentiation growth. Second, `L11` and `L20` are asymmetric between the two
mouths, the `y → -1` limits carrying factors of `e^(πc)` that the `y → +1` limits
do not. That asymmetry is the fingerprint of the wormhole mass: it disappears
when `c → 0`, where the massless Ellis solution is symmetric about the throat.

---

## Scalar versus axial perturbations

The two assemblers differ in exactly **one line**. Eight of the nine coefficient
functions are shared; only `L00`, the frequency-independent coefficient of `Φ`,
which is where the effective potential enters, changes between the sectors:

```text
scalar:  L00(y) = -(1/4) cos²(πy/2) · [ (4 - c²) cos²(πy/2) + 2c sin(πy) + 4λ ]
                  / (1 - y²)²

axial:   L00(y) = -(1/4) cos²(πy/2) (y² - 1)²
                  · [ 3(c² - 4) cos²(πy/2) - 6c sin(πy) + 4λ ]
                  / (1 - y²)⁴
```

The angular momentum enters both sectors only through `λ = L(L+1)`, which is why
both expressions share the same endpoint limit `-π²λ/16`.

---

## Requirements

| Component | Version | Notes |
| --- | --- | --- |
| [Maple](https://www.maplesoft.com/products/Maple/) | 18 or newer | The worksheets are saved in the Maple 18 format. The `.mpl` assemblers use only `LinearAlgebra`, `ChebyshevT` and `ExportMatrix`, so they run on any modern release. |
| [MATLAB](https://www.mathworks.com/products/matlab.html) | R2015b or newer | Only `polyeig`, `legend` and standard plotting are used. |
| [Advanpix Multiprecision Computing Toolbox](https://www.advanpix.com/) | any | **Required.** `chaseeigs.m` depends on `mp.Digits`, `mp.read` and multiprecision `polyeig`. Double precision is far too coarse for this problem. |

Why multiprecision is not optional: the compactified operator concentrates its
resolution near the mouths, and the Chebyshev pencil at `n` in the low hundreds
has a condition number well beyond what 16 significant digits can absorb. In
double precision the physical eigenvalues are swamped by rounding noise long
before the spectrum converges.

---

## Quick start

The workflow has two stages: assemble in Maple, diagonalize in MATLAB.

### 1. Clone the repository and create the data directory

```bash
git clone https://github.com/dutykh/EllisBronnikov.git
cd EllisBronnikov
mkdir -p data
```

The `data` subdirectory is where Maple writes the matrices and where MATLAB
looks for them. It is not tracked, so it must exist before the first run.

### 2. Assemble the matrices in Maple

Start Maple in the repository root and load the assembler for the sector you
want:

```maple
read "matrixassembler-scalar.mpl":     # or "matrixassembler-axial.mpl"

MatrixAssembler(320, 200, 2, 0.7, "/absolute/path/to/EllisBronnikov"):
MatrixAssembler(320, 250, 2, 0.7, "/absolute/path/to/EllisBronnikov"):
MatrixAssembler(320, 300, 2, 0.7, "/absolute/path/to/EllisBronnikov"):
```

The arguments are, in order, the number of working digits `d`, the number of
Chebyshev modes `n`, the angular momentum `L`, the mass parameter `c`, and the
absolute path to the repository. Each call writes three ASCII files:

```text
data/M0_<n>.mat
data/M1_<n>.mat
data/M2_<n>.mat
```

Choose `d` at least as large as the precision you intend to use in MATLAB;
`chaseeigs.m` requests `n` digits, so `d ⩾ n` with a comfortable margin is the
safe setting.

### 3. Solve the eigenvalue problem in MATLAB

Open `chaseeigs.m` and set the toolbox path on the `addpath` line:

```matlab
addpath('/path/to/advanpix/toolbox/');
```

Then adjust the resolutions to match what Maple produced,

```matlab
list = [300, 250, 200];
```

start MATLAB in the repository root, and run:

```matlab
>> chaseeigs
```

The script reports its progress, then draws the three spectra on a single
complex plane, `Re ω` horizontally and `Im ω` vertically, with one marker style
per resolution.

---

## Parameter reference

### `MatrixAssembler(d, n, L, c, p)`

| Argument | Type | Meaning | Guidance |
| --- | --- | --- | --- |
| `d` | integer | Working precision, in decimal digits, for the whole assembly. | Set it at or above the MATLAB precision. Assembly cost grows with `d`, so do not inflate it beyond need. |
| `n` | integer | Number of Chebyshev modes, and hence the matrix size. | Run at least three values (for example 200, 250, 300) so that convergence can be assessed. |
| `L` | integer | Angular momentum of the multipole under study. | `L = 0` is admissible in the scalar sector; the axial sector starts at `L = 2`. |
| `c` | numeric | Dimensionless mass parameter, the ratio of the Schwarzschild radius to the throat radius. | The stability transitions of interest live near `c ≈ 1` in the scalar sector and below it in the axial sector. `c = 0` recovers the massless Ellis wormhole. |
| `p` | string | Absolute path to the directory that contains `data`. | Passed without a trailing slash; the procedure appends `/data/` itself. |

### Knobs in `chaseeigs.m`

| Variable | Meaning |
| --- | --- |
| `addpath(...)` | Location of the Advanpix toolbox. Must be edited before the first run. |
| `list` | Resolutions to load and compare. Every entry needs a matching triple of files in `data`. |
| `Marks` | Marker and colour specification, one cell per entry of `list`. The script stops with an explicit error if there are fewer markers than resolutions. |
| `mp.Digits(n)` | Working precision for the eigenvalue solve, tied by default to the mode count. |
| `axis([-6 6 -5 1])` | Viewing window in the complex frequency plane. Widen it when hunting for modes far from the origin. |

Inside `chaseeigs.m` the symbol `L` denotes the number of resolutions in `list`,
not the angular momentum; the angular momentum is fixed earlier, in Maple.

After the run, the eigenvectors, eigenvalues and condition numbers computed at
the highest resolution remain in the workspace as `V0`, `e0` and `k0` for
interactive inspection.

---

## Interpreting and validating the output

Spectral discretizations of non-self-adjoint operators return two populations of
eigenvalues mixed together: the genuine quasi-normal frequencies, and numerical
artefacts produced by the truncation. The overlay plot is the instrument that
separates them.

- **Converged modes** sit at the same point for every resolution in `list`. The
  markers land on top of one another, and refining `n` moves them by an amount
  that shrinks rapidly.
- **Spurious modes** drift visibly as `n` changes. They typically populate the
  outer parts of the plotted window and the far lower half-plane, and they should
  be discarded.
- **Stability verdict.** Modes with `Im ω > 0` grow in time. Their presence for a
  given pair `(L, c)` is the signature of instability. Because the window
  `axis([-6 6 -5 1])` includes a strip above the real axis, such modes are
  visible in the default view.
- **Purely imaginary modes**, the central finding of the paper, appear on the
  vertical axis at `Re ω = 0`. They are easy to overlook when only the
  oscillatory branches are examined, which is why the full complex window
  matters.

A practical convergence test: recompute at a larger `n`, or at the same `n` with
more digits `d`, and count how many leading significant digits of a candidate
frequency are stable. Physical modes gain digits as `n` grows; spurious ones do
not.

---

## Practical notes and troubleshooting

**`ExportMatrix` fails or produces nothing.** The `data` directory must exist
before `MatrixAssembler` is called, and the path `p` must be absolute and free of
a trailing slash. Maple does not create the directory for you.

**MATLAB cannot find the matrices.** `chaseeigs.m` uses the relative path
`data/M0_<n>.mat`, so MATLAB must be started in, or changed to, the repository
root. Every entry of `list` needs all three files.

**Why `mp.read` and not `load`.** Maple exports with
`ExportMatrix(..., target = MATLAB, mode = ascii)`, which writes the full decimal
expansion of every entry as text. Reading such a file with MATLAB's `load` would
truncate each entry to double precision and discard exactly the digits the method
depends on. `mp.read` parses the entire expansion.

**Assembly time.** `MatrixAssembler` differentiates the Chebyshev expansion
symbolically inside a doubly nested loop, so its cost grows roughly as the cube
of `n` and linearly in the digit count. At `n = 300` with several hundred digits
the assembly is measured in hours, not minutes. Assemble once, save the files,
and iterate on the MATLAB side.

**Error `Vector of markers is insufficient.`** `Marks` has fewer entries than
`list`. Add one marker specification per resolution.

**Changing `L` or `c`.** Both are baked into the exported matrices. Any change
requires re-running `MatrixAssembler`, and it is worth keeping separate `data`
directories per parameter pair to avoid overwriting a previous sweep, since the
file names encode only `n`.

---

## QNMs Hall of Fame

This work belongs to a longer series on the quasi-normal ringing of black holes
and wormholes, gathered at the
**[Quasi-Normal Modes Hall of Fame](https://qnms.denys-dutykh.com/)**.

Subtitled *Celestial Music*, the site is a curated showcase of the peer-reviewed
manuscripts produced by Davide Batic, Denys Dutykh and their collaborators across
the UAE, the Czech Republic and Italy. Each entry links the published article,
its arXiv preprint, and the companion repository of codes, so that the analytical
derivation, the numerical implementation and the resulting spectra can be
inspected side by side. The present article and this repository are among the
manuscripts listed there.

The framing is deliberate: a compact object perturbed and left to relax rings
down at a discrete set of complex frequencies determined by nothing but its own
geometry, the closest thing a spacetime has to a characteristic timbre. The Hall
of Fame collects those spectra, one geometry at a time.

---

## Citation

If these codes or the method they implement contribute to your research, please
cite the article:

> D. Batic and D. Dutykh, *Instability analysis of massive static phantom
> wormholes via the spectral method*, **Eur. Phys. J. C 85**, 144 (2025).
> [doi:10.1140/epjc/s10052-025-13867-x](https://doi.org/10.1140/epjc/s10052-025-13867-x)

```bibtex
@article{Batic2025Instability,
  author        = {Batic, Davide and Dutykh, Denys},
  title         = {Instability analysis of massive static phantom wormholes
                   via the spectral method},
  journal       = {The European Physical Journal C},
  volume        = {85},
  number        = {2},
  pages         = {144},
  year          = {2025},
  doi           = {10.1140/epjc/s10052-025-13867-x},
  eprint        = {2502.05486},
  archivePrefix = {arXiv},
  primaryClass  = {gr-qc}
}
```

**Publication record**

| Stage | Date |
| --- | --- |
| Received | 25 November 2024 |
| Accepted | 24 January 2025 |
| Published | 05 February 2025 |
| arXiv preprint | 08 February 2025 (`arXiv:2502.05486` [gr-qc]) |

Links: [published article](https://link.springer.com/article/10.1140/epjc/s10052-025-13867-x)
· [DOI](https://doi.org/10.1140/epjc/s10052-025-13867-x)
· [arXiv:2502.05486](https://arxiv.org/abs/2502.05486)

---

## License

Released under the [GNU General Public License, version 3](LICENSE). You are free
to use, modify and redistribute this material under the terms of that license.
Academic use additionally calls for the citation above.

---

## Authors

**Dr. Davide Batic**
Department of Mathematics, Khalifa University of Science and Technology,
Abu Dhabi, United Arab Emirates

**Dr. Denys Dutykh**
Department of Mathematics, Khalifa University of Science and Technology,
Abu Dhabi, United Arab Emirates
[Quasi-Normal Modes Hall of Fame](https://qnms.denys-dutykh.com/)

Repository: [github.com/dutykh/EllisBronnikov](https://github.com/dutykh/EllisBronnikov)
