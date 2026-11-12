# Technical reference: code structure and computation pipeline

Two workbenches share this code:

* `metamaterial-workbench.html` (from `index.html`): the **FE workbench**, Bloch finite elements only (parametric beam lattices and drawn cells with beams, shell walls and plates).

Rebuild both single files with `python3 tools/build_single.py`.


Band-structure workbench for locally resonant beam lattices: Bloch bands, complete and full band
gaps, complex band (phase and attenuation), mass/stiffness sensitivity and edge waves, for
out-of-plane (flexural) and in-plane waves, with Timoshenko or Euler–Bernoulli beams.

There is no server and no external library. Every computation is plain JavaScript that runs in
the browser (or in Node.js), and nothing is loaded from the internet.

## Files

| File | What it is |
|---|---|
| `metamaterial-workbench.html` | The complete workbench as one self-contained file (what you use day to day). |
| `index.html` | The same page, loading the separate files below. Open it in a browser; identical behaviour. |
| `src/engine.js` | **The computational engine.** All numerics: geometry, finite elements, Bloch solver, eigen-solvers, gap detection, complex band, edge-wave strip. No browser code. Works in Node.js. |
| `src/app.js` | The interface: reads the inputs, calls the engine, runs the analyses, draws the plots, writes the assessment text. |
| `src/style.css` | Page styling. |
| `examples/run_analysis.js` | Command-line runner: the same computation without a browser, results written as CSV. |
| `examples/design.json` | Your sample (14 mm hexagon, block resonator) as an input file for the runner. |
| `src/homog.js` | Homogenised plate: static energy homogenisation of any periodic 3D beam frame (A, B, D, D_free), resonator modal reduction, plane-wave expansion, full 3D Bloch FE of frames. |
| `src/viz.js` | Shared solve (FEM / PWE / both), band and complex-band plots, gap table, animated mode shapes. |
| `src/draw.js` | Drawn cells: meshing (beams, shell walls, shell plates with holes), flat shell element (CST + DKT + drilling), welding and periodicity, shift-invert Bloch solver, homogenisation and resonators of drawn cells. |
| `src/drawapp.js` | FE workbench cell editor (cell type “Drawn cell”): drawing tools, hexagonal templates, 3D view, FE results. |
| `tests/validate_homog.js` | 24 checks for homogenisation, complex bands, plane waves, shells and drawn cells (`node tests/validate_homog.js`). |
| `tools/build_single.py` | Rebuilds the single offline HTML from the source files. |
| `tests/validate.js` | Validation against analytical solutions (13 checks). |

## Running without the browser

Needs Node.js (any version from 16 on), no other installation.

    node tests/validate.js                                    # checks the engine, about 5 s
    node examples/run_analysis.js examples/design.json out    # bands, complex band, gaps -> out/*.csv

`run_analysis.js` writes, for out-of-plane (`oop`) and in-plane (`ip`) waves:

* `bands_<fam>.csv` frequency of every band at every k-point on the zone path, and the share of
  kinetic energy in the resonators for each point;
* `complex_<fam>.csv` every wave at every frequency along Γ–X: phase per cell Re(μ) and attenuation
  per cell Im(μ), in Np and dB;
* `gaps.csv` complete gaps of each polarisation with their type (local resonance or Bragg), and the
  full gap where both polarisations are stopped.

Edit `design.json` for other designs. For a cell the parameters cannot describe, export it from the
workbench ("Copy this design into the custom editor"), save the text as a file, and give its name as
`"custom_cell_file"` in `design.json`.

## Computation pipeline

    inputs (mm, g, GPa)
      -> baseCell + addResonators          cell specification: nodes, beams, sections, masses, springs
      -> buildModel                        finite-element model of one unit cell, per polarisation
      -> kPath                             wave vectors along the irreducible Brillouin-zone boundary
      -> solveK at every k                 Bloch eigenproblem [K(k) - w^2 M(k)] u = 0
      -> findGaps                          complete gaps; intersect oop and ip -> full gap
      -> buildWFE + wfeAt                  complex band: propagation constants at fixed frequency
      -> ribbonSpec + solveK               edge waves: strip N cells wide, free edges
      -> sensitivity: rebuild with scaled mass/stiffness, central differences of gap edges

### 1. Unit cell (`baseCell`, `addResonators`, `specFromJSON`)

A cell is a list of nodes (x, y), beams between nodes (with the cell offset `n = [n1, n2]` of the
far node, so a beam can cross into the neighbouring cell), named cross-sections, point masses with
rotary inertia, and mass-on-spring resonators. `baseCell` builds square, triangular, kagome and
hexagon (regular, stretched or re-entrant) cells; `addResonators` adds arms (straight or L-shaped,
with a point mass or a real block at the tip), springs, or a central mass on ligaments, splitting
walls where an arm is attached. Inside the wall, the arm root is a stiff stub of wall height.

### 2. Finite elements (`buildModel`, `bendBlock`, `beamKM`, `beamKMip`)

Each beam is meshed into elements in proportion to its length. Two independent polarisations,
three degrees of freedom per node each:

* **out-of-plane** (`beamKM`): w and the two cross-section rotations; bending about the in-plane
  axis (I = b h³/12) plus St-Venant torsion (GJ, thin-rectangle series for J).
* **in-plane** (`beamKMip`): u, v and the rotation θz; axial bar (EA) plus bending in the plane
  (Iz = h b³/12).

Bending uses `bendBlock`, a **Timoshenko** element with interdependent interpolation (exact in
statics, no shear locking) and the consistent mass and rotary inertia of Friedman and Kosmatka:
Φ = 12EI/(κGAℓ²), κ is Cowper's shear coefficient. With `beam: 'eb'`, Φ = 0 and the element is
Euler–Bernoulli. Element matrices are rotated to global axes (`tAT`) and assembled.

### 3. Bloch–Floquet periodicity (`assembleBloch`, `addElem`)

A DOF in the cell displaced by R = n1·a1 + n2·a2 equals the reference DOF times e^(i k·R). Elements
crossing the cell boundary are assembled with that phase on their far node, giving the Hermitian
matrices K(k), M(k).

### 4. Eigen-solution (`solveK`, `tred2`, `tqli`, `tridiagInv`)

The complex Hermitian problem of size n is written as a real symmetric one of size 2n (each
eigenvalue appears twice), reduced to standard form with the Cholesky factor of M, tridiagonalised
by Householder reflections (`tred2`) and solved by implicit QL (`tqli`). Eigenvectors, when
needed for mode shapes or energy shares, come from inverse iteration on the tridiagonal matrix
(`tridiagInv`). The energy share of a mode is uᴴ M_res u / uᴴ M u over the resonator DOFs.

### 5. Zone path and gaps (`recip`, `bzPolygon`, `symmetryInfo`, `makePath`, `kPath`, `findGaps`)

The Brillouin zone is built from the reciprocal vectors (`bzPolygon`). `symmetryInfo` tests which
mirrors and rotations map the cell onto itself, and `makePath` picks Γ–X–K–Γ (hexagonal), Γ–X–M–Γ
(square), the quarter zone or the half zone accordingly. A **complete gap** is reported where the
maximum of band n along the path lies below the minimum of band n + 1. The **full gap** is the
intersection of the out-of-plane and in-plane complete gaps (done in `app.js`, `intersectGaps`).

### 6. Complex band (`buildWFE`, `wfeAt`, `eigNonsym`, `luSolve`)

Wave finite element method along a1. Beams crossing the a1 boundary are cut at mid-element; the
interior DOFs are condensed out of the dynamic stiffness D = K − ω²M (`luSolve`), and continuity
u_R = λ u_L with force balance gives the quadratic eigenproblem
(D_RL + λ(D_LL + D_RR) + λ² D_LR) u_L = 0. It is solved through λ = (1 + s)/(1 − s), which keeps it
regular when D_LR is rank-deficient, as a linear companion problem (`eigNonsym`: balancing,
Hessenberg reduction, shifted QR). with λ = e^(−iμ): Re(μ) is the phase per cell, Im(μ) the attenuation
per cell. Known blind spot: standing modes whose cut faces do not move (a wall vibrating between
stationary points) carry nothing across the boundary and appear only in the Bloch bands.

### 7. Edge waves (`ribbonSpec`, `nodeDepth`, `massDiag`, `edgeShares`)

`ribbonSpec` stacks N cells across the strip, keeps periodicity only along the edge, and leaves
both long sides free. Strip bands are computed against the edge wave number; bulk bands are
projected on the same wave number by sweeping the other component through a full period. A strip
mode is an edge wave when it lies outside every projected bulk band and has at least 60 % of its
kinetic energy (`edgeShares`) in the outer cells.

### 8. Sensitivity (`runSweep` in `app.js`)

Scales one quantity (resonator mass, resonator stiffness, host density, host modulus, cell size),
rebuilds the model, recomputes the gaps, and reports S = ∂ ln f_edge / ∂ ln p by central
differences at the nominal design.

## Engine API (src/engine.js)

| Function | Inputs | Returns |
|---|---|---|
| `baseCell(G)` | `{lattice, a, l, h, theta}` in m, deg | cell spec without resonators |
| `addResonators(base, R)` | resonator description (see `examples/run_analysis.js`) | cell spec with resonators |
| `specFromJSON(text)` / `specToJSON(spec)` | custom-cell text (mm, g) | cell spec / text |
| `buildModel(spec, mat, opt)` | `mat = {E, nu, rho}` SI; `opt = {ne, path, family: 'oop'|'ip', beam: 'tim'|'eb', scale}` | model: matrices, lattice, masses |
| `kPath(md.Lat, nper)` | k-points per segment | `{ks, xs, ticks, labels, seg}` |
| `solveK(md, kx, ky, nb, wantVec, wantFrac)` | wave vector in rad/m | `{f (Hz), vecs, frac}` |
| `findGaps(F, nb)` | `F[k][band]` in Hz | `[{lo, hi, below}]` |
| `buildWFE(md)`, `wfeAt(W, f)` | frequency in Hz | `[{att (Np/cell), ph (rad/cell)}]` |
| `ribbonSpec(spec, N, 'a1'|'a2')` | strip width in cells | strip spec |
| `sectionProps(sec)`, `shearCoef(shape, nu)`, `bendBlock(...)` | section, element data | properties, element matrices |

All engine quantities are SI (m, kg, Pa, Hz); the interface and the runner convert from mm and g.

## Assumptions and limits

* Members are beams along their centrelines (Timoshenko: shear and rotary inertia included; no
  warping or distortion of the section). Joint regions and the solid blocks are represented by
  centreline beams with their real section and mass.
* Joints are rigid. Material is linear elastic, without damping.
* Out-of-plane and in-plane waves are independent only for a flat lattice symmetric about its
  mid-plane.
* Gaps are checked along the boundary of the irreducible zone, not inside it (the usual practice;
  the half-zone path can be forced for cells with low symmetry).
* For final geometry, confirm the gap edges with a 3D solid-element Bloch model (COMSOL, ANSYS);
  expect differences of a few percent.
