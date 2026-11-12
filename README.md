# Metamaterial Workbench

**Band structures, band gaps and edge waves of locally resonant metamaterials, in your browser.**
No installation, no server, no licence fee. Light and dark themes. the whole tool is plain JavaScript that runs locally in any modern browser (or in Node.js).

▶ **Use it online:** https://guruvasawankumar-rgb.github.io/metamaterial-workbench/
⬇ **Use it offline:** download [`metamaterial-workbench.html`](metamaterial-workbench.html) and open it in a browser. It is a single self-contained file.
📘 **User manual:** [`docs/User_Manual.docx`](docs/User_Manual.docx) (40 pages, step by step, with screenshots)

Author: **Guruva Sawan Kumar**, Material and Structural Dynamics Laboratory, Sapienza University of Rome, in the group of **Prof. Walter Lacarbonara** · Licence: [MIT](LICENSE) · How to cite: see [below](#how-to-cite)

---

## What it does

Define a unit cell, a lattice plus a resonator, and the workbench computes how elastic waves travel through the periodic structure.

**Geometry**
- Ready-made hexagonal cells: (A) cantilever with tip mass, (B) centre mass on six spokes with an inner ring and six blocks ("spider web"), (C) beam across with a centre mass.
- *Build your own*: hexagon, square, triangle, kagome or free lattice, plus a resonator chosen from a list (cantilever with tip mass, centre mass on spokes, beam across, membrane with centre mass, mass on spring).
- Lattice vectors fitted automatically to a drawn frame (square, rectangle, parallelogram or hexagon), and kept in step when the frame is resized.
- A cell editor for anything else. You can add standard shapes (any regular polygon, rectangle or circle) as plates, beam frames, shell walls or holes. You can also add free beams, walls, plates and masses. Parts can be dragged with the mouse or have their numbers typed exactly.
- Separate materials for the host, the resonator and a membrane. Each part's thickness can be centred on, above or below its joining plane.

**Analyses**
- Bloch band structure (beams: Timoshenko or Euler–Bernoulli; plates and walls: flat shell elements). Out-of-plane, in-plane and complete band gaps, with mode shapes.
- Complex band structure (attenuation inside the gaps).
- Acoustic and optical branches labelled on every dispersion plot. Acoustic branches start at 0 Hz at Γ, where the whole cell moves together. Optical branches start above zero, where parts of the cell move against each other.
- Every gap labelled: which resonator mode opens it (local resonance) or whether it comes from the periodicity of the host (Bragg), and the bands at its edges.
- Zoom on the frequency axis of the band structure and complex band (scroll, drag, double-click to reset) to inspect narrow gaps.
- Animated mode shapes over any number of cells (1 × 1 to 15 × 15, square or M × N), with motion arrows: red for out-of-plane motion (up/down), blue for in-plane motion (sideways).
- A written summary after every run: each gap's range, width and type (local resonance or Bragg), and its peak attenuation.
- Sensitivity sweeps: resonator or host mass and stiffness, cell size, or any dimension of any part.
- Edge waves:
  - a strip periodic along one edge;
  - a finite M × N sample with free edges (parallelogram or staggered rectangle);
  - a driven response with a shaker on an edge, showing where the energy goes.
- Topology check: the Chern number (zero for passive structures) and the Zak phase of the bands below each gap.
- Zone path chosen automatically from the cell's symmetry, or entered as your own points.

## Quick start

1. Open the online version, or the downloaded `metamaterial-workbench.html`.
2. Under **1 Geometry**, choose a ready-made cell, or *Build your own*.
3. Set the material under **2 Material**, then press **Run analysis**.
4. Read the gaps in the **Dispersion** tab. Then try **Sensitivity** and **Edge waves**.

From the command line (Node.js 16 or later, nothing else to install):

    node tests/validate.js          # analytical checks of the beam engine
    node tests/validate_homog.js    # homogenisation, shells, drawn cells, plane waves, solvers
    node examples/run_analysis.js examples/design.json results    # bands, complex band, gaps -> CSV

After editing anything in `src/`, rebuild the single-file pages with `python3 tools/build_single.py`.

## Validation

The numerical engine is checked automatically against:
- analytical beam and plate solutions;
- patch tests of the shell element;
- the composite-beam stiffness of offset (stacked) parts;
- a dense solver, for the banded strip and finite-sample solvers;
- the exact half-cell shift of the Zak phase (see `tests/`).

Homogenised-plate results agree with full Bloch finite elements in the long-wave limit.

**Experimental validation.** The predictions were compared by the author with measurements on 3D-printed PLA samples (hexagonal beam–mass and spider-web resonator lattices):
- Band gaps and edge-wave frequency ranges agree once the stiffness of the printed material is calibrated.
- Details: for the multi-bandgap honeycomb lattice see *"A multi-bandgap metamaterial with multi-frequency resonators"*, Journal of Composite Materials 57(4), 783–804 (2023), [doi:10.1177/00219983231151578](https://doi.org/10.1177/00219983231151578).

## How to cite

If you use the workbench in your work, please cite the software (the DOI is created by Zenodo; see [CITATION.cff](CITATION.cff)):

> Guruva, S. K. (2026). *Metamaterial Workbench* (version 1.1.0) [Computer software]. https://doi.org/10.5281/zenodo.XXXXXXX

## How it was made

The physics, the requirements, the test cases and the experimental validation are the author's. Most of the code was written with the help of an AI assistant (Claude, by Anthropic), directed and checked by the author against analytical solutions, earlier numerical studies and measurements. Please report any error you find through the Issues tab.

## Acknowledgements

The workbench was conceived and developed independently by the author, in his own time, out of the difficulties he met during his PhD research. It builds on his PhD thesis, *Experimental Investigations Into Analysis of Wave Propagation in Acoustic Metamaterials* (Sapienza University of Rome, 2026), supervised by **Prof. Walter Lacarbonara**; the measurements from that thesis were used to validate the tool.

## Contributing

Questions, bug reports and suggestions are welcome as GitHub Issues. Pull requests should keep `node tests/validate.js` and `node tests/validate_homog.js` passing.

## Licence

MIT: free to use, change and share, including commercially, as long as the copyright notice is kept. See [LICENSE](LICENSE).
