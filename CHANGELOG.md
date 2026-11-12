# Changelog

## 1.1.0 (2026-11-11)
- Mode shapes: choose how many cells to show (1 × 1 to 15 × 15, or any M × N).
- Motion arrows on mode shapes: red for out-of-plane (up/down), blue for in-plane (sideways).
- Acoustic and optical branches labelled on the dispersion curves and named in the mode description.
- Gap labels: resonator mode that opens each gap (or Bragg), bands at its edges; "Opened by" column in the gap table.
- Zoom and pan on the frequency axis of the band structure and complex band.
- "Fit lattice to frame shape": a₁ and a₂ worked out from a square, rectangle, parallelogram or hexagon frame, optionally kept in step when the frame changes.
- The first frame shape drawn in an empty cell sets the lattice automatically; a warning with a one-click fix appears when frame and lattice do not match.
- Branch labels explained in plain words (acoustic: whole cell moves together; optical: parts move against each other).
- New look: light GitHub-style panels, colour-blind-safe palette, light/dark theme switch.

## 1.0.0 (2026-11-11)
First public release.

- Ready-made hexagonal cells (cantilever with tip mass, spider web, beam across) and *Build your own* (lattice + resonator).
- Cell editor with standard shapes, mouse and numeric editing, roles, materials and thickness reference.
- Bloch bands, complex bands, written gap summary, sensitivity sweeps, zone path from symmetry or user-defined.
- Edge waves: strip, finite M × N sample (parallelogram or staggered rectangle), driven response with a shaker.
- Topology check: Chern number and Zak phase.
- User manual (docs/User_Manual.docx).
