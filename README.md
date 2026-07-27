# Gas-Coupled Laser Acoustic Detection Simulation

MATLAB simulations for studying acoustic wave propagation in air and directional sensing concepts inspired by gas-coupled laser acoustic detection (GCLAD).

This repository focuses on simulation and visualization workflows for a related but different goal than the source paper: exploring how laser-based, non-contact acoustic measurement ideas can be modeled in 2D and 3D using `k-Wave`.

## Project Summary

The code in this repo was built around a project on measuring and visualizing sound in 2D using a laser-based setup. The simulations here are used to:

- model acoustic propagation in air with `k-Wave`
- compare source geometries and directional response
- render 2D pressure fields as video
- inspect 3D fields through slice-based visualization
- test how runtime changes with grid size, spacing, and frequency

The work is inspired by the paper `Gas-coupled laser acoustic detection as a non-contact line detector for photoacoustic and ultrasound imaging`, which is included in `docs/references/`.

## Repository Layout

```text
.
├── analysis/                Parameter sweeps and convergence studies
├── docs/
│   └── references/          Paper and reference material
├── k-wave-toolbox-version-1.4/
├── legacy/                  Older experiments and superseded variants
├── notebooks/               MATLAB live scripts and exploratory notebooks
├── results/                 Suggested output location for generated media
├── src/                     Main simulation and visualization functions
├── LICENSE.txt
├── README.md
└── startup.m                Adds project folders to the MATLAB path
```

## Main Entry Points

- `src/simulate_plane.m`: baseline 2D line-source simulation with optional video output
- `src/simulate_plane_upscale.m`: 2D simulation with directional gain analysis and optional upscaled rendering
- `src/simulate_speakers_upscale.m`: two-source / rotated-source variant for directional studies
- `src/simulate_plane_3d.m`: 3D simulation that saves the pressure field to a `.mat` file
- `src/visualize_plane_3d.m`: animated orthogonal slice viewer for saved 3D data
- `src/visualize_isoslices.m`: static 3D slice visualization
- `analysis/convergence_analysis.m`: parameter sweep for runtime and scaling behavior
- `analysis/convergence_analysis2.m`: alternate sweep implementation using a combined parallel loop

## Setup

This project depends on the bundled `k-Wave` toolbox.

In MATLAB, from the repository root, run:

```matlab
startup
```

That adds the main project folders and `k-Wave` to your MATLAB path for the current session.

## Quick Start

Run a simple 2D simulation:

```matlab
startup
simulate_plane('render', false)
```

Render a 2D video:

```matlab
startup
simulate_plane('render', true, 'mp4_filename', 'results/plane_simulation.mp4')
```

Run directional analysis:

```matlab
startup
simulate_plane_upscale('gainAnalysis', true, 'plots', {'linear', 'dB', 'interp'})
```

Run the two-source variant:

```matlab
startup
simulate_speakers_upscale('src_spacing', 0.02, 'alpha', 15, 'gainAnalysis', true)
```

Run a 3D simulation and inspect it:

```matlab
startup
simulate_plane_3d(0.10, 1e-3, 40e3, 'results/kwave_3d_data.mat')
visualize_plane_3d('results/kwave_3d_data.mat')
```

Run a convergence study:

```matlab
startup
[results, stats] = convergence_analysis()
```

## Notes

- The main medium settings currently model an air-like homogeneous medium.
- Several files in `legacy/` were kept intentionally to show the evolution of the project.
- `notebooks/` contains exploratory live scripts rather than polished entry points.
- Output files such as `.mp4`, `.gif`, and `.mat` are best written into `results/`.

## Reference

Primary inspiration:

J. L. Johnson, K. van Wijk, J. N. Caron, and M. Timmerman, `Gas-coupled laser acoustic detection as a non-contact line detector for photoacoustic and ultrasound imaging`, Journal of Optics, 2016.

The PDF is included here for context in `docs/references/`.

## License

The repository includes the original `k-Wave` license and this project's root `LICENSE.txt`. Review both before redistribution.
