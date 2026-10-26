# GEDI waveform simulator validation

The focused simulator in `R/simulator.R` was tested against Steven Hancock's
original `gedisimulator` source at commit
`49ad4f26b03083b61ee0bd424e8a7e4a8f2c301d` (2026-05-27). A standalone C
fixture generator in `tools/hancock_reference_fixture.c` directly executes the
supported equations and C float operations from `setGediPulse()`,
`waveFromPointCloud()`, and `applyPulseShape()` in upstream `gediIO.c`.

| Upstream behavior | rGEDI implementation |
|---|---|
| Gaussian footprint, `fSigma = 5.5` m | Gaussian return weights with the same default |
| Ignore footprint weights below 0.0006 | `intensity_threshold = 0.0006` |
| 15 ns pulse FWHM | Converted to a 0.955 m Gaussian sigma when `pSigma < 0` |
| Pulse sampled at waveform resolution unless `-pRes` is supplied | `pRes = NULL` uses `res`; an explicit `pRes` is also supported |
| Bin first, then convolve with the pulse | Same pulse support, bin offsets, edge divisor, and pulse-after ordering |
| Normalize for spatially varying ALS return density | Same 1.5 m default grid, origin, C integer conversion, and last-return counts |
| Separate classified ground returns | LAS class 2 returns are written to `grxwaveform` |
| Normalize waveform integral | The clean waveform integrates to one at the requested resolution |
| Ground-relative RH metrics | Uses classified ground when available and a lowest-peak fallback otherwise |
| Reflectance-corrected canopy cover | Uses canopy and ground energies with configurable reflectance ratio |
| Add detector noise | Optional relative Gaussian noise via `noise` |

The fixture contains a deterministic multi-return point cloud and tests count,
intensity, and fractional-return weighting, each with and without last-return
density normalization. Every waveform bin, every classified-ground bin,
elevation, and total energy is compared. The R and C arrays agree within
`5e-7`; this tolerance covers double-versus-float arithmetic while preserving
the C program's bin decisions.

The R implementation is deliberately smaller than the upstream command-line
simulator. It implements the symmetric Gaussian, discrete-return path needed by
rGEDI to create full waveforms and derive RH, canopy-cover, energy,
ground-height, and foliage-height-diversity metrics. `gediWFMetrics()` exposes
a focused set of direct waveform metrics; it is not a port of every Gaussian
fitting, denoising, detector, photon-counting, or LAI option in `gediMetric`.
Custom wavefront files, side-lobe models, and deconvolution are also outside
this focused implementation.

The bin accumulator uses a vectorized grouped sum rather than an R loop. This
improves the earlier portable rGEDI implementation, while the upstream native C
program remains the appropriate performance reference for very large ALS
collections. `maxBins` trims the simulated elevation window while preserving
the requested vertical sample spacing.

Validation is covered by `tests/testthat/test-sampling-modeling-simulator.R`.
The tests check the locked C reference arrays, HDF5 structure, sample order,
fixed waveform resolution, unit-integral normalization, classified-ground
recovery, RH ordering, canopy cover, noise behavior, and ASCII output. Rebuild
the C fixture with the command recorded at the top of its source file. Rebuild
the visual comparison with:

```r
source("readme/build-local-examples.R")
```

Upstream source: <https://bitbucket.org/StevenHancock/gedisimulator/src/master/>
