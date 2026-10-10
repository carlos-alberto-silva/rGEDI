# GEDI waveform simulator validation

The focused simulator in `R/simulator.R` was compared with Steven Hancock's
original `gedisimulator` source at commit
`49ad4f26b03083b61ee0bd424e8a7e4a8f2c301d` (2026-05-27). The comparison used
the upstream implementations in `gediIO.c`, `gediMetric.c`, and `gediNoise.c`.

| Upstream behavior | rGEDI implementation |
|---|---|
| Gaussian footprint, `fSigma = 5.5` m | Gaussian return weights with the same default |
| Ignore footprint weights below 0.0006 | `intensity_threshold = 0.0006` |
| 15 ns pulse FWHM | Converted to a 0.955 m Gaussian sigma when `pSigma < 0` |
| Bin first, then convolve with the pulse | Fixed-resolution vectorized binning followed by Gaussian convolution |
| Normalize for spatially varying ALS return density | Enabled by default with last-return density cells |
| Separate classified ground returns | LAS class 2 returns are written to `grxwaveform` |
| Normalize waveform integral | The clean waveform integrates to one at the requested resolution |
| Ground-relative RH metrics | Uses classified ground when available and a lowest-peak fallback otherwise |
| Reflectance-corrected canopy cover | Uses canopy and ground energies with configurable reflectance ratio |
| Add detector noise | Optional relative Gaussian noise via `noise` |

The R implementation is deliberately smaller than the upstream command-line
simulator. It implements the path needed by rGEDI to create full waveforms and
derive RH, canopy-cover, energy, ground-height, and foliage-height-diversity
metrics. It does not implement custom wavefront files, side-lobe models,
deconvolution, or the upstream detector and photon-counting noise modes.

The bin accumulator uses a vectorized grouped sum rather than an R loop. This
improves the earlier portable rGEDI implementation, while the upstream native C
program remains the appropriate performance reference for very large ALS
collections. `maxBins` trims the simulated elevation window while preserving
the requested vertical sample spacing.

Validation is covered by `tests/testthat/test-sampling-modeling-simulator.R`.
The tests check HDF5 structure, sample order, fixed waveform resolution,
unit-integral normalization, classified-ground recovery, RH ordering, canopy
cover, noise behavior, and ASCII output. Rebuild the visual comparison with:

```r
source("readme/build-local-examples.R")
```

Upstream source: <https://bitbucket.org/StevenHancock/gedisimulator/src/master/>
