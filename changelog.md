# Changelog

## [Unreleased]

### Added

### Changed

### Removed

### Fixed

## [1.3.0] - 2024-11-10

### Changed
- Code cleanup and minor optimizations
- Re-import select chemical metadata from Common Chemistry and PubChem. This includes being more consistent about proper capitalization of chemical synonyms. No work was done on the main massive pubchem database (large/small files) at this time.
- Move some vapor pressure extrapolation logic from thermo to proper functions in chemicals
- Add Hansen solubility parameter lookup functions and data files
- Fix issue with sqlite and multithreading https://github.com/CalebBell/thermo/issues/156
- Improve DIPPR equations integral over T calculation mode to use magic number -10 for calculation order instead of -1j which seemed clever at the time but was really a horrible idea and has caused lots of grief with pypy, and numba over the years
- Improve Github Actions to test more platforms
- Fluids version dependency now >= 1.0.27

## [1.2.0] - 2024-07-26

### Changed
- Compatibility with NumPy 2.0 and SciPy 1.14
- Fluids version dependency now >= 1.0.26

## [1.1.4] - 2023-06-04

### Changed
- Code cleanup with ruff (experiment)
- Add {critical, triple, boiling} data from REFPROP and other high-precision EOSs. This is the preferred method where available.