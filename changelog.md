# Changelog

## [Unreleased]

### Added

### Changed

### Removed

### Fixed

## [1.3.3] - 2025-03-16

### Added

- Added Baumé scale conversion function: 8 new functions Baume_heavy_to_SG, SG_to_Baume_heavy, Baume_light_to_SG, SG_to_Baume_light, Baume_heavy_to_rho, rho_to_Baume_heavy, Baume_light_to_rho, rho_to_Baume_light
- Handle division by zero at Wagner T == Tc condition - return mathmatically correct result, infinity
- Updated several Joback predictions following correction of one SMARTS pattern
- Fix four duplicate/inconsistent entries in Zabransky which was previously being handled accidentally by overwritting data due to duplicates in the original source material - now the better data is intentionally loaded

### Changed
- Chemical Metadata work


## [1.3.2] - 2024-12-08

### Changed
- Chemical Metadata work
- Improved performance of T_converter and viscosity_converter, which previously build splines on the fly; they now have already-built conversion routines


## [1.3.1] - 2024-11-10

### Changed
- Added Dipalmitin back as synonym for 502-52-3

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