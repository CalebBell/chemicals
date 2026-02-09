# Changelog

## [1.5.1] - 2026-02-08

### Added

- New `COSTALD_mixture_compressed` function for compressed liquid mixture density using the COSTALD-Tait CSP method (API Procedure 6A3.4)
- New `COSTALD_compressed` function for compressed liquid density of pure components

### Changed

- CI: Skip numba tests on macOS Intel (no llvmlite wheels available)
- CI: Enable numba tests on macOS ARM
- CI: Remove `zoneinfo` from cx_Freeze excludes (required by pandas)
- CI: Remove unnecessary numba dependency from Nuitka test recipe

### Removed

### Fixed

- Fixed E coefficient in heat capacity table from Perry 8E for 1-heptene
- Corrected spelling of "Demmler" viscosity cups (previously misspelled as "Demmier"). Added citation to Hydraulic Institute's Fluid Properties new website which shows the correct spelling.

## [1.5.0] - 2025-10-26

### Added

- Project is now PEP 517 compliant and doesn't use deprecated setup.py commands anymore
- GitHub Actions workflow for publishing to PyPI using environment protection
- `uv` package manager integration across CI workflows for faster dependency resolution
- New consolidated `quality.yml` workflow for linting and testing
- New `build_third_party_packagers.yml` workflow consolidating cx_Freeze, PyInstaller, Nuitka, and py2exe testing
- Pre-commit hooks configuration
- Justfile with extensive development automation commands - github actions have been refactored to use this where possible, making them locally debuggable

### Changed

- **Breaking**: Dropped Python 3.8 support; minimum Python version is now 3.9
- Migrated from setup.py to pyproject.toml-only configuration (PEP 517)
- Consolidated multiarch CI workflows with reduced test matrix for improved performance
- Optimized multi-numpy/scipy testing workflow using `uv`
- Simplified third-party packager testing into single consolidated workflow
- Moved coverage configuration to pyproject.toml
- Moved pytest configuration from pytest.ini to pyproject.toml
- Moved mypy configuration from mypy.ini to pyproject.toml
- Moved Ruff configuration from .ruff.toml to pyproject.toml
- Reorganized development requirements into pyproject.toml optional dependencies
- Updated PyInstaller spec file configuration

### Removed

- Removed setup.py (replaced by pyproject.toml)
- Removed standalone configuration files: pytest.ini, mypy.ini, .ruff.toml
- Removed separate workflow files: build_cxfreeze_library.yml, build_nuitka_library.yml, build_py2exe_library.yml, build_pyinstaller_library.yml
- Removed pre-commit.yml and security.yml workflows (consolidated into quality.yml)
- Removed separate requirements files (test, docs, multiarch) - now in pyproject.toml

### Fixed

- Fixed test_identifiers.py compatibility issues
- Fixed test_numba.py compatibility issues

### Security

- Implemented PyPI publishing workflow with manual approval gate

## [1.4.0] - 2025-10-19

### Added

- Python 3.13 and 3.13t (free-threaded) support with PYTHON_GIL=0 configuration
- Pre-commit configuration with Ruff, mdformat, and file validators
- New GitHub Actions workflows for pre-commit checks and security scanning
- Packaging compatibility workflows for cx_Freeze, PyInstaller, and py2exe
- Standalone test scripts and demo builders for verifying packaged distributions
- Coverage HTML artifact uploads to all test workflows
- Concurrency controls to workflows to cancel redundant builds
- Justfile for streamlined development tasks (setup, docs, test, typecheck, lint)
- Security scanning with pip-audit and bandit

### Changed

- Minimum Python version raised from 3.6 to 3.8
- Updated actions to latest versions (setup-qemu v3, run-on-arch v3)
- Updated macOS CI runners (macos-13 → macos-15-intel, added macos-latest for ARM)
- Extensive code quality improvements with Ruff linting across entire codebase:
  - String quote normalization to double quotes
  - Removed unused imports and variables
  - Improved code formatting and PEP 8 compliance
  - Better type hints compatibility
- Merged type hints across the codebase with improved accuracy
- Renamed workflow files to use .yml extension consistently
- Updated copyright year to 2025
- Fixed numerous typos across documentation files
- Improved Sphinx configuration for Python 3.13 compatibility
- Enhanced docstring and markdown formatting
- Updated README to reflect Python 3.8+ requirement

### Removed

- Dropped Python 3.6 and 3.7 support
- Removed obsolete platform-specific exclusions

### Security

- Added automated security scanning workflow documented in SECURITY.md
- Added .claude to .gitignore

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
