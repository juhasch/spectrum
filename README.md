# SPECTRUM: Spectral Analysis in Python

[![PyPI version](https://badge.fury.io/py/spectrum.svg)](https://pypi.python.org/pypi/spectrum)
[![CI](https://github.com/cokelaer/spectrum/actions/workflows/main.yml/badge.svg?branch=master)](https://github.com/cokelaer/spectrum/actions/workflows/main.yml)
[![Coverage Status](https://coveralls.io/repos/cokelaer/spectrum/badge.png?branch=master)](https://coveralls.io/r/cokelaer/spectrum?branch=master)
[![Conda-Forge License](https://anaconda.org/conda-forge/spectrum/badges/license.svg)](https://anaconda.org/conda-forge/spectrum)
[![Conda-Forge Version](https://anaconda.org/conda-forge/spectrum/badges/version.svg)](https://anaconda.org/conda-forge/spectrum/badges/version.svg)
[![Conda-Forge Downloads](https://anaconda.org/conda-forge/spectrum/badges/downloads.svg)](https://anaconda.org/conda-forge/spectrum)
[![JOSS](http://joss.theoj.org/papers/e4e34e78e4a670f2ca9a6a97ce9d3b8e/status.svg)](http://joss.theoj.org/papers/e4e34e78e4a670f2ca9a6a97ce9d3b8e)

<img src="http://www.thomas-cokelaer.info/software/spectrum/html/_images/psd_all.png" alt="Spectrum PSDs" width="50%" align="right" />

**Spectrum** contains tools to estimate Power Spectral Densities using methods based on Fourier transform, parametric methods, or eigenvalue analysis:

- **Fourier methods**: correlogram, periodogram, and Welch estimates. Standard tapering windows (Hann, Hamming, Blackman) and more exotic ones (DPSS, Taylor, ...).
- **Parametric methods**: Yule–Walker, Burg, MA and ARMA, covariance and modified covariance.
- **Non‑parametric methods** based on eigen analysis (e.g., MUSIC) and minimum variance analysis.
- **Multitapering**.

The targeted audience is diverse. Although the use of power spectrum of a signal is fundamental in electrical engineering (e.g. radio communications, radar), it has a wide range of applications from cosmology (e.g., detection of gravitational waves in 2016), to music (pattern detection) or biology (mass spectroscopy).

---

## Quick Installation

Spectrum is available on PyPI:

```bash
pip install spectrum
```

and conda (via conda-forge):

```bash
conda config --append channels conda-forge
conda install spectrum
```

To install the conda executable itself, see `https://www.continuum.io/downloads`.

---

## Links and Contributions

- **Repository**: `https://github.com/cokelaer/spectrum`
- **Contributors**: `https://github.com/cokelaer/spectrum/graphs/contributors`
- **Issues**: `https://github.com/cokelaer/spectrum/issues`
- **Documentation**: `http://pyspectrum.readthedocs.io/`
- **Citation**: Cokelaer et al. (2017), "Spectrum": Spectral Analysis in Python, Journal of Open Source Software, 2(18), 348, doi: `10.21105/joss.00348`

---

## Changelog (summary)

### 0.9.0

- Handle new NumPy API (keeping backward compatibility).
- Include [PR #73](https://github.com/cokelaer/spectrum/pull/73) (thanks to @butala) to speed up FFT.
- Fix rho calculation in Burg algorithm thanks to [PR #82](https://github.com/cokelaer/spectrum/pull/82) from @cl445.
- Remove warnings/deprecations related to pkgresources, NumPy, and SciPy.
- Run Black through entire codebase.

### 0.8.1

- Move CI to GitHub Actions.
- Include Python 3.9 support.
- Include contribution from @tikuma-lshhsc to speed up `eigenfre` module.
- Fix deprecation warnings.

---

## Notebooks (external contributions)

- `http://nbviewer.ipython.org/gist/juhasch/5182528`


