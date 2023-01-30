# pymoog

[![DOI](https://zenodo.org/badge/10543633.svg)](https://zenodo.org/badge/latestdoi/10543633) [![Documentation Status](https://readthedocs.org/projects/pymoog/badge/?version=latest)](https://pymoog.readthedocs.io/en/latest/?badge=latest)

`pymoog` is a python3 wrapper for running the LTE spectrum synthesis part of the code [MOOG](https://www.as.utexas.edu/~chris/moog.html) written by Chris Sneden.
It wraps up the (a bit) teidous steps for generating a synthetic spectra into four python commands, while retaining the functions provided by MOOG.
Besides, it also provides some other functions for analysing the MOOG result, mainly contribution function and fitting stellar parameters.

Here you can:
- generate a synthetic spectra
- alter the stellar parameters for the spectra, such as $T_\mathrm{eff}$, metallicity, abundance ratios or resolution.
- determine some stellar parameters, such as microturbulance velocity or abundance ratios.

without considering:
- where to grab stellar atmosphere models and line lists
- how to interpolate the models.

Documentation (including the installation) is available in [here](https://pymoog.readthedocs.io/en/latest/).

## Acknowledgement

This package has made use of the VALD database, operated at Uppsala University, the Institute of Astronomy RAS in Moscow, and the University of Vienna.