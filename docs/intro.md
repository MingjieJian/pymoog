# Introduction to pymoog

`pymoog` is a python3 wrapper for running the LTE spectrum synthesis part of the code [MOOG](https://www.as.utexas.edu/~chris/moog.html) written by Chris Sneden.
It wraps up the (a bit) teidous steps for generating a synthetic spectra into four python commands, while retaining the functions provided by MOOG.

Here you can:
- generate a synthetic spectra
- alter the stellar parameters for the spectra, such as $T_\mathrm{eff}$, metallicity, abundance ratios or resolution.
- determine some stellar parameters, such as microturbulance velocity or abundance ratios.

without considering:
- where to grab stellar atmosphere models and line lists
- how to interpolate the models.


```{warning}
We cannot use any code blindly, so it is important to read the output files of MOOG, i.e., MOOG.out1/2/3, and the [user manual](https://www.as.utexas.edu/~chris/codes/WRITEnov2019.pdf) of MOOG when you find something may gone wrong.
```
