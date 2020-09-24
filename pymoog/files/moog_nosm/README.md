# moog_nosm
MOOG without the dependent of SM.

[MOOG](http://www.as.utexas.edu/~chris/moog.html) is a software of generating synthestic spectra and fitting it with observed spectra to get abundances from stellar model and line list. It is written by Chris Sneden.

Original MOOG is related to plotting software [SM](http://www.astro.princeton.edu/~rhl/sm/), which is (maybe) hard to install and call by MOOG. The MOOG in this repository is undependent on SM, so the command MOOG and MOOGSILENT is basically the same.

Beside this some lines are also modified:

- Pi value inside Smooth.f is changed to 3.1415927.
- Some judgement is added to make makefile of MOOGSILENT run normally.
- Input and output format of pixel value of synthetic spectrum inside Binplotprep.f, Smooth.f, Synpop.f, Synspec.f and Vargauss.f is change to f8.4, adding a space between numbers when they are small than 0.
- The output format of Lineinfo.f are modified to output the opacities when `atmosphere` and `lines` are set to 2 and 4.

The modified MOOG of versions NOV2019 and FEB2017 are provided, and their installation methods are identical.

## Install

Here the install method of version NOV2019 are present as example. For the installation of version NOV2019, change `FEB2017` to `NOV2019`.

1. Clone or download the respository.
2. Enter the MOOG folder by `cd PATH_TO_REPO/moog_nosm_NOV2019`.
3. Run the bash file `./install.sh` or `bash ./install.sh` 
4. There will be a prompt asking which kind of machine the installation is on; choose either `Linux` (or simply hit Enter) or `Mac`. I test it in Ubuntu but not for Mac so please be cautious.
5. Then installation should go by itself.
6. If you see some prompt of `make sure that you have entered the proper parameters for MOOG into the FORTRAN source driver routine Moog.f !!!!!!!!!!!!` then you should see a `MOOG` and `MOOGSILENT` file present in the current folder.
7. Done.

## Warnings

Some warnings appear during the make command.
Most of them are no-harm, such as:
- `Warning: Line truncated at (1) [-Wline-truncation]` and
- `Warning: Possible change of value in conversion from REAL(8) to INTEGER(4) at (1) [-Wconversion]`.

They are intended and tried to get rid of some of those.

However there are two may worth a slight change:
- `Warning: ‘ikount’ may be used uninitialized in this function [-Wmaybe-uninitialized]` in `Molquery.f`
- `Warning: Rank mismatch in argument ‘xnew’ at (1) (rank-1 and scalar) [-Wargument-mismatch]` in `OpacHydrogen.f`

I leave them as in the initial code.
