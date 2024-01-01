# Installation

## Install gfortran compilar

`pymoog` use [gfortran](https://gcc.gnu.org/wiki/GFortran), a free Fortran compiler as a part of GCC.
Although it is possible to use other Fortran compiler to compile MOOG, using gfortran is (to my experience) the mose easy one and thus pymoog requires gcc to be installed before installation. 

First, we check whether gfortran is installed in your computer.

Type `gfortran` in the terminal; if you get:

```bash
gfortran: fatal error: no input files
compilation terminated.
```

then gfortran is already installed and please go to next section. 

If you get:

```bash
command not found: gfortran
```

then you need to install gfortran:

- Linux: `sudo apt install gcc`
-  Windows: please use [Windows Subsystem for Linux](https://docs.microsoft.com/en-us/windows/wsl/) and follow the instruction for Linux.
- Mac: please refer to [this post](https://discussions.apple.com/thread/8336714); the main point is using `brew` to install (`brew install gcc`).
I have no experience on using Mac so, finger crossed. 

## Install pymoog

- Using pypi (recommended)
    - `pip install pymoog`
- From github
    - `clone` this repository and `cd` into the corresponding folder;
    - `pip install .`

Note that `pymoog` requires some large files (atmosphere models and line lists, stored in [here](https://zenodo.org/record/7495246#.Y9Vh_cmSljE)) to run , and they will be downloaded during `pip install`.
Please test your conntection to Zenodo by accessing the link first.
Thus the installation may takes 10 or 20 minutes without anything prompting in the command line.
When upgrading `pymoog`, the program check if there is any newer version of the large files.
If the local files are up-to-date, then the download will be skipped to speed up the installation process. 

## Uninstall pymoog

- Uninstall in pip: `pip uninstall pymoog`.
- Remove all the files in pymoog running folder: `rm -r ~/.pymoog`