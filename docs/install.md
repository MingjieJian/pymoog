# Installation

## Install gfortran compilar

`pymoog` use [gfortran](https://gcc.gnu.org/wiki/GFortran), a free Fortran compiler as a part of GCC.
Although it is possible to use other Fortran compiler to compile MOOG, using gfortran is (to my experience) the mose easy one and thus pymoog requires gcc to be installed before installation. 

- Check whether gfortran is installed in your computer

Type `gfortran` in the terminal; if you get:

```
gfortran: fatal error: no input files
compilation terminated.
```

then gfortran is already installed and please go to [Install pymoog](#install-pymoog). 

If you get:

```
command not found: gfortran
```

then you need to install gfortran.

- Linux: `sudo apt install gcc`
-  Windows: please use [Windows Subsystem for Linux](https://docs.microsoft.com/en-us/windows/wsl/) and follow the instruction for Linux.
- Mac: please refer to [this post](https://discussions.apple.com/thread/8336714); the main point is using `brew` to install (`brew install gcc`).
I have no experience on using Mac so finger crossed. 

## <a name="install-pymoog"></a>Install pymog

- Using pypi (recommended)
    - `pip install pymoog`
- From github
    - `clone` this repository and `cd` into the corresponding folder;
    - `pip install .`