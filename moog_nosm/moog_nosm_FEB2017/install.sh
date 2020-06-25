#!/bin/bash

# Bash program to modify the USER SETUP AREA of Moog.f and Moogsilent.f.

# Change the moogpath from '/Users/chris/CODES/moogfeb2017/' to `pwd` or user specified path 

path=`pwd`'/'
sed "22s?'.*'?'$path'?" Moog_bak.f > Moog.f
sed "22s?'.*'?'$path'?" Moogsilent_bak.f > Moogsilent.f

sed -i "_old" "29s/'.*'/'$machine'/" Moog.f
sed -i "_old" "29s/'.*'/'$machine'/" Moogsilent.f

# Install MOOG and MOOGSILENT
make -f Makefile.rh64 
make -f Makefile.rh64silent