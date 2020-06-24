#!/bin/bash

# Bash program to modify the USER SETUP AREA of Moog.f and Moogsilent.f.

# Change the moogpath from '/Users/chris/CODES/moogfeb2017/' to `pwd` or user specified path 

path=`pwd`'/'
sed "22s?'.*'?'$path'?" Moog_bak.f > Moog.f
sed "22s?'.*'?'$path'?" Moogsilent_bak.f > Moogsilent.f

# Change the machine to user specified type.
machine='None'
while [ $machine == 'None' ]
do
	echo 'Please specify the machine type: either Linux, or Mac. (Default Linux)'
	read machine
	if [[ -z "$machine" ]] || [ $machine == 'Linux' ] 
	then
		machine='pcl'
		break
	elif [ $machine == 'Mac' ]
	then
		machine='mac'
		break
	else 
		echo 'Please choose between Linux and Mac, or you need to modify the files and install MOOG by yourself.'
		machine='None'
	fi
done

sed -i "29s/'.*'/'$machine'/" Moog.f
sed -i "29s/'.*'/'$machine'/" Moogsilent.f

# Install MOOG and MOOGSILENT

if [ $machine == 'pcl' ]
then
	make -f Makefile.rh64 
	make -f Makefile.rh64silent 
elif [ $machine == 'mac' ]
then
	make -f Makefile.maclap 
	make -f Makefile.maclapsilent 
fi
