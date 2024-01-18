#!/bin/bash

# Bash program to modify the USER SETUP AREA of Moog.f and Moogsilent.f.

#!/bin/bash

# Bash program to modify the USER SETUP AREA of Moog.f and Moogsilent.f.

# Change the moogpath from '/Users/chris/CODES/moogfeb2017/' to `pwd` or user specified path 
# Note that the code will break the path to multiple lines if it is longer than 60.

path=`pwd`'/'
path='asdfadsfasdfweadfasfcasdfasdfaewfawsefasdfawefagasdfgasdfawefawefawefawfdasdgasgawsgewfaew'

chunk_length=55

if [ ${#path} -gt $chunk_length ]; then
    
    num_lines=$(((${#path} + $chunk_length - 1) / $chunk_length + 1))

    for ((i=1; i<num_lines; i++)); do
        start=$((($i - 1) * $chunk_length))
        replacement_string+="     .  '${path:start:$chunk_length}'"
		replacement_string+=" //"
        replacement_string+=$'\\n'
    done
	replacement_string=${replacement_string::-4}
else 
	replacement_string="     .  '${path}'"
fi

sed -i "22s#.*#$replacement_string#" Moog.f
# sed -i "22s#.*#$replacement_string#" Moogsilent.f

# # Change the machine to user specified type.
# machine='None'
# while [ $machine == 'None' ]
# do
# 	machine=`uname`
# 	if [[ -z "$machine" ]] || [ $machine == 'Linux' ] 
# 	then
# 		machine='pcl'
# 		break
# 	elif [ $machine == 'Darwin' ]
# 	then
# 		machine='mac'
# 		break
# 	else 
# 		machine='pcl'
# 	fi
# done

# # machine='pcl'
# sed -i "29s/'.*'/'$machine'/" Moog.f
# sed -i "29s/'.*'/'$machine'/" Moogsilent.f

# # Install MOOG and MOOGSILENT

# if [ $machine == 'pcl' ]
# then
# 	make -f Makefile.rh64 
# 	make -f Makefile.rh64silent 
# elif [ $machine == 'mac' ]
# then
# 	make -f Makefile.maclap 
# 	make -f Makefile.maclapsilent 
# fi
