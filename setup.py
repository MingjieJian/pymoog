import subprocess
import os

# Define MOOGMODELING_path. This path will store the code of moog and any other temporary files in the calculation.
MOOGMODELING_path = '{}/.pymoog/'.format(os.environ['HOME'])

# Create the folder according to MOOGMODELING_path
if not(os.path.isdir(MOOGMODELING_path)):
    os.mkdir(MOOGMODELING_path)

# Create the folder for calculation
if not(os.path.isdir(MOOGMODELING_path + 'rundir')):
    os.mkdir(MOOGMODELING_path + 'rundir')

# Copy the moog_nosm folder to MOOGMODELING_path; if the folder already exist it will be removed first.
if os.path.isdir(MOOGMODELING_path + 'moog_nosm'):
    rm_status = subprocess.run(['rm', '-r', MOOGMODELING_path + 'moog_nosm'], stdout=subprocess.PIPE)
cp_status = subprocess.run(['cp', '-r', 'moog_nosm', MOOGMODELING_path + 'moog_nosm'], stdout=subprocess.PIPE)

# Check the permission of ./install
chmod_subp = subprocess.run(['chmod', '775', './install.sh'], cwd=MOOGMODELING_path+'moog_nosm/moog_nosm_FEB2017/', stdout=subprocess.PIPE)

# Check if gfortran is installed in the system.
# Not done yet

# Run the MOOGMODELING_path/moog_nosm/moog_nosm_FEB2017/install.sh script
install = subprocess.run(['bash', './install.sh', 'Linux'], cwd=MOOGMODELING_path+'moog_nosm/moog_nosm_FEB2017/', encoding='UTF-8', stdout=subprocess.PIPE)
print(install.stdout)

# Check if MOOG and MOOGSILENT is in the folder
if not(os.path.isfile(MOOGMODELING_path+'moog_nosm/moog_nosm_FEB2017/MOOG')) or not(os.path.isfile(MOOGMODELING_path+'moog_nosm/moog_nosm_FEB2017/MOOG')):
    raise ValueError("MOOG is not installed correctly!")
else:
    print('Successfully installed MOOG!')