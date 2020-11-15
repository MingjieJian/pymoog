import subprocess
import os
import setuptools

if os.environ.get('READTHEDOCS') != 'True':
    # Define MOOGMODELING_path. This path will store the code of moog and any other temporary files in the calculation.
    MOOGMODELING_path = '{}/.pymoog/'.format(os.environ['HOME'])

    # Create the folder according to MOOGMODELING_path
    if not(os.path.isdir(MOOGMODELING_path)):
        os.mkdir(MOOGMODELING_path)

    # Create the folder for calculation
    if not(os.path.isdir(MOOGMODELING_path + 'rundir')):
        os.mkdir(MOOGMODELING_path + 'rundir')
        
    # Copy the files folder into working directory
    if not(os.path.isdir(MOOGMODELING_path + 'files')):
        os.mkdir(MOOGMODELING_path + 'files')
    rm_status = subprocess.run(['rm', '-r', MOOGMODELING_path + 'files/'], stdout=subprocess.PIPE)
    cp_status = subprocess.run(['cp', '-r', 'pymoog/files', MOOGMODELING_path + 'files'], stdout=subprocess.PIPE)

    # untar Kurucz model files
    tar_status = subprocess.run(['tar', '-xzvf', MOOGMODELING_path + 'files/model/kurucz/standard/single.tar.gz', '-C', MOOGMODELING_path + 'files/model/kurucz/standard/'], stdout=subprocess.PIPE)

    # untar the kurucz line list.
    tar_status = subprocess.run(['tar', '-xzvf', MOOGMODELING_path + 'files/linelist/kurucz/kurucz.list.tar.gz', '-C', MOOGMODELING_path + 'files/linelist/kurucz/'], stdout=subprocess.PIPE)

    # Copy the moog_nosm folder to MOOGMODELING_path; if the folder already exist it will be removed first.
    if os.path.isdir(MOOGMODELING_path + '/moog_nosm'):
        rm_status = subprocess.run(['rm', '-r', MOOGMODELING_path + 'moog_nosm'], stdout=subprocess.PIPE)
    mv_status = subprocess.run(['mv', MOOGMODELING_path + 'files/moog_nosm', MOOGMODELING_path + 'moog_nosm'], stdout=subprocess.PIPE)

    # Check the permission of ./install
    chmod_subp = subprocess.run(['chmod', '775', './install.sh'], cwd=MOOGMODELING_path+'moog_nosm/moog_nosm_NOV2019/', stdout=subprocess.PIPE)

    # Check if gfortran is installed in the system.
    # Not done yet

    # Run the MOOGMODELING_path/moog_nosm/moog_nosm_FEB2017/install.sh script
    install = subprocess.run(['bash', './install.sh', 'Linux'], cwd=MOOGMODELING_path+'moog_nosm/moog_nosm_NOV2019/', encoding='UTF-8', stdout=subprocess.PIPE)
    print(install.stdout)

    # Check if MOOG and MOOGSILENT is in the folder
    if not(os.path.isfile(MOOGMODELING_path+'moog_nosm/moog_nosm_NOV2019/MOOG')) or not(os.path.isfile(MOOGMODELING_path+'moog_nosm/moog_nosm_NOV2019/MOOG')):
        raise ValueError("MOOG is not installed correctly!")
    else:
        print('Successfully installed MOOG!')
    
    
with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
      name='pymoog',
      version='0.0.7',
      description='The python wrapper to run LTE spectra synthesis code MOOG.',
      long_description=long_description,
      long_description_content_type="text/markdown",
      url='https://github.com/MingjieJian/pymoog',
      author='Mingjie Jian, Pranav Satheesh and Kruthi Krishna',
      author_email='ssaajianmingjie@gmail.com',
      classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Framework :: IPython",
        "Operating System :: OS Independent",
        "Development Status :: 2 - Pre-Alpha",
        "Topic :: Scientific/Engineering :: Astronomy"
      ],
      python_requires=">=3.5",
      packages=setuptools.find_packages(),
      install_requires=[
          'numpy >= 1.18.0',
          'pandas >= 1.0.0',
          'matplotlib >= 3.1.0',
          'mendeleev >= 0.6.0',
          'scipy >= 1.4.0'
      ],
      include_package_data=True,  
    #   package_data={'': ['moog_nosm/moog_nosm_FEB2017/']},
      zip_safe=False)