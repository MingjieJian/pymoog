import subprocess
import os
import setuptools
import pandas as pd
import numpy as np

if os.environ.get('READTHEDOCS') != 'True':
    # Define MOOGMODELING_path. This path will store the code of moog and any other temporary files in the calculation.
    MOOGMODELING_path = '{}/.pymoog/'.format(os.environ['HOME'])

    # Create the folder according to MOOGMODELING_path
    if not(os.path.isdir(MOOGMODELING_path)):
        os.mkdir(MOOGMODELING_path)
        
    # Copy the files folder into working directory
    if not(os.path.isdir(MOOGMODELING_path + 'files')):
        os.mkdir(MOOGMODELING_path + 'files')
    cp_status = subprocess.run(['cp', '-r', 'pymoog/files', MOOGMODELING_path], stdout=subprocess.PIPE)

    # Download large files from Zenodo 
    zenodo_status = subprocess.run(['zenodo_get', '10.5281/zenodo.7495246', '-o', '~/.pymoog/files/'])
    # Find the latest version of pymoog_lf
    version_list = [i for i in os.listdir('/home/mingjie/.pymoog/files/') if '.tar.gz' in i]
    remove_list = [i for i in version_list if i[-7:] != '.tar.gz']
    version_list = [i for i in version_list if i[-7:] == '.tar.gz']
    version_split_list = [i.split('_')[2].split('.')[0:3] for i in version_list]

    version_list_temp = []
    for version in version_split_list:
        version_list_temp.append([int(version[0][1:]), int(version[1]), int(version[2])])
        
    tar_index = pd.DataFrame(version_list_temp).sort_values([0, 1, 2], ascending=False).index
    remove_list += (list(np.array(version_list)[tar_index[1:]]))
    latest_version = version_list[tar_index[0]]
    # Clear old lf versions
    for version in remove_list:
        rm_status = subprocess.run(['rm', MOOGMODELING_path + '/files/' + version], stdout=subprocess.PIPE)
    # unzip latest version
    tar_status = subprocess.run(['tar', '-xzvf', MOOGMODELING_path + 'files/' + latest_version, '-C', MOOGMODELING_path + 'files/'], stdout=subprocess.PIPE)

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
    if not(os.path.isfile(MOOGMODELING_path+'moog_nosm/moog_nosm_NOV2019/MOOG')) or not(os.path.isfile(MOOGMODELING_path+'moog_nosm/moog_nosm_NOV2019/MOOGSILENT')):
        raise ValueError("MOOG is not installed correctly!")
    else:
        print('Successfully installed MOOG!')
    
with open("README.md", "r") as fh:
    long_description = fh.read()

if os.environ.get('READTHEDOCS') != 'True':
    setuptools.setup(
        name='pymoog',
        version='0.1.0',
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
            'scipy >= 1.4.0',
            'astropy >= 4.0',
            'spectres',
            'tqdm',
            'zenodo_get'
        ],
        include_package_data=True,  
        #   package_data={'': ['moog_nosm/moog_nosm_FEB2017/']},
        zip_safe=False)
else:
        setuptools.setup(
        name='pymoog',
        version='0.1.0',
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
            'scipy >= 1.4.0',
            'spectres',
            'tqdm',
            'zenodo_get'
        ],
        include_package_data=True,  
        #   package_data={'': ['moog_nosm/moog_nosm_FEB2017/']},
        zip_safe=False)