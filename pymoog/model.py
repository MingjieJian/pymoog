#!/usr/bin/python
# This is the python model dealing with the stellar atmosphere model.
import os
from PyAstronomy import pyasl
import subprocess
import numpy as np
import re
import pandas as pd
from pymoog import line_data

MOOG_path = '{}/.pymoog/moog_nosm/moog_nosm_FEB2017/'.format(os.environ['HOME'])
MOOGrun_path = '{}/.pymoog/rundir/'.format(os.environ['HOME'])
MOOG_file_path = '{}/.pymoog/files/'.format(os.environ['HOME'])
      
def KURUCZ_download(teff, logg, m_h, save_name=None):
    '''
    Download the Kurucz ATLAS9 model using pyasl.

    Parameters
    ----------
    teff : int
        The effective temperature of the model
    logg : float
        logg value of the model
    m_h : float
        [M/H] value (overall metallicity) of the model
    save_name : str
        The path to save the model file (including the name).

    Returns
    ----------
    self.model_path : str
        The path and name of saved model.

    '''

    # The model is invoked here
    model = pyasl.getKuruczModel(teff, logg, m_h)

    # Save the model file to specified position
    if save_name == None:
        output_file_name = MOOGrun_path + 'model.mod'
    else:
        output_file_name = save_name
    output_file = open(output_file_name,"w")

    for line in model:
        output_file.write(line + "\n")

    output_file.close()

    model_path = output_file_name
    return model_path


def KURUCZ_convert(model_path, vmicro=2.0, abun_change=None, converted_model_path=None):
    '''
    Convert the model file from Kurucz format in to MOOG format.

    Parameters
    ----------
    model_path : str, optional
        The path of donloaded model file.
    v_micro : float, default 2.0
        Microturbulance velocity of the spectra.
    abun_change : list of pairs [int, float]
        Abundance change, have to be a list of pairs of atomic number and [X/Fe] values.
    converted_model_path : str, optional
        The name of converted model. Default will be saved into MOOG working folder.
    '''

    model_file = open(model_path)
    # Convert the model files into MOOG format.

    # Read and save the first two lines (except 'TITLE ') into header.
    header = ['Kurucz model: ' + model_file.readline()]
    header = header + [model_file.readline()]
    try:
        m_h = re.findall(r'\[(.*)\]', header[1])[0]
    except IndexError:
        m_h = 0
        print(header)

    # Read the abundance change as well as model lines.
    temp = model_file.readline() + model_file.readline() + model_file.readline()

    abun_list = ''
    temp = model_file.readline()
    while temp[:17] == ' ABUNDANCE CHANGE':
        abun_list = abun_list + temp[17:]
        temp = model_file.readline()
    abun = np.array(abun_list.split(), dtype='f').reshape(int(len(abun_list.split())/2), 2)

    # Load the abundances from Asplund 2009 (which MOOG used; hard-coded 20180531)
    xabu = [12.00,10.93, 1.05, 1.38, 2.70, 8.43, 7.83, 8.69, 4.56, 7.93,
            6.24, 7.60, 6.45, 7.51, 5.41, 7.12, 5.50, 6.40, 5.03, 6.34,
            3.15, 4.95, 3.93, 5.64, 5.43, 7.50, 4.99, 6.22, 4.19, 4.56,
            3.04, 3.65, 2.30, 3.34, 2.54, 3.25, 2.52, 2.87, 2.21, 2.58,
            1.46, 1.88,-5.00, 1.75, 0.91, 1.57, 0.94, 1.71, 0.80, 2.04,
            1.01, 2.18, 1.55, 2.24, 1.08, 2.18, 1.10, 1.58, 0.72, 1.42,
            -5.00, 0.96, 0.52, 1.07, 0.30, 1.10, 0.48, 0.92, 0.10, 0.84,
            0.10, 0.85,-0.12, 0.85, 0.26, 1.40, 1.38, 1.62, 0.92, 1.17,
            0.90, 1.75, 0.65,-5.00,-5.00,-5.00,-5.00,-5.00,-5.00, 0.02,
            -5.00,-0.54,-5.00,-5.00,-5.00]

    # Read the model lines
    temp = temp.split()
    model_lines = []
    model_linen = int(temp[2])
    for i in range(int(temp[2])):
        model_lines.append(model_file.readline().split()[:7])

    # Prepare the microtrubulance value.
    vmicro = '{}E00'.format(vmicro)

    # Write the model file.
    # header, abun09, model_lines and vmicro
    if converted_model_path == None:
        c_model_path = MOOGrun_path + '/model.mod'
    else:
        c_model_path = converted_model_path
    c_model_file = open(c_model_path, 'w')

    # Header part
    c_model_file.writelines(header[0])
    c_model_file.writelines(header[1])

    # Model part
    c_model_file.writelines('ntau=       ' + str(model_linen) + '\n')
    for i in model_lines:
        c_model_file.writelines(' ' + ' '.join(i) + '\n')

    # Microturbulant velocity part
    c_model_file.writelines('    ' + vmicro + '\n')

    # Element shift part
    if abun_change != None:
        abun_change_num = len(abun_change)
        c_model_file.writelines('NATOMS      {}   {}\n'.format(abun_change_num, m_h))
        for abun in abun_change:
            c_model_file.writelines('      {:.2f}    {:.2f}\n'.format(abun[0], xabu[abun[0]-1]+abun[1]+float(m_h)))
    else:
        c_model_file.writelines('NATOMS      0   {}\n'.format(m_h))

    # Molecular line switches (temporary closed)
    c_model_file.writelines('NMOL        0')
    c_model_file.close()

    # cv_situation = os.path.isfile(c_model_path)
    # return c_model_file, cv_situation  
    
def read_model(model_path):
    '''
    Read the model into pandas.DataFrame for interpolation

    Parameters
    ----------
    model_path : str, optional
        The path of downloaded model file.
    '''
    
    model_df = pd.read_csv(model_path)
    
def value2pm(value):
    '''
    Transform the metallicity value to Kurucz format.
    Example: -1.0 -> m10
    
    Parameters
    ----------
    value : float
        The value of metallicity.
    '''
    if value < 0:
        return 'm{:02.0f}'.format(np.abs(value)*10)
    else:
        return 'p{:02.0f}'.format(np.abs(value)*10)
    
def split_kurucz_model():
    '''
    Split the Kurucz model into single. Internal function.
    ''' 
    grid_kurucz = pd.read_csv('files/grid_points_kurucz.csv')
    for m_h in grid_kurucz.groupby('m_h').size().index:
        file = open('files/model/kurucz/standard/a{}k2.dat'.format(value2pm(m_h)))
        content = file.readlines()

        is_first = True
        for line in content:
            if 'TEFF' in line:
                if not(is_first):
                    with open('files/model/kurucz/standard/single/teff{:.0f}logg{:.1f}m_h{:+.1f}.dat'.format(teff, logg, m_h), 'w') as file:
                        file.writelines(single_model)
                teff = float(line[5:13])
                logg = float(line[21:30])
                is_first = False
                single_model = [line]
            else:
                single_model.append(line)