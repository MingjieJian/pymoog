#!/usr/bin/python
# The deprecate functions are stored here. 
# WARNING: the depandence of the functions here are not supported.
from PyAstronomy import pyasl
import os

MOOG_path = '{}/.pymoog/moog_nosm/moog_nosm_NOV2019/'.format(os.environ['HOME'])
MOOG_run_path = '{}/.pymoog/rundir/'.format(os.environ['HOME'])
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
        output_file_name = MOOG_run_path + 'model.mod'
    else:
        output_file_name = save_name
    output_file = open(output_file_name,"w")

    for line in model:
        output_file.write(line + "\n")

    output_file.close()

    model_path = output_file_name
    return model_path