#!/usr/bin/python
from re import sub
import numpy as np
import pandas as pd
import os
import platform
from . import private

MOOG_path = '{}/.pymoog/moog_nosm/moog_nosm_NOV2019/'.format(os.environ['HOME'])
MOOG_run_path = '{}/.pymoog/rundir/'.format(os.environ['HOME'])
MOOG_file_path = '{}/.pymoog/files/'.format(os.environ['HOME'])

## Convert the element column to element specics

def save_linelist(linelist_all, sub_ll_name, wav_start=None, wav_end=None, header=None, negative=False):
    '''
    Save the linelist in MOOG format into specified position.
    
    Parameters
    ----------
    linelist_all : pandas.Dataframe
        The Dataframe of linelist in MOOG format
    sub_ll_name : str
        The name of the line list to be saved into.
    wav_start : float
        Start wavelength of the line list.
    end_start : float
        End wavelength of the line list.
    type : str, = 'vald'
        Type of the line list. Now only 'vald' is supported.
    negative : bool
        Switch to permit negative wavelength. 
    '''
    
    # Crop the line list according to wavelength, if needed.
    if not(negative):
        index = linelist_all['wavelength'] > 0
    else:
        index = np.abs(linelist_all['wavelength']) >= 0
    if wav_start != None:
        index = index & (linelist_all['wavelength'] > wav_start)
    if wav_end != None:
        index = index & (linelist_all['wavelength'] < wav_end) 
           
    sub_linelist = linelist_all[index]
    sub_linelist.reset_index(drop=True, inplace=True)
    
    # Judge if the length of the line list is 0; if so raise an error.
    if len(sub_linelist) == 0:
        raise ValueError('The length of line list is 0. Consider enalrge the wavelength or check the input line list.')
    
    # Decidcde which format to save the linelist according to C6 value.
    if np.any(abs(sub_linelist['C6'].values) > 1e-25):
        output_format = '%10.3f%10.5f%10.4f%10.3f%10.3f%10.3f%10.3f'
    elif np.any(abs(sub_linelist['C6'].values) < 1e-25):
        output_format = '%10.3f%10.5f%10.4f%10.3f%10.2E%10.3f%10.3f'
    
    # Remove the last column if no EW values.
    if len(sub_linelist.columns) == 6:
        output_format = output_format[:-6]
    np.savetxt(sub_ll_name, np.array(sub_linelist), fmt=output_format)
    if 'linux' in platform.system().lower():
        run_status = private.subprocess.run(['sed', '-i', 's/nan/   /g', sub_ll_name], capture_output=True)
    elif 'darwin' in platform.system().lower():
        run_status = private.subprocess.run(['sed', '-i', "''", 's/nan/   /g', sub_ll_name], capture_output=True)
    else:
        # Same as Linux
        run_status = private.subprocess.run(['sed', '-i', 's/nan/   /g', sub_ll_name], capture_output=True)
    if run_status.returncode != 0:
        raise ValueError('NaN may not be removed correctly in the line list. The stderr text is: {}'.format(run_status.stderr)) 
    if header == None:
        header = 'Linelist'
    if 'linux' in platform.system().lower():
        run_status = private.subprocess.run(['sed', '-i', '1 i\{}'.format(header), sub_ll_name], capture_output=True)
    elif 'darwin' in platform.system().lower():
        run_status = private.subprocess.run(['sed', '-i', "''", '1 i\{}'.format(header), sub_ll_name], capture_output=True)
    else:
        # Same as linux
        run_status = private.subprocess.run(['sed', '-i', '1 i\{}'.format(header), sub_ll_name], capture_output=True)

def read_linelist(linelist_name, loggf_cut=None, mode='default'):
    '''
    Read the post-processed linelist.
    
    Parameters
    ----------
    linelist_name : str
        The MOOG format line list
    loggf_cut : float, optional
        Cut on loggf (only save for the lines with loggf > loggf_cut)
    mode : str, default 'default'
        Reading mode for reading line-list. 'default' will first try to read using 'npy' mode then 'ascii' mode if the corresponding .npy file does not exist. Note that the efficiency of 'npy' mode is much higher than 'ascii' mode.
    '''
    
    available_line_list = ['ges', 'ges_hfs_iso', 'ges_nohfs_noiso', 'vald_3000_24000', 'vald_winered', 'mb99_j', 'mb99_k', 'apogee', 'kurucz', 'kurucz_winered']

    if (linelist_name[-5:] != '.list' and linelist_name[-4:] != '.npy') and linelist_name in available_line_list:
        # Read built-in line list
        if linelist_name == 'ges':
            linelist_name = 'ges_hfs_iso'
        if mode == 'default':
            linelist_name_full = MOOG_file_path + '/pymoog_lf/linelist/{}/{}.npy'.format(linelist_name.split('_')[0], linelist_name)
            mode = 'npy'
            if not(os.path.exists(linelist_name_full)):
                linelist_name_full = MOOG_file_path + '/pymoog_lf/linelist/{}/{}.list'.format(linelist_name.split('_')[0], linelist_name)
                mode = 'ascii'
                if not(os.path.exists(linelist_name_full)):
                    raise ValueError('Neither npy nor ascii format of internal line list exists.')
        elif mode == 'npy':
            linelist_name_full = MOOG_file_path + '/pymoog_lf/linelist/{}/{}.npy'.format(linelist_name.split('_')[0], linelist_name)
        elif mode == 'ascii':
            linelist_name_full = MOOG_file_path + '/pymoog_lf/linelist/{}/{}.list'.format(linelist_name.split('_')[0], linelist_name)
        else:
            raise ValueError('mode must be "default", "npy" or "ascii".')
    elif linelist_name[-5:] == '.list':
        linelist_name_full = linelist_name
        mode = 'ascii'
    elif linelist_name[-4:] == '.npy':
        linelist_name_full = linelist_name
        mode = 'npy'
    else:
        raise ValueError("Built in line list type not recognized. Please use one of the following:\n              'ges', 'ges_hfs_iso', 'ges_nohfs_noiso', 'vald_3000_24000', 'vald_winered', 'mb99_j', 'mb99_k', 'kurucz', 'kurucz_winered' or 'apogee'.")
    
    if mode == 'npy':
        linelist_array = np.load(linelist_name_full, allow_pickle=True)
        linelist = pd.DataFrame(linelist_array, columns=['wavelength', 'id', 'EP', 'loggf', 'C6', 'D0', 'EW'])
    elif mode == 'ascii':
        linelist = pd.read_fwf(linelist_name_full, colspecs=[(0,11), (11,21), (21,31), (31,41), (41,51), (51,61), (61,71)], names=['wavelength', 'id', 'EP', 'loggf', 'C6', 'D0', 'EW'], skiprows=1)
        
    # MOOG seems to crash if there is line with EP larger than 50eV, so they are removed.
    # Need to be test for other line lists
    linelist = linelist[(linelist['EP'] <= 50)]
    if loggf_cut != None:
        linelist = linelist[(linelist['loggf'] >= loggf_cut)]
        linelist.reset_index(drop=True, inplace=True)
    return linelist