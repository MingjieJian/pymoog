#!/usr/bin/python
import numpy as np
import pandas as pd
import mendeleev as md
import re
import sys
import subprocess
import pkg_resources
import os

MOOG_file_path = '{}/.pymoog/files/'.format(os.environ['HOME'])

## Convert the element column to element specics

    
def save_linelist(linelist_all, sub_ll_name, wav_start=None, wav_end=None, header=None):
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
    '''
    
    # Crop the line list according to wavelength, if needed.
    index = linelist_all['wavelength'] > 0
    if wav_start != None:
        index = index & (linelist_all['wavelength'] > wav_start)
    if wav_end != None:
        index = index & (linelist_all['wavelength'] < wav_end) 
           
    sub_linelist = linelist_all[index]
    sub_linelist.reset_index(drop=True, inplace=True)
    if np.any(abs(sub_linelist['C6'].values) > 1e-25):
        output_format = '{:10.4f}{:10.5f}{:10.4f}{:10.3f}{:10.3f}{:10.3f}{:10.3f}\n'
    elif np.any(abs(sub_linelist['C6'].values) < 1e-25):
        output_format = '{:10.4f}{:10.5f}{:10.4f}{:10.3f}{:10.2E}{:10.3f}{:10.3f}\n'
    # Remove the last column if no EW values.
    if len(sub_linelist.columns) == 6:
        output_format = output_format[:-9] + '\n'
    with open(sub_ll_name, 'w') as file:
        if header == None:
            file.write('Linelist\n')
        else:
            file.write('{}\n'.format(header))
        for i in range(len(sub_linelist)):
            line = output_format.format(*sub_linelist.iloc[i].values).replace('nan', '   ')
            file.write(line)

def read_linelist(linelist_name, loggf_cut=None):
    '''
    Read the post-processed linelist.
    
    Parameters
    ----------
    linelist_name : str
        The MOOG format line list
    loggf_cut : float, optional
        Cut on loggf (only save for the lines with loggf > loggf_cut)
    '''
    
    available_line_list = ['ges', 'ges_hfs_iso', 'ges_nohfs_noiso', 'vald_3000_11000', 'vald_11000_24000', 'mb99_j', 'mb99_k', 'apogee', 'kurucz']
    
    if linelist_name[-5:] != '.list' and linelist_name in available_line_list:
        # Read built in line list
        if linelist_name == 'ges':
                    linelist_name = 'ges_hfs_iso'
        linelist_name = MOOG_file_path + 'linelist/{}/{}.list'.format(linelist_name.split('_')[0], linelist_name)
    elif linelist_name[-5:] == '.list':
        pass
    else:
        raise ValueError("Built in line list type not recognized. Please use one of the following:\n              'ges', 'ges_hfs_iso', 'ges_nohfs_noiso', 'vald_3000_11000', 'vald_11000_24000', 'mb99_j', 'mb99_k', 'kurucz' or 'apogee'.")
    linelist = pd.read_fwf(linelist_name,
            colspecs=[(0,11), (11,21), (21,31), (31,41), (41,51), (51,61), (61,71)],
            names=['wavelength', 'id', 'EP', 'loggf', 'C6', 'D0', 'EW'],
            skiprows=1)
    # MOOG seems to crash if there is line with EP larger than 50eV, so they are removed.
    # Need to be test for other line lists
    linelist = linelist[(linelist['EP'] <= 50)]
    if loggf_cut != None:
        linelist = linelist[(linelist['loggf'] >= loggf_cut)]
        linelist.reset_index(drop=True, inplace=True)
    return linelist
