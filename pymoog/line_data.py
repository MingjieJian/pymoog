#!/usr/bin/python
import numpy as np
import pandas as pd
import mendeleev as md
import re
import sys
import subprocess
import pkg_resources
import os
from . import synth
from . import weedout

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
    if np.any(abs(sub_linelist['C6'].values) > 1e-25):
        output_format = '{:10.3f}{:10.5f}{:10.4f}{:10.3f}{:10.3f}{:10.3f}{:10.3f}\n'
    elif np.any(abs(sub_linelist['C6'].values) < 1e-25):
<<<<<<< HEAD
        output_format = '{:10.3f}{:10.5f}{:10.4f}{:10.3f}{:10.2E}{:10.3f}{:10.3f}\n'
=======
        output_format = '{:10.4f}{:10.5f}{:10.4f}{:10.3f}{:10.2E}{:10.3f}{:10.3f}\n'
>>>>>>> master
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

def find_lines(linelist_keep, linelist_all):
    line_index_keep = []
    for i in range(len(linelist_keep)):
        indice = (np.abs(linelist_all['wavelength'] - linelist_keep.loc[i, 'wavelength']) < 0.001)
        for col in ['id', 'EP', 'loggf']:
            indice = indice & (np.abs(linelist_all[col] - linelist_keep.loc[i, col]) < 0.001)
        if len(linelist_all[indice]) == 0:
            raise ValueError('No match line found.')
        line_index_keep.append(linelist_all[indice].index.values[0])
    return line_index_keep


def find_single_dominant_line(line_wav_input, teff, logg, fe_h, resolution, r_d_blend_thre=0.1, line_list='ges', weedout_switch=False, search_half_width=0.5, linelist_serach=False, abun_change=None):

    # Establish the linelist
    linelist_all = read_linelist(line_list)
    save_linelist(linelist_all, MOOG_run_path + 'all.list', wav_start=line_wav_input-search_half_width, wav_end=line_wav_input+search_half_width)
    linelist_all = read_linelist(MOOG_run_path + 'all.list')

    # Calculate the blending ratio
    s = synth.synth(teff, logg, fe_h, line_wav_input-search_half_width-1, line_wav_input+search_half_width+1, resolution, line_list=MOOG_run_path + 'all.list')

    # Whole spectra 
    if abun_change is not None:
        s.prepare_file(abun_change=abun_change)
    else:
        s.prepare_file()
    s.run_moog()
    s.read_spectra()
    wav_all, flux_all = s.wav, s.flux

    # weedout lines
    if weedout_switch != False:
        w = weedout.weedout(teff, logg, fe_h, line_wav_input-search_half_width, line_wav_input+search_half_width, line_list=line_list)
        w.prepare_file()
        w.run_moog()
        
    # Target line exclude
    if weedout_switch:
        linelist_keep = read_linelist(MOOG_run_path + 'keep.list')
    else:
        linelist_keep = read_linelist(MOOG_run_path + 'all.list')
    
    line_index_keep = find_lines(linelist_keep, linelist_all)

    r_blend_depth_list = []
    for line_index in line_index_keep:
        s = synth.synth(teff, logg, fe_h, line_wav_input-search_half_width-1, line_wav_input+search_half_width+1, resolution, line_list=MOOG_run_path+'line.list')
        if abun_change is not None:
            s.prepare_file(abun_change=abun_change)
        else:
            s.prepare_file()
        linelist_exclude = linelist_all.drop(line_index).reset_index(drop=True)
        save_linelist(linelist_exclude, MOOG_run_path + 'line.list')
        s.run_moog()
        s.read_spectra()
        wav_exclude, flux_exclude = s.wav, s.flux

        # Target line only
        linelist_target = linelist_all.loc[line_index:line_index].reset_index(drop=True)
        line_wavlength = linelist_target.loc[0, 'wavelength']
        line_loggf = linelist_target.loc[0, 'loggf']
        line_EP = linelist_target.loc[0, 'EP']
        if abun_change is not None:
            s.prepare_file(abun_change=abun_change)
        else:
            s.prepare_file()
        save_linelist(linelist_target, MOOG_run_path + 'line.list')
        s.run_moog()
        s.read_spectra()
        wav_target, flux_target = s.wav, s.flux
        # plt.plot(wav_exclude, flux_exclude)
        # print(linelist_target)

        # Calculate the EW and blending fraction
        EW = (np.sum(1-flux_all)*0.02 - np.sum(1-flux_exclude)*0.02) * 1000
        depth = 1 - np.min(flux_all[np.abs(wav_all-line_wavlength) <= 0.03])
        r_blend_depth = (1-flux_exclude[np.argmin(np.abs(wav_exclude-line_wavlength))]) / (1-flux_all[np.argmin(np.abs(wav_all-line_wavlength))])

        r_blend_depth_list.append(r_blend_depth)

    linelist_keep['r_blend_depth'] = r_blend_depth_list

    if len(line_index_keep) > 0:
        try:
            target_line_index = np.abs(linelist_keep.loc[linelist_keep['r_blend_depth'] < 0.1, 'wavelength'] - line_wav_input).sort_values().index[0]
            target_line = linelist_keep.loc[target_line_index:target_line_index].reset_index(drop=True)
        except IndexError:
            # No dominant line is found
            target_line = pd.DataFrame(np.array([np.nan]*8)).T
            target_line.columns = ['wavelength', 'id', 'EP', 'loggf', 'C6', 'D0', 'EW', 'r_blend_depth']
    else:
        # No line is found
        target_line = pd.DataFrame(np.array([np.nan]*8)).T
        target_line.columns = ['wavelength', 'id', 'EP', 'loggf', 'C6', 'D0', 'EW', 'r_blend_depth']

    if linelist_serach:
        return target_line, linelist_keep
    else:
        return target_line