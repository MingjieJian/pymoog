#!/usr/bin/python
from re import sub
import numpy as np
import pandas as pd
import os
import platform
from . import private
from . import synth
from . import weedout
from . import rundir_num

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

def find_lines(linelist_keep, linelist_all, max_del_wav=0.05):
    line_index_keep = []
    for i in linelist_keep.index:
        indice = (np.abs(linelist_all['wavelength'] - linelist_keep.loc[i, 'wavelength']) < max_del_wav)
        for col in ['id', 'EP']:
            # Note: some difference in loggf may appear in different version of the line list; so it is better to keep the version same, and here loggf is not used as the criteria for distinguishing lines.
            indice = indice & (np.abs(linelist_all[col] - linelist_keep.loc[i, col]) < 0.001)
        if len(linelist_all[indice]) == 0:
            raise ValueError('No match line found.')
        else:
            line_index_keep.append(linelist_all[indice].index.values[0])
    return line_index_keep

def find_single_dominant_line(line_wav_input, teff, logg, m_h, resolution, line_list='ges', weedout_switch=False, search_half_width=0.5, include_strong_linelist=False, r_blen_thres=0.1, abun_change=None):

    '''
    Find the dominant line from a line list.

    Parameters
    ----------
    line_wav_input : float
        Central wavelength for the searching.
    teff : float
        The effective temperature of the model
    logg : float
        logg value of the model
    m_h : float
        [M/H] value (overall metallicity) of the model
    resolution : float
        Resolution of the synthetic spectra; this will passed to MOOG and convolute with initial spectra.
    line_list : str or pd.DataFrame, default vald_3000_24000 
        The name of the linelist file.
    weedout_switch : bool or float, default False
        The switch for running weedout driver before synth. If False then weedout is not run; if True the weedout is run with kappa_ratio=0.01, and if a float (> 0 and < 1) is given then weedout is run with the kappa_ratio set as the number.
    search_half_width : float, default 0.5
        The +- width for searching the dominant line.
    include_strong_linelist : bool, default False
        Whether include all the linelist after weedout as a separate output.
    r_blen_thres : float, default 0.1
        The threshold of blending ratio. Only the line with blending ratio smaller than r_blen_thres can be selected as dominant line.
    abun_change : dict of pairs {int:float, ...}
            Abundance change, have to be a dict of pairs of atomic number and [X/Fe] values.

    Returns
    ----------
    dominant_line : pandas.DataFrame 
        The dataframe containing the dominant line.
    linelist_keep : pandas.DataFrame, optional
        The line list after weedout. Only appear when include_strong_linelist is True.
    '''

    # Establish the linelist
    linelist_all = read_linelist(line_list)
    linelist_all = linelist_all[np.abs(linelist_all['wavelength']-line_wav_input) < search_half_width]

    # Calculate the blending ratio
    s = synth.synth(teff, logg, m_h, line_wav_input-search_half_width-1, line_wav_input+search_half_width+1, resolution, line_list=line_list)
    s.prepare_file(abun_change=abun_change)
    # Whole spectra 
    s.run_moog()
    s.read_spectra()
    wav_all, flux_all = s.wav, s.flux

    # weedout lines
    if weedout_switch != False:
        w = weedout.weedout(teff, logg, m_h, line_wav_input-search_half_width, line_wav_input+search_half_width, line_list=line_list, kappa_ratio=weedout_switch)
        w.prepare_file()
        w.run_moog()
        
    # Target line exclude
    if weedout_switch:
        w.read_linelist()
        linelist_keep = w.keep_list
    else:
        linelist_keep = linelist_all
    
    line_index_keep = find_lines(linelist_keep, linelist_all)

    r_blend_ratio_list = []
    for line_index in line_index_keep:
        s = synth.synth(teff, logg, m_h, line_wav_input-search_half_width-1, line_wav_input+search_half_width+1, 
                        resolution, line_list=line_list)
        s.prepare_file(abun_change=abun_change)
        linelist_exclude = linelist_all.drop(line_index).reset_index(drop=True)
        save_linelist(linelist_exclude, s.rundir_path + 'line.list')
        s.run_moog()
        s.read_spectra(remove=False)
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
        save_linelist(linelist_target, s.rundir_path + 'line.list')
        s.run_moog()
        s.read_spectra()
        wav_target, flux_target = s.wav, s.flux

        # Calculate the EW and blending fraction
        EW = (np.sum(1-flux_all)*0.02 - np.sum(1-flux_exclude)*0.02) * 1000
        depth = 1 - np.min(flux_all[np.abs(wav_all-line_wavlength) <= 0.03])
        r_blend_ratio = (1-flux_exclude[np.argmin(np.abs(wav_exclude-line_wavlength))]) / (1-flux_all[np.argmin(np.abs(wav_all-line_wavlength))])

        r_blend_ratio_list.append(r_blend_ratio)

    linelist_keep['r_blend_depth'] = r_blend_ratio_list

    if len(line_index_keep) > 0:
        try:
            dominant_line_index = np.abs(linelist_keep.loc[linelist_keep['r_blend_depth'] < r_blen_thres, 'wavelength'] - line_wav_input).sort_values().index[0]
            dominant_line = linelist_keep.loc[dominant_line_index:dominant_line_index].reset_index(drop=True)
        except IndexError:
            # No dominant line is found
            dominant_line = pd.DataFrame(np.array([np.nan]*8)).T
            dominant_line.columns = ['wavelength', 'id', 'EP', 'loggf', 'C6', 'D0', 'EW', 'r_blend_depth']
    else:
        # No line is found
        dominant_line = pd.DataFrame(np.array([np.nan]*8)).T
        dominant_line.columns = ['wavelength', 'id', 'EP', 'loggf', 'C6', 'D0', 'EW', 'r_blend_depth']

    if include_strong_linelist:
        return dominant_line, linelist_keep
    else:
        return dominant_line