#!/usr/bin/python
from re import sub
import numpy as np
import pandas as pd
import os
import platform
from . import private
from . import synth, line_data
from . import weedout

MOOG_path = '{}/.pymoog/moog_nosm/moog_nosm_NOV2019/'.format(os.environ['HOME'])
MOOG_run_path = '{}/.pymoog/rundir/'.format(os.environ['HOME'])
MOOG_file_path = '{}/.pymoog/files/'.format(os.environ['HOME'])

## Convert the element column to element specics

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
    linelist_all = line_data.read_linelist(line_list)
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
        line_data.save_linelist(linelist_exclude, s.rundir_path + 'line.list')
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
        line_data.save_linelist(linelist_target, s.rundir_path + 'line.list')
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

def cal_d_blending_ratio(teff, logg, m_h, start_wav, end_wav, resolution, linelist_all='vald_3000_24000', weedout_switch=0.01,  abun_change=None):
    '''
    Calculate the depth blending ratio of the lines in the line list.
    The depth blending ratio is defined as d(other lines) / d(other lines + target line).  

    Parameters
    ----------
    teff : float
        The effective temperature of the model
    logg : float
        logg value of the model
    m_h : float
        [M/H] value (overall metallicity) of the model
    start_wav : float
        The start wavelength of the line list
    end_wav : float
        The end wavelength of the line list
    resolution : float
        Resolution of the synthetic spectra; this will passed to MOOG and convolute with initial spectra.
    line_list : str or pd.DataFrame, default vald_3000_24000 
        The name of the linelist file.
    weedout_switch : bool or float, default 0.01
        The switch for running weedout driver before synth. If False then weedout is not run; if True the weedout is run with kappa_ratio=0.01, and if a float (> 0 and < 1) is given then weedout is run with the kappa_ratio set as the number.
    abun_change : dict of pairs {int:float, ...}
            Abundance change, have to be a dict of pairs of atomic number and [X/Fe] values.

    Returns
    ----------
    linelist_all : pandas.DataFrame, optional
        The line list with blending ratio stored in d_blend_ratio.
    '''

    # Establish the linelist
    if isinstance(linelist_all, str):
        linelist_all = line_data.read_linelist(linelist_all)
    elif isinstance(linelist_all, private.pd.DataFrame):
        pass
    else:
        raise TypeError('Type of input linelist have to be either str or pandas.DataFrame.')
    
    linelist_all = linelist_all[(linelist_all['wavelength'] >= start_wav) & (linelist_all['wavelength'] <= end_wav)]

    # Calculate the blending ratio
    s = synth.synth(teff, logg, m_h, start_wav, end_wav, resolution, line_list=linelist_all)
    s.prepare_file(abun_change=abun_change)
    # Whole spectra 
    s.run_moog()
    s.read_spectra()
    wav_all, flux_all = s.wav, s.flux

    # weedout lines
    if weedout_switch != False:
        w = weedout.weedout(teff, logg, m_h, start_wav, end_wav, line_list=linelist_all, kappa_ratio=weedout_switch)
        w.prepare_file()
        w.run_moog()

    # Target line exclude
    if weedout_switch:
        w.read_linelist()
        linelist_keep = w.keep_list
        line_index_keep = find_lines(linelist_keep, linelist_all)
    else:
        line_index_keep = linelist_all.index

    d_blend_ratio_list = []
    for line_index in linelist_all.index:
        
        if line_index in line_index_keep:
            target_wav = linelist_all.loc[line_index, 'wavelength']
            linelist_exclude = linelist_all.drop(line_index).reset_index(drop=True)
            s = synth.synth(teff, logg, m_h, target_wav-1, target_wav+1, resolution, line_list=linelist_exclude)
            s.prepare_file(abun_change=abun_change)
            s.run_moog()
            s.read_spectra(remove=False)
            wav_exclude, flux_exclude = s.wav, s.flux
            
            # Calculate the blending fraction
            r_blend_ratio = (1-flux_exclude[np.argmin(np.abs(wav_exclude-target_wav))]) / (1-flux_all[np.argmin(np.abs(wav_all-target_wav))])
            if r_blend_ratio > 1 and np.abs((1-flux_exclude[np.argmin(np.abs(wav_exclude-target_wav))])-(1-flux_all[np.argmin(np.abs(wav_all-target_wav))])) < 0.001:
                r_blend_ratio = 1
            d_blend_ratio_list.append(r_blend_ratio)
        else:
            d_blend_ratio_list.append(np.nan)

    linelist_all['d_blend_ratio'] = d_blend_ratio_list

    return linelist_all