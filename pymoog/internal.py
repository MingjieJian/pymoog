#!/usr/bin/python
# Internal functions for renewing the database of stellar atmosphere model and linlist.
# This file is split to two part: line list part and model part.
# WARNING: the dependene in this module may not be completly satisified, and functions may can only run on Mingjie's computer.

import numpy as np
import pandas as pd
import os, pickle
from . import model
import matplotlib.pyplot as plt
from scipy.spatial import Delaunay
from pymoog import line_data
import mendeleev
import astropy.units as u
import re
from astropy import constants

MOOG_path = '{}/.pymoog/moog_nosm/moog_nosm_NOV2019/'.format(os.environ['HOME'])
MOOG_run_path = '{}/.pymoog/rundir/'.format(os.environ['HOME'])
MOOG_file_path = '{}/.pymoog/files/'.format(os.environ['HOME'])

# Line list part

element2index_dict = {'TiO':[22,8], 'CH':[6,1], 'OH':[8,1], 'MgH':[12,1], 'SiH':[14,1], 'C2':[6,6], 'CN':[6,7], 'CO':[6,8], 'NH':[7,1]}
atoms_number = pd.read_csv(MOOG_file_path + '/atoms.csv')
if os.environ.get('READTHEDOCS') != 'True':
    atoms = pd.read_csv(MOOG_file_path + '/atoms.csv')
    atoms_dict = dict(zip(atoms['symbol'], atoms['mass_number']))
    diss_energy = pd.read_csv(MOOG_file_path + '/dissociation_energy_list.csv')

def get_isotope_list(string):
    '''
    Get the isotope list of element from the last column of VALD line list. 
    
    Example:  (48)TiO -> [48, 16]
    
    Parameters
    ----------
    string : str
        The string in the format of "(\d*)[A-Z][a-z]*". This is the last part in VALD linelist.
    '''
    a = re.findall(r'\(\d*\)|[A-Z][a-z]*', string)
    isotope_list = []
    i = 0
    while i < len(a):
        if a[i][0] == '(':
            isotope_list.append(int(a[i].strip('()')))
            i += 1
        else:
            isotope_list.append(atoms_dict[a[i]])
        i += 1
    return isotope_list

def element2index(string_all):
    '''
    Convert element string to index in VALD format. 
    
    Example: TiO 1, ... (48)TiO -> 822.01648; Fe 1, ... Fe -> 26.0.
    
    Parameters
    ----------
    string_all : str
        The string in containing element index in VALD linelist. Combination of the first and last column.
    '''   
    
    # Split the string
    string, isotope_string = string_all.split(',')
    # isotope_string = isotope_string[-12:]
    isotope_string = isotope_string.split()[-1]
    element_string, ion_stage = string.split(' ')

    if element_string in element2index_dict.keys():
        element_indices = element2index_dict[element_string]
    else:
        element_indices = [atoms.loc[atoms['symbol'] == element_string, 'atom_number'].values[0]]
        ion_stage = int(ion_stage) - 1
        # p = re.compile(r"[A-Z][a-z]*")
        # p_num = re.compile(r"\d")
        # ele_loca = []
        # ele_name = []
        # num_loca = []
        # num = []
        # for m in p.finditer(element_string):
        #     ele_loca.append(m.start())
        #     ele_name.append(m.group())

        # for m in p_num.finditer(element_string):
        #     num_loca.append(m.start())
        #     num.append(m.group())
        # element_string_list = []
        # for i in range(len(ele_name)):
        #     if ele_loca[i]+1 in num_loca:
        #         add_list = [ele_name[i]] * int(num[num_loca.index(ele_loca[i]+1)])
        #     else:
        #         add_list = [ele_name[i]]
        #     element_string_list = element_string_list + add_list

        # element_indices = []
        # for ele in element_string_list:
        #     element_indices.append(mendeleev.element(ele).atomic_number)
   
    if len(element_indices) == 1:
        return '{}.{:05.0f}'.format(element_indices[0], ion_stage*10000)
    else:
        isotope_list = get_isotope_list(isotope_string)
        # isotope_list = [x for _,x in sorted(zip(element_indices,isotope_list))]
        element_indices.sort()
        isotope_list.sort()
        element_indices_string = '{:2.0f}' + '{:02.0f}'*(len(element_indices)-1) + '.0' + '{:02.0f}'*len(isotope_list)
        element_indices_num = float(element_indices_string.format(*element_indices, *isotope_list))
        return '{:.5f}'.format(element_indices_num)
        # return element_indices_string.format(*element_indices, *isotope_list)

def get_diss_energy(ele_index):
    '''
    Get dissociation for an molecular particle from ele_index.
    Source: https://labs.chem.ucsb.edu/zakarian/armen/11---bonddissociationenergy.pdf
    
    Only support those in VALD linelist.
    
    Parameters
    ----------
    ele_index : str or float
        The element index in MOOG format.
    '''
    diss_energy['diss_energy(eV)'] = diss_energy['dissociation_energy (kJ/mol)'] / 96.485
    diss_energy_pd = diss_energy
    ele_index = np.floor(float(ele_index))
    try:
        diss_energy_value = diss_energy_pd.loc[diss_energy_pd['element_index'] == ele_index, 'diss_energy(eV)'].values[0]
        return diss_energy_value
    except: 
        return np.nan

def vald2moog_format(init_linelist_name, out_linelist_name, loggf_cut=None):
    '''
    Transform VALD linelist download from VALD website into MOOG format.
    
    Parameters
    ----------
    init_linelist_name : str
        The VALD format line list.
    out_linelist_name : str
        Output line list name
    head : int, optional
        If specified then only save the first `head` number of lines.
    loggf_cut : float, optional
        Cut on loggf (only save for the lines with loggf > loggf_cut)
    '''
    # Find the footer index of VALD line pair
    with open(init_linelist_name) as file:
        contents = file.readlines()
        try:
            footer_index = len(contents) - contents.index('* oscillator strengths were scaled by the solar isotopic ratios.\n')
        except ValueError:
            footer_index = 0
 
    # Delete all the '.
    file = open(init_linelist_name)
    file_content = file.readlines()
    for i in range(len(file_content)):
        file_content[i] = file_content[i].replace("'", '')
    file.close()
    file = open(init_linelist_name, 'w')
    file.writelines(file_content)
    file.close()
    
    # subprocess.run(['sed', "s/'//g", init_linelist_name, '>', 'temp'])
    # subprocess.run(['mv', "temp", init_linelist_name])    
    vald_init = pd.read_csv(init_linelist_name, skiprows=2, skipfooter=footer_index, usecols=range(9), engine = 'python', names=['element', 'wavelength', 'EP', 'loggf', 'rad_damp', 'Stark_damp', 'Walls_damp', 'Lande_factor', 'Comment'])

    if loggf_cut != None:
        vald_init = vald_init[vald_init['loggf'] >= loggf_cut]
    vald_init['element_all'] = vald_init[['element', 'Comment']].apply(lambda x: ', '.join(x), axis=1)
    vald_init['element_index'] = vald_init['element_all'].map(element2index)
    vald_init['diss_energy'] = vald_init['element_index'].map(get_diss_energy)

    vald_out = vald_init[['wavelength', 'element_index', 'EP', 'loggf', 'Walls_damp', 'diss_energy']]
    vald_out.columns = ['wavelength', 'element_index', 'EP', 'loggf', 'C6', 'diss_energy']
    vald_out = vald_out.astype(np.float64)
    
    # Remove triple or higher ionized lines; MOOG cannot do this.
    vald_out = vald_out[np.around(np.mod(vald_out['element_index'],1), decimals=1) < 0.3]
    
    line_data.save_linelist(vald_out, out_linelist_name)

def ges2moog(ges_path, save_path):
    GES = pd.read_csv(ges_path, sep='\t')
    GES['diss_energy'] = np.nan
    GES = GES[GES['moog_support'] == 'T']

    GES_moog = GES[['wave_A', 'spectrum_moog_species', 'lower_state_eV', 'loggf', 'waals', 'diss_energy', 'theoretical_ew']]
    GES_moog.columns = ['wavelength', 'element_index', 'EP', 'loggf', 'C6', 'D0', 'theoretical_ew']
    
    line_data.save_linelist(GES_moog, save_path, header='MB99 linelist')
    

def ele2ele_num(string):
    ion_dict = {'I':1, 'II':2}
    str_split = string.split(' ')
    return mendeleev.element(str_split[0]).atomic_number + ion_dict[str_split[1]] / 10

def mb992moog(mb99_path, save_path):
    mb99_j = pd.read_fwf(mb99_path, colspecs=[(0,7), (8,16), (24,29), (32,37), (48,56)], names=['ele', 'wavelength', 'EP', 'loggf', 'C6'])
    mb99_j['ele'] = mb99_j['ele'].map(ele2ele_num)
    mb99_j['D0'] = np.nan
    mb99_j['EW'] = np.nan

    mb99_j_out = mb99_j[['wavelength', 'ele', 'EP', 'loggf', 'C6', 'D0', 'EW']]
    line_data.save_linelist(mb99_j_out, save_path, header='MB99 linelist')
    
def kurucz2moog(kurucz_path):
    # gfall08oct17.dat
    kurucz_all = pd.read_fwf('files/linelist/kurucz/gfall08oct17.dat', colspecs=[(0,11), (11,18), (18,24), (24,36), (52,64), (93,98),  (109, 116)], names=['wavelength', 'loggf_init', 'ele', 'E(cm-1)_1', 'E(cm-1)_2', 'C6', 'hpf_frac'])
    kurucz_all['ele'] = kurucz_all['ele'] // 1 + kurucz_all['ele'] % 1 * 10 
    kurucz_all['EP_1'] = kurucz_all['E(cm-1)_1'] / (1/constants.h / constants.c).to(u.cm**-1 / u.eV).value
    kurucz_all['EP_2'] = kurucz_all['E(cm-1)_2'] / (1/constants.h / constants.c).to(u.cm**-1 / u.eV).value
    kurucz_all['EP'] = np.where((kurucz_all['EP_1'] <= kurucz_all['EP_2']), kurucz_all['EP_1'], kurucz_all['EP_2'])
    kurucz_all['loggf'] = kurucz_all['loggf_init'] + kurucz_all['hpf_frac']
    indices = (kurucz_all['ele'] % 1 <= 0.2) & ~np.isnan(kurucz_all['loggf']) & (kurucz_all['wavelength'] >= 200) & (kurucz_all['EP'] <= 50)    
    kurucz_use = kurucz_all.loc[indices, ['wavelength', 'ele', 'EP', 'loggf', 'C6']].reset_index(drop=True)
    kurucz_use['wavelength'] = kurucz_use['wavelength'] * 10
    kurucz_use['D0'] = np.nan
    kurucz_use['EW'] = np.nan
    line_data.save_linelist(kurucz_use, 'files/linelist/kurucz/kurucz.list', wav_start=2000, wav_end=7e5)
    return kurucz_use
    
def get_species(num):
    if num <= 100:
        return num
    else:
        num_str = '{:011.6f}'.format(num)
        atom_1, atom_2, iso_1, iso_2 = int(num_str[:2]), int(num_str[2:4]), int(num_str[6:8]), int(num_str[9:])
        if iso_1 == 0:
            iso_1 = atoms_number.loc[atoms_number['atom_number'] == atom_1, 'mass_number'].values[0]
        if iso_2 == 0:
            iso_2 = atoms_number.loc[atoms_number['atom_number'] == atom_2, 'mass_number'].values[0]
    return atom_1*100 + atom_2 + iso_1 / 1000 + iso_2 / 100000


# Model part
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
            if 'EFF ' in line:
                if not(is_first):
                    with open('files/model/kurucz/standard/single/teff{:.0f}logg{:.1f}m_h{:+.1f}.dat'.format(teff, logg, m_h), 'w') as w_file:
                        w_file.writelines(model_line)
                teff = float(line[5:13])
                logg = float(line[21:29])  
                model_line = [line]
                is_first = False
            else:
                model_line.append(line)
        with open('files/model/kurucz/standard/single/teff{:.0f}logg{:.1f}m_h{:+.1f}.dat'.format(teff, logg, m_h), 'w') as w_file:
            w_file.writelines(model_line)

def search_grid_point_kurucz():
    '''
    The function to search all the grid points of Kurucz model and save the list to grid_path.
    The search is limit to standard model with microturbulent = 2.
    Internal use.
    '''
    teff_range = np.arange(3500, 50001, 250)
    logg_range = np.arange(0, 5.1, 0.5)
    m_h_range = np.concatenate([np.arange(-5, -0.4, 0.5), np.arange(-0.3, 0, 0.1), [0], np.arange(0.1, 0.35, 0.1) ,[0.5, 1]])

    grid_point_kurucz = []
    for m_h in m_h_range:
        for teff in teff_range:
            for logg in logg_range:
                if os.path.isfile('files/model/kurucz/standard/single/teff{:.0f}logg{:.1f}m_h{:+.1f}.dat'.format(teff, logg, m_h)):
                    _, b, _ = model.read_Kurucz_model('files/model/kurucz/standard/single/teff{:.0f}logg{:.1f}m_h{:+.1f}.dat'.format(teff, logg, m_h))
                    length = b.shape[0]
                    column = b.shape[1]
                    if len(grid_point_kurucz) == 0:
                        grid_point_kurucz = np.array([[teff, logg, m_h, length, column]])
                    else:
                        grid_point_kurucz = np.concatenate([grid_point_kurucz, np.array([[teff, logg, m_h, length, column]])])
    grid_kurucz = pd.DataFrame(grid_point_kurucz, columns=['Teff', 'logg', 'm_h', 'length', 'column'])
    return grid_kurucz

def plot_model_grid():
    '''
    Plot the grid of models in each metallicity.
    Internal use.
    '''
    
    grid_kurucz = pd.read_csv('files/grid_points_kurucz.csv')
    for m_h in grid_kurucz.groupby('m_h').size().index:
        plt.figure(figsize=(13,4))
        index = grid_kurucz['m_h'] == m_h
        grid_matrix = np.array(grid_kurucz.loc[index, ['Teff', 'logg']])
        tri = Delaunay(grid_matrix)
        for i in range(len(tri.simplices)-1, -1, -1):
            if min(grid_matrix[tri.simplices[i]][:,0]) >= 35000:
                teff_gap = 5000
            else:
                teff_gap = 1500
            if  np.ptp(grid_matrix[tri.simplices[i]][:,0]) >= teff_gap or np.ptp(grid_matrix[tri.simplices[i]][:,1]) > 0.5:
                tri.simplices = np.concatenate([tri.simplices[:i], tri.simplices[i+1:]])

        plt.triplot(grid_matrix[:,0], grid_matrix[:,1], tri.simplices, zorder=0, lw=1, color='gray',alpha=0.5)
        if m_h < 0.5:
            plt.plot([50000, 42500], [5, 5], color='gray', zorder=0, alpha=0.5, lw=1)
        elif m_h == 0.5:
            plt.plot([45000, 40000], [5, 5], color='gray', zorder=0, alpha=0.5, lw=1)
        elif m_h == 1:
            plt.plot([40000, 37500], [5, 5], color='gray', zorder=0, alpha=0.5, lw=1)

        
        plt.scatter(grid_kurucz.loc[index & (grid_kurucz['length']==72), 'Teff'], grid_kurucz.loc[index & (grid_kurucz['length']==72), 'logg'], s=5, label='Model length: 72')
        plt.scatter(grid_kurucz.loc[index & (grid_kurucz['length']==64), 'Teff'], grid_kurucz.loc[index & (grid_kurucz['length']==64), 'logg'], s=5, c='C3', label='Model length: 64')

        plt.legend()
        plt.xlim((1175, 52325))
        plt.title('[Fe/H] = {:.1f}'.format(m_h))
        plt.xlabel(r'$T_\mathrm{{eff}}$'); plt.ylabel('logg')
        plt.gca().invert_xaxis(); plt.gca().invert_yaxis()
        plt.tight_layout()
        plt.savefig('../docs/img/grid_points_kurucz/m_h{:+.1f}.png'.format(m_h), dpi=250)
        plt.close()

# Here we calculate the Delaunay triangulation for MRACS models, and sotre them into corredponding folders.
def save_marcs_delaunay(folder_path, tri_columns=['teff', 'g', '[M/H]', 'vmicro', 'mass']):
    '''
    Generate the DataFrame of grid points and Delaunay triangulation for MARCS models, and save it to the same folder. For internal use.
    
    Parameters
    ----------
    folder_path : str
        The folder path for generating the grid points DataFrame and Delaunay triangulation.
    tri_columns : list
        List of columns to be used for Delaunay triangulation.
    '''
    marcs_name_list = os.listdir(folder_path)
    marcs_name_list = [i for i in marcs_name_list if 'mod' in i]
    marcs_grid = []

    for marcs_model in marcs_name_list:
        marcs_grid.append(list(model.marcs_filename2para(marcs_model)))

    marcs_grid = pd.DataFrame(marcs_grid, columns=['chem', 'geo', 'teff', 'logg', 'mass', 'vmicro', '[M/H]', '[alpha/Fe]', 
                                                           '[C/Fe]', '[N/Fe]', '[O/Fe]', 'r', 's'])    
    marcs_grid['g'] = 10**marcs_grid['logg']
    
    # Sanity check: chem and geo have to be the same in the dataframe.
    if len(marcs_grid.groupby(['chem', 'geo']).size()) > 1:
        raise ValueError('Chemical composition or model geometry is not unique in the folder. Please limit them to only one.')
        
    grid_matrix = np.array(marcs_grid[tri_columns])
    tri = Delaunay(grid_matrix)
    
    marcs_grid.to_csv(folder_path + '/grid_points.csv', index=False)
    pickle.dump(tri, open(folder_path + '/tri.pkl', 'wb'))
    
    pass

def main():
   init_linelist_name = sys.argv[1]
   out_linelist_name = sys.argv[2]
   vald2moog_format(init_linelist_name, out_linelist_name)

if __name__ == "__main__":
    main()