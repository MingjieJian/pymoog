import numpy as np
import pandas as pd
import mendeleev as md
import re
import sys

# Convert the element column to element specics

element2index_dict = {'TiO':[22,8], 'CH':[6,1], 'OH':[8,1], 'MgH':[12,1], 'SiH':[14,1], 'C2':[6,6], 'CN':[6,7], 'CO':[6,8]}
atoms = pd.read_csv('files/atoms.csv')
atoms_dict = dict(zip(atoms['symbol'], atoms['mass_number']))

def get_isotope_list(string):
    '''
    Get the isotope list of element from the last column of VALD line list.
    (48)TiO -> [48, 16]
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
    string, isotope_string = string_all.split(',')
    isotope_string = isotope_string[-12:]
    element_string, ion_stage = string.split(' ')

    if element_string in element2index_dict.keys():
        element_indices = element2index_dict[element_string]
    else:
        p = re.compile(r"[A-Z][a-z]*")
        p_num = re.compile(r"\d")
        ele_loca = []
        ele_name = []
        num_loca = []
        num = []
        for m in p.finditer(element_string):
            ele_loca.append(m.start())
            ele_name.append(m.group())

        for m in p_num.finditer(element_string):
            num_loca.append(m.start())
            num.append(m.group())
        element_string_list = []
        for i in range(len(ele_name)):
            if ele_loca[i]+1 in num_loca:
                add_list = [ele_name[i]] * int(num[num_loca.index(ele_loca[i]+1)])
            else:
                add_list = [ele_name[i]]
            element_string_list = element_string_list + add_list

        ion_stage = int(ion_stage)
        element_indices = []
        for ele in element_string_list:
            element_indices.append(md.element(ele).atomic_number)
   
    if len(element_indices) == 1:
        return '{}.{}'.format(element_indices[0], ion_stage*10000)
    else:
        isotope_list = get_isotope_list(isotope_string)
        # isotope_list = [x for _,x in sorted(zip(element_indices,isotope_list))]
        element_indices.sort()
        isotope_list.sort()
        element_indices_string = '{:2.0f}' + '{:02.0f}'*(len(element_indices)-1) + '.0' + '{:02.0f}'*len(isotope_list)
        element_indices_num = float(element_indices_string.format(*element_indices, *isotope_list))
        return element_indices_num
        # return element_indices_string.format(*element_indices, *isotope_list)

def get_diss_energy(ele_index):
    '''
    Get dissociation for an molecular particle from ele_index.
    '''
    
    diss_energy = pd.read_csv('./files/dissociation_energy_list.csv')
    diss_energy['diss_energy(eV)'] = diss_energy['dissociation_energy (kJ/mol)'] / 96.485
    diss_energy_pd = diss_energy
    ele_index = np.floor(float(ele_index))
    try:
        diss_energy_value = diss_energy_pd.loc[diss_energy_pd['element_index'] == ele_index, 'diss_energy(eV)'].values[0]
        return diss_energy_value
    except: 
        return np.nan
    
def save_linelist(linelist_all, sub_ll_name, wav_start=None, wav_end=None, type='vald'):
    if type != 'vald':
        raise Exception('Only vald linelist is supported!')
    
    # Crop the line list according to wavelength, if needed.
    index = linelist_all['wavelength'] > 0
    if wav_start != None:
        index = index & (linelist_all['wavelength'] > wav_start)
    if wav_end != None:
        index = index & (linelist_all['wavelength'] < wav_end) 
           
    sub_linelist = linelist_all[index]
    sub_linelist.reset_index(drop=True, inplace=True)
    with open(sub_ll_name, 'w') as file:
        file.write('VALD linelist\n')
        for i in range(len(sub_linelist)):
            if np.isnan(sub_linelist.iloc[i].values[-1]):
                file.write('{:10.4f}{:10.5f}{:10.4f}{:10.3f}{:10.3f}\n'.format(*sub_linelist.iloc[i].values[:-1]))
            else:
                file.write('{:10.4f}{:10.5f}{:10.4f}{:10.3f}{:10.3f}{:10.3f}\n'.format(*sub_linelist.iloc[i].values))

def read_linelist(linelist_path, loggf_cut=None):
    '''
    Read the post-processed linelist.
    '''
    linelist = pd.read_fwf(linelist_path,
            colspecs=[(0,11), (11,21), (21,31), (31,41), (41,51), (51,61)],
            names=['wavelength', 'id', 'EP', 'loggf', 'C6', 'D0'],
            skiprows=1)
    if loggf_cut != None:
        linelist = linelist[linelist['loggf'] >= loggf_cut]
        linelist.reset_index(drop=True, inplace=True)
    return linelist
    
def vald2moog_format(init_linelist_name, out_linelist_name, head=None, loggf_cut=None):
    '''
    Transform VALD linelist into MOOG format for a given wavlength range.
    '''
    # Find the footer index of VALD line pair
    with open(init_linelist_name) as file:
        contents = file.readlines()
        try:
            footer_index = len(contents) - contents.index('* oscillator strengths were scaled by the solar isotopic ratios.\n')
        except ValueError:
            footer_index = 0
 
    vald_init = pd.read_csv(init_linelist_name,skiprows=2, skipfooter=footer_index, usecols=range(9), engine = 'python' ,
                                names=['element', 'wavelength', 'EP', 'loggf', 'rad_damp', 'Stark_damp', 'Walls_damp', 'Lande_factor', 'Comment'])
    
    if head != None:
        vald_init = vald_init[:head]
    if loggf_cut != None:
        vald_init = vald_init[vald_init['loggf'] >= loggf_cut]
    vald_init['element_all'] = vald_init[['element', 'Comment']].apply(lambda x: ', '.join(x), axis=1)
    vald_init['element_index'] = vald_init['element_all'].map(element2index)
    vald_init['diss_energy'] = vald_init['element_index'].map(get_diss_energy)

    vald_out = vald_init[['wavelength', 'element_index', 'EP', 'loggf', 'Walls_damp', 'diss_energy']]
    vald_out = vald_out.astype(np.float64)
    
    # Remove triple or higher ionized lines; MOOG cannot do this.
    vald_out = vald_out[np.around(np.mod(vald_out['element_index'],1), decimals=1) < 0.3]
    
    save_linelist(vald_out, out_linelist_name)

def main():
   init_linelist_name = sys.argv[1]
   out_linelist_name = sys.argv[2]
   vald2moog_format(init_linelist_name, out_linelist_name)

if __name__ == "__main__":
    main()
