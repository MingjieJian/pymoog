import numpy as np
import pandas as pd
import mendeleev as md
import re
import sys

def main():
   init_linelist_name = sys.argv[1]
   out_linelist_name = sys.argv[2]

if __name__ == "__main__":
    main()

# Convert the element column to element specics

def get_isotope_list(string):
    '''
    Get the isotope list of element from the last column of VALD line list.
    (48)TiO -> 01648
    '''
    a = re.findall('\(\d*\)|[A-Z][a-z]*', string)
    isotope_list = []
    i = 0
    while i < len(a):
        if a[i][0] == '(':
            isotope_list.append(int(a[i].strip('()')))
            i += 1
        else:
            isotope_list.append(md.element(a[i]).mass_number)
        i += 1
    return isotope_list

def element2index(string_all):
    # print(string_all)
    string, isotope_string = string_all.split(',')
    isotope_string = isotope_string[-12:]
    element_string, ion_stage = string.split(' ')

    p = re.compile("[A-Z][a-z]*")
    p_num = re.compile("\d")
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
    sorted_index = sorted(range(len(element_indices)), key=lambda k: element_indices[k])
   
    if len(element_indices) == 1:
        return '{}.{}'.format(element_indices[0], ion_stage*10000)
    else:
        isotope_list = get_isotope_list(isotope_string)
        isotope_list = [x for _,x in sorted(zip(element_indices,isotope_list))]
        element_indices.sort()
        isotope_list = list(dict.fromkeys(isotope_list))
        element_indices_string = '{:2.0f}' + '{:02.0f}'*(len(element_indices)-1) + '.0' + '{:02.0f}'*len(isotope_list)
        return element_indices_string.format(*element_indices, *isotope_list)

def get_diss_energy(ele_index, diss_energy_pd=diss_energy):
    ele_index = np.floor(float(ele_index))
    try:
        diss_energy_value = diss_energy_pd.loc[diss_energy_pd['element_index'] == ele_index, 'diss_energy(eV)'].values[0]
        return diss_energy_value
    except: 
        return np.nan
    
    
def vald2moog_format(init_linelist_name):
    '''
    Transform VALD linelist into MOOG format.
    '''

    vald_init = pd.read_csv(init_linelist_name,skiprows=2, skipfooter=107, usecols=range(9), names=['element', 'wav', 'EP', 'loggf', 'rad_damp', 'Stark_damp', 'Walls_damp', 'Lande_factor', 'Comment'])
    vald_init = vald_init[:50]
    vald_init['element_all'] = vald_init[['element', 'Comment']].apply(lambda x: ', '.join(x), axis=1)

    vald_init['element_index'] = vald_init['element_all'].map(element2index)
    vald_init['diss_energy'] = vald_init['element_index'].map(get_diss_energy)

    vald_out = vald_init[['wav', 'element_index', 'EP', 'loggf', 'Walls_damp', 'diss_energy']]

    return vald_out

diss_energy = pd.read_csv('pymoog/files/dissociation_energy_list.csv')
diss_energy['diss_energy(eV)'] = diss_energy['dissociation_energy (kJ/mol)'] / 96.485