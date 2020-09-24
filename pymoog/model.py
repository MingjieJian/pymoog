#!/usr/bin/python
# This is the python model dealing with the stellar atmosphere model.
import os
import subprocess
import numpy as np
import re
import pandas as pd
from pymoog import line_data
from scipy.spatial import Delaunay

MOOG_path = '{}/.pymoog/moog_nosm/moog_nosm_NOV2019/'.format(os.environ['HOME'])
MOOG_run_path = '{}/.pymoog/rundir/'.format(os.environ['HOME'])
MOOG_file_path = '{}/.pymoog/files/'.format(os.environ['HOME'])

if os.environ.get('READTHEDOCS') != 'True':
    #directory_path = os.path.dirname(os.path.abspath(__file__))
    grid_kurucz = pd.read_csv(MOOG_file_path + '/grid_points_kurucz.csv')
    grid_matrix = np.array(grid_kurucz[['Teff', 'logg', 'm_h']])
    tri = Delaunay(grid_matrix)
    
def read_Kurucz_model(model_path):
    '''
    Read the Kurucz model and save it as np.array.
    
    Parameters
    ----------
    model_path : str
        The path of Kurucz model.
    
    Returns
    ----------
    abun : np.array
        The 'ABUNDANCE CHANGE' part of Kurucz model. The first column is element number and second column is the abundance change.
    model_lines : np.array
        The array of Kurucz model. 
    pradk : float
        The pradk number in Kurucz model.
    '''
    
    model_file = open(model_path)
    # Convert the model files into MOOG format.

    # Read and save the first two lines (except 'TITLE ') into header.
    header = ['Kurucz model: ' + model_file.readline()]
    header = header + [model_file.readline()]

    # Read the abundance change as well as model lines.
    temp = model_file.readline() + model_file.readline()

    abun_list = ''
    temp = model_file.readline()
    abun_list = abun_list + temp[42:]
    temp = model_file.readline()
    while 'ABUNDANCE CHANGE' in temp:
        abun_list = abun_list + temp[temp.index('ABUNDANCE CHANGE')+16:]
        temp = model_file.readline()
    abun = np.array(abun_list.split(), dtype='f').reshape(int(len(abun_list.split())/2), 2)

    # Read the model lines
    temp = temp.split()
    model_lines = []
    for _ in range(int(temp[2])):
        model_lines.append(model_file.readline().split())
    model_lines = np.array(model_lines, dtype=np.float64)

    # Read the PRADK value
    pradk = float(model_file.readline()[5:])
    
    return abun, model_lines, pradk

def save_interpo_model(teff, logg, m_h, abun, model_line, pradk, to_path):
    '''
    Save the array of kurucz model (e.g., the output of `read_Kurucz_model`) into a file with Kurucz format. The given stellar parameters will be written into the file.
    
    Parameters
    ----------
    teff : int
        The effective temperature of the model
    logg : float
        logg value of the model
    m_h : float
        [M/H] value (overall metallicity) of the model
    abun : np.array
        The 'ABUNDANCE CHANGE' part with the same format in the output of read_Kurucz_model.
    model_line : np.array
        The model_line part with the same format in the output of read_Kurucz_model.
    pradk : float
        pradk value.
    to_path : str
        The path to save the model.
        
    '''
    if to_path == None:
        to_path = MOOG_run_path + 'model.mod'
    else:
        pass

    content = ['Kurucz model: ' + 'TEFF   {:.1f}  GRAVITY {:.5f} LTE\n'.format(teff, logg)]
    content = content + ['TITLE SDSC GRID  [{:+.1f}]   VTURB 2.0 KM/S    L/H 1.25\n'.format(m_h)]    
    content = content + [' OPACITY IFOP 1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 0 0 0 0 0\n']
    content = content + [' CONVECTION ON   1.25 TURBULENCE OFF  0.00  0.00  0.00  0.00\n']
    content = content + [' ABUNDANCE SCALE   {:.5f} ABUNDANCE CHANGE 1 {:.5f} 2 {:.5f}\n'.format(10**m_h, *abun[0:2,1])]

    i = 2
    while i <= 98:
        if i != 98: 
            content = content + [' ABUNDANCE CHANGE  {0:2.0f} {6:6.2f}  {1:2.0f} {7:6.2f}  {2:2.0f}  {8:6.2f}  {3:2.0f} {9:6.2f}  {4:2.0f} {10:6.2f}  {5:2.0f} {11:6.2f}\n'.format(*abun[i:i+6,0], *abun[i:i+6,1])]
            i += 6
        else:
            content = content + [' ABUNDANCE CHANGE  {0:2.0f} {1:6.2f}\n'.format(abun[i,0], abun[i,1])]
            i += 1

    content = content + ['READ DECK6 {:2.0f} RHOX,T,P,XNE,ABROSS,ACCRAD,VTURB\n'.format(len(model_line))]
    
    # model_line part
    for line in model_line:
        if len(line) == 9:
            content.append(' {:.8E}  {:7.1f} {:9.3E} {:9.3E} {:9.3E} {:9.3E} {:9.3E} {:9.3E} {:9.3E}\n'.format(*line))
        elif len(line) == 7:
            content.append(' {:.8E}  {:7.1f} {:9.3E} {:9.3E} {:9.3E} {:9.3E} {:9.3E}\n'.format(*line))

    # End part
    content.append('PRADK {:.4E}\n'.format(pradk))
    content.append('BEGIN                    ITERATION  15 COMPLETED\n')
    
    with open(to_path, 'w') as file:
        file.writelines(content)

def interpolate_model(teff, logg, m_h, to_path=None, abun_change=None, kurucz_format=False, molecules=None):
    '''
    Interpolate the model in Kurucz format according to given stellar paraeters when necessary.
    
    Parameters
    ----------
    teff : int
        The effective temperature of the model
    logg : float
        logg value of the model
    m_h : float
        [M/H] value (overall metallicity) of the model
    to_path : str, optional
        The path of Kurucz model. If not given then it will be in MOOG_run_path + 'model.mod'
    abun_change : dict of pairs {int:float, ...}
        Abundance change, have to be a dict of pairs of atomic number and [X/Fe] values.   
    kurucz_format : bool, default False
        If False then the model in MOOG format will be saved; if True then the initial Kurucz format  will be saved.
    '''
    
    if to_path == None:
        to_path = MOOG_run_path + 'model.mod'
    
    p = np.array([teff, logg, m_h])
    
    # Find the grid point for interpolation and their coefficients.
    tri_simplex = tri.find_simplex(p)
    if tri_simplex == -1:
        raise ValueError('The given stellar parameters are outside grid points, failed to interpolate.')
    else:
        tri_index = tri.simplices[tri_simplex]
    grid_kurucz_sub = grid_kurucz.loc[tri_index].reset_index(drop=True)

    b = tri.transform[tri_simplex][:3].dot(np.transpose(p - tri.transform[tri_simplex][3]))
    b = np.concatenate([b, [1-sum(b)]])

    grid_kurucz_use = grid_kurucz_sub[b != 0].reset_index(drop=True)

    # Judge if the grid space is too large for interpolation
    if max(grid_kurucz_use['Teff'] >= 35000):
        teff_space_lim = 5000
    else:
        teff_space_lim = 1500
    teff_space_bad = np.ptp(grid_kurucz_use['Teff']) > teff_space_lim
    logg_space_bad = np.ptp(grid_kurucz_use['logg']) > 0.5
    m_h_space_bad = np.ptp(grid_kurucz_use['m_h']) > 0.5

    if np.any([teff_space_bad, logg_space_bad, m_h_space_bad]):
        raise ValueError('The separation between grid points is too large, failed to interpolate.')

    b = b[b != 0]

    if len(grid_kurucz_use) == 1:
        # No interpolation
        model_path = MOOG_file_path + 'model/kurucz/standard/single/teff{:.0f}logg{:.1f}m_h{:+.1f}.dat'.format(*np.array(grid_kurucz_use.loc[0]))
        subprocess.run(['cp', model_path, to_path])
        if not kurucz_format:
            KURUCZ_convert(model_path=to_path, abun_change=abun_change)
    else:
        # Interpolation
        short_64 = np.any(grid_kurucz_use['length'] == 64)
        column_7 = np.any(grid_kurucz_use['column'] == 7)
        for i in range(len(grid_kurucz_use)):
            model_path = MOOG_file_path + 'model/kurucz/standard/single/teff{:.0f}logg{:.1f}m_h{:+.1f}.dat'.format(*np.array(grid_kurucz_use.loc[i]))
            abun_single, model_line_single, pradk_single = read_Kurucz_model(model_path)

            # Cut the long model (72) into short (64) if one of the grid points model is short.
            if short_64 and len(model_line_single) == 72:
                model_line_single = model_line_single[8:]
            # Cut the large column (9) into small (7) if one of the grid points model is small.
            if column_7 and model_line_single.shape[1] == 9:
                model_line_single = model_line_single[:,:7]
            if i == 0:
                abun = abun_single * b[i]
                model_line = model_line_single * b[i]
                pradk = pradk_single * b[i]
            else:
                abun = abun + abun_single * b[i]
                model_line = model_line + model_line_single * b[i]
                pradk = pradk + pradk_single * b[i]

        # Output the interpolated model
        save_interpo_model(teff, logg, m_h, abun, model_line, pradk, to_path)
        if not kurucz_format:
            KURUCZ_convert(model_path=to_path, abun_change=abun_change, molecules=molecules)
        
def KURUCZ_convert(model_path=None, vmicro=2.0, abun_change=None, converted_model_path=None, model_type='atlas9', molecules=None):
    '''
    Convert the model file from Kurucz format in to MOOG format.

    Parameters
    ----------
    model_path : str, optional
        The path of donloaded model file. If not given then it will be MOOG_run_path + 'model.mod'
    v_micro : float, default 2.0
        microturbulance velocity of the spectra.
    abun_change : dict of pairs {int:float, ...}
        Abundance change, have to be a dict of pairs of atomic number and [X/Fe] values.
    converted_model_path : str, optional
        The name of converted model. Default will be saved into MOOG working folder.
    type : str, default 'atlas9'
        The type if input model, either 'atlas9' or 'atlas12'.
    '''
    if model_path == None:
        model_path = MOOG_run_path + 'model.mod'
    model_file = open(model_path)
    # Convert the model files into MOOG format.

    # Read and save the first two lines (except 'TITLE ') into header.
    header = model_file.readline() + model_file.readline()
    teff, logg, m_h, vmicro_model, l_h = [float(s) for s in re.findall(r'[-+]?[0-9]*\.?[0-9]+', header)]

    # Read the abundance change as well as model lines.
    temp = model_file.readline() + model_file.readline() + model_file.readline()

    abun_list = ''
    temp = model_file.readline()
    while 'ABUNDANCE CHANGE' in temp[:17]:
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
    if model_type == 'atlas12':
        for _ in range(22):
            temp = model_file.readline()
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
        c_model_path = MOOG_run_path + '/model.mod'
    else:
        c_model_path = converted_model_path
    c_model_file = open(c_model_path, 'w')

    # Header part
    c_model_file.writelines('KURUCZ\n')
    c_model_file.writelines('TEFF = {:.1f}, LOGG = {:.1f}, M/H = {:.1f}, VTURB = {:.1f}, L/H = {:.2f}\n'.format(teff, logg, m_h, vmicro_model, l_h))

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
        ele_names = list(abun_change.keys())
        ele_names.sort()
        for element in ele_names:
            if element == 2:
                c_model_file.writelines('      {:.2f}    {:.2f}\n'.format(element, xabu[element-1]+abun_change[element]))
            else:
                c_model_file.writelines('      {:.2f}    {:.2f}\n'.format(element, xabu[element-1]+abun_change[element]+float(m_h)))
    else:
        c_model_file.writelines('NATOMS      0   {}\n'.format(m_h))
    
    # Molecular line switches
    if molecules != None:
        molecules_num = len(molecules)
        c_model_file.writelines('NMOL        {}\n'.format(molecules_num))
        molecules_str = ['{:11.1f}'.format(i) for i in molecules]
        molecules_str = ''.join(molecules_str)
        c_model_file.writelines(molecules_str)
    else:
        c_model_file.writelines('NMOL        0')
    c_model_file.close()

    # cv_situation = os.path.isfile(c_model_path)
    # return c_model_file, cv_situation  