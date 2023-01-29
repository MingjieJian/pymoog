#!/usr/bin/python
# This is the python model dealing with the stellar atmosphere model.
import os
import subprocess
import numpy as np
import re
import pandas as pd
from pymoog import line_data
from scipy.spatial import Delaunay
import pickle

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
    abun_list = abun_list + temp[42:].replace('E', '')
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

def save_kurucz_interpo_model(teff, logg, m_h, abun, model_line, pradk, to_path):
    '''
    Save the array of kurucz model (e.g., the output of `read_Kurucz_model`) into a file in Kurucz format. The given stellar parameters will be written into the file.
    
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
    content = content + ['TITLE SDSC GRID  [{:+.2f}]   VTURB 2.0 KM/S    L/H 1.25\n'.format(m_h)]    
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

def interpolate_kurucz_model(teff, logg, m_h, vmicro=2, abun_change=None, molecules_include=None, kurucz_format=False, abun_scale='kurucz', to_path=None):
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
    vmicro : float, default 2
        The microtrubulance velocity of the synthesized spectra (this is different from the v_micro of the atmosphere model which is always 2)
    abun_change : dict of pairs {int:float, ...}
        Abundance change, have to be a dict of pairs of atomic number and [X/Fe] values.   
    molecules_include : list
        Molecules to be included to molecular calculation. Follows the MOOG notation.
    kurucz_format : bool, default False
        If False then the model in MOOG format will be saved; if True then the initial Kurucz format will be saved.
    abun_scale : str, default "kurucz"
        The abundance scale of metallicity. If not 'kurucz', then will use Asplund 2009.
    to_path : str, optional
        The path of Kurucz model. If not given then it will be in MOOG_run_path + 'model.mod'
    '''
    
    if to_path == None:
        to_path = MOOG_run_path + 'model.mod'
    if abun_scale == 'kurucz':
        m_h_input = m_h
        m_h = m_h - 0.17
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
        model_path = MOOG_file_path + '/pymoog_lf/model/kurucz/standard/single/teff{:.0f}logg{:.1f}m_h{:+.1f}.dat'.format(*np.array(grid_kurucz_use.loc[0]))
        subprocess.run(['cp', model_path, to_path])
        if not kurucz_format:
            kurucz2moog(model_path=to_path, vmicro=vmicro, abun_change=abun_change, converted_model_path=to_path)
    else:
        # Interpolation
        short_64 = np.any(grid_kurucz_use['length'] == 64)
        column_7 = np.any(grid_kurucz_use['column'] == 7)
        for i in range(len(grid_kurucz_use)):
            model_path = MOOG_file_path + '/pymoog_lf/model/kurucz/standard/single/teff{:.0f}logg{:.1f}m_h{:+.1f}.dat'.format(*np.array(grid_kurucz_use.loc[i]))

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
        if to_path == False:
            return abun, model_line, pradk
        
        # Output the interpolated model
        if abun_scale == 'kurucz':
            save_kurucz_interpo_model(teff, logg, m_h_input, abun, model_line, pradk, to_path)
        else:
            save_kurucz_interpo_model(teff, logg, m_h, abun, model_line, pradk, to_path)
        if not kurucz_format:
            if abun_scale == 'kurucz':
                kurucz2moog(model_path=to_path, vmicro=vmicro, abun_change=abun_change, molecules_include=molecules_include, m_h_model=m_h, converted_model_path=to_path)
            else:
                kurucz2moog(model_path=to_path, vmicro=vmicro, abun_change=abun_change, molecules_include=molecules_include, converted_model_path=to_path)
        
def kurucz2moog(model_path=None, vmicro=2.0, abun_change=None, converted_model_path=None, model_format='atlas9', molecules_include=None, m_h_model=None):
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
    molecules_include : list
        Molecules to be included to molecular calculation. Follows the MOOG notation.
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
    if model_format == 'atlas12':
        for _ in range(22):
            temp = model_file.readline()
    temp = temp.split()
    model_lines = []
    model_linen = int(temp[2])
    for i in range(int(temp[2])):
        model_lines.append(model_file.readline().split()[:7])

    # Prepare the microtrubulance value.
    vmicro = '{:.2E}'.format(vmicro)

    # Write the model file.
    # header, abun09, model_lines and vmicro
    if converted_model_path == None:
        c_model_path = MOOG_run_path + '/model.mod'
    else:
        c_model_path = converted_model_path
    c_model_file = open(c_model_path, 'w')

    # Header part
    c_model_file.writelines('KURUCZ\n')
    if m_h_model is None:
        c_model_file.writelines('TEFF = {:.1f}, LOGG = {:.1f}, M/H = {:.2f}, VTURB = {:.1f}, L/H = {:.2f}\n'.format(teff, logg, m_h, vmicro_model, l_h))
    else:
        c_model_file.writelines('TEFF = {:.1f}, LOGG = {:.1f}, M/H = {:.2f} ({:.2f}), VTURB = {:.1f}, L/H = {:.2f}\n'.format(teff, logg, m_h, m_h_model, vmicro_model, l_h))
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
    if molecules_include != None:
        molecules_num = len(molecules_include)
        c_model_file.writelines('NMOL        {}\n'.format(molecules_num))
        molecules_str = ['{:11.1f}'.format(i) for i in molecules_include]
        molecules_str = ''.join(molecules_str)
        c_model_file.writelines(molecules_str)
    else:
        c_model_file.writelines('NMOL        0')
    c_model_file.close()

def read_marcs_model(model_path):
    '''
    Read the MARCS model file and save it as dict.
    
    Parameters
    ----------
    model_path : str
        Name of the MARCS model file.
        
    Returns
    ----------
    marcs_model : dict
        The dict of input MARCS model.
    '''
    
    geo_style = {'s':'spherical', 'p':'plane-parallal'}
    marcs_model = {}
    
    with open(model_path) as file:
        model_contents = file.readlines()
    model_contents = [i[:-1] for i in model_contents]

    # Read the model header
    marcs_model['model name'] = model_contents[0]
    marcs_model['geometry'] = geo_style[marcs_model['model name'][0]]

    marcs_model['teff'] = float(model_contents[1].split()[0])
    marcs_model['last iteration'] = model_contents[1].split()[-1][-8:]
    marcs_model['flux'] = float(model_contents[2].split()[0])
    marcs_model['g'] = float(model_contents[3].split()[0])
    marcs_model['vmicro'] = float(model_contents[4].split()[0])
    marcs_model['mass'] = float(model_contents[5].split()[0])
    marcs_model['[M/H]'] = float(model_contents[6].split()[0])
    marcs_model['[alpha/Fe]'] = float(model_contents[6].split()[1])
    marcs_model['radius'] = float(model_contents[7].split()[0])
    marcs_model['luminosity'] = float(model_contents[8].split()[0])
    marcs_model['conv:alpha'], marcs_model['conv:nu'], marcs_model['conv:y'], marcs_model['conv:beta'] = [float(i) for i in model_contents[9].split()[:4]]
    marcs_model['X'], marcs_model['Y'], marcs_model['Z'] = [float(i) for i in model_contents[10].split()[:3]]
    marcs_model['12C/13C'] = float([i for i in model_contents[10].split() if i[:7] == '12C/13C'][0].split('=')[1])
    marcs_model['log abundance'] = np.fromstring(''.join(model_contents[12:22]), sep=' ')
    marcs_model['number of depth points'] = int(model_contents[22].split()[0])

    ms_1 = pd.read_csv(model_path, skiprows=24, skipfooter=57*4+1, sep=' +', engine='python')
    ms_2 = pd.read_csv(model_path, skiprows=24+57, skipfooter=57*3+1, sep=' +', engine='python').drop(['k', 'lgTauR'], axis=1)
    marcs_model['model structure'] = pd.concat([ms_1, ms_2], axis=1)

    pp_1 = pd.read_csv(model_path, skiprows=24+57*2+1, skipfooter=57*2, sep=' +', engine='python')
    pp_1.columns = ['k', 'lgPgas', 'H I', 'H-', 'H2', 'H2+', 'H2O', 'OH', 'CH', 'CO', 'CN', 'C2', 'remove']
    pp_1 = pp_1.drop('remove', axis=1)
    pp_2 = pd.read_csv(model_path, skiprows=24+57*3+1, skipfooter=57, sep=' +', engine='python').drop('k', axis=1)
    pp_3 = pd.read_csv(model_path, skiprows=24+57*4+1, skipfooter=0, sep=' +', engine='python').drop('k', axis=1)
    marcs_model['logarithmic partial pressures'] = pd.concat([pp_1, pp_2, pp_3], axis=1)
    
    return marcs_model

def save_marcs_model(marcs_model, save_name):
    
    '''
    Save the dict of MARCS model to file.
    
    Parameters
    ----------
    marcs_model : dict
        The MARCS model dict.
    save_name : str
        The name to save the model.
    '''
    
    marcs_model_content = [marcs_model['model name']]
    marcs_model_content.append('  {:4.0f}.      Teff [K].         Last iteration; yyyymmdd={}'.format(marcs_model['teff'], marcs_model['last iteration']))
    marcs_model_content.append('  {:10.4E} Flux [erg/cm2/s]'.format(marcs_model['flux']))
    marcs_model_content.append('  {:10.4E} Surface gravity [cm/s2]'.format(marcs_model['g']))
    marcs_model_content.append('  {:<10.1f} Microturbulence parameter [km/s]'.format(marcs_model['vmicro']))
    marcs_model_content.append('  {:<10.1f} Mass [Msun]'.format(marcs_model['mass']))
    marcs_model_content.append(' {:+5.2f} {:+5.2f} Metallicity [Fe/H] and [alpha/Fe]'.format(marcs_model['[M/H]'], marcs_model['[alpha/Fe]']))
    if marcs_model['geometry'] == 'spherical':
        marcs_model_content.append('  {:10.4E} Radius [cm] at Tau(Rosseland)=1.0'.format(marcs_model['radius']))
        marcs_model_content.append('  {:10.5f} Luminosity [Lsun]'.format(marcs_model['luminosity']))
    elif marcs_model['geometry'] == 'plane-parallal':
        marcs_model_content.append('  {:10.4E} 1 cm radius for plane-parallel models'.format(marcs_model['radius']))
        marcs_model_content.append('  {:<10.1f} Luminosity [Lsun] FOR A RADIUS OF 1 cm!'.format(marcs_model['luminosity']))
    marcs_model_content.append('  {:4.2f} {:4.2f} {:5.3f} {:4.2f} are the convection parameters: alpha, nu, y and beta'.format(marcs_model['conv:alpha'], marcs_model['conv:nu'], marcs_model['conv:y'], marcs_model['conv:beta']))
    marcs_model_content.append('  {:7.5f} {:7.5f} {:8.2E} are X, Y and Z, 12C/13C={:.0f}'.format(marcs_model['X'], marcs_model['Y'], marcs_model['Z'], marcs_model['12C/13C']))
    marcs_model_content.append('Logarithmic chemical number abundances, H always 12.00')
    marcs_model_content.append(' ' + np.array2string(marcs_model['log abundance'], formatter={'float_kind':'{:6.2f}'.format}, max_line_width=74)[1:-1])
    marcs_model_content.append('{:4.0f} Number of depth points'.format(marcs_model['number of depth points']))
    marcs_model_content.append('Model structure')

    ms_1_formatters = ['{:3.0f}'.format, '{:5.2f}'.format, '{:7.4f}'.format, '{:10.3E}'.format, '{:>7.1f}'.format, 
                       '{:11.4E}'.format, '{:11.4E}'.format, '{:11.4E}'.format, '{:11.4E}'.format]
    marcs_model_content.append(marcs_model['model structure'][marcs_model['model structure'].columns[:9]].to_string(formatters=ms_1_formatters, index=False))
    ms_2_formatters = ['{:3.0f}'.format, '{:4.2f}'.format, '{:11.4E}'.format, '{:11.4E}'.format, '{:5.3f}'.format, 
                       '{:10.3E}'.format, '{:7.5f}'.format, '{:13.6E}'.format]
    marcs_model_content.append(marcs_model['model structure'][marcs_model['model structure'].columns[[0, 1, 9, 10, 11, 12, 13, 14]]].to_string(formatters=ms_2_formatters, index=False))
    marcs_model_content.append('Assorted logarithmic partial pressures')
    pp_1_formatters = ['{:3.0f}'.format, '{:6.3f}'.format, '{:6.2f}'.format, '{:6.2f}'.format, '{:6.2f}'.format, '{:6.2f}'.format, 
                       '{:6.2f}'.format, '{:6.2f}'.format, '{:6.2f}'.format, '{:6.2f}'.format, '{:6.2f}'.format, '{:6.2f}'.format]
    marcs_model_content.append(marcs_model['logarithmic partial pressures'][marcs_model['logarithmic partial pressures'].columns[:12]].to_string(formatters=pp_1_formatters, index=False))
    pp_2_formatters = ['{:3.0f}'.format, '{:6.2f}'.format, '{:6.2f}'.format, '{:6.2f}'.format, '{:6.2f}'.format, '{:6.2f}'.format, 
                       '{:6.2f}'.format, '{:6.2f}'.format, '{:6.2f}'.format, '{:6.2f}'.format, '{:6.2f}'.format, '{:6.2f}'.format]
    marcs_model_content.append(marcs_model['logarithmic partial pressures'][marcs_model['logarithmic partial pressures'].columns[[0] + list(range(12, 23))]].to_string(formatters=pp_2_formatters, index=False))
    marcs_model_content.append(marcs_model['logarithmic partial pressures'][marcs_model['logarithmic partial pressures'].columns[[0] + list(range(23, 34))]].to_string(formatters=pp_2_formatters, index=False))

    marcs_model_content = [i+'\n' for i in marcs_model_content]
    
    with open(save_name, 'w') as file:
        file.writelines(marcs_model_content)
    
    pass

def marcs2moog(marcs_model, save_name, abun_change=None, molecules_include=None):
    '''
    Convert MARCS format model to MOOG format.
    
    Parameters
    ----------
    marcs_model : dict
        The dict of MARCS model to be converted.
    save_name : str
        The name to save the model.
    abun_change : dict of pairs {int:float, ...}
        Abundance change, have to be a dict of pairs of atomic number and [X/Fe] values.
    molecules_include : list
        Molecules to be included to molecular calculation. Follows the MOOG notation.
    '''
    
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
    
    moog_model_content = ['KURUCZ']
    moog_model_content.append('TEFF={:.0f}, LOGG={:.1f}, [M/H]={:.2f}, [alpha/Fe]={:.2f}, Vmicro={:.1f}, mass={:.2f}'.format(marcs_model['teff'], np.log10(marcs_model['g']), marcs_model['[M/H]'], marcs_model['[alpha/Fe]'], marcs_model['vmicro'], marcs_model['mass']))
    moog_model_content.append('ntau=       ' + str(marcs_model['number of depth points']))
    ms_formatters = ['{:15.8E}'.format, '{:8.1f}'.format, '{:10.3E}'.format, '{:10.3E}'.format, '{:10.3E}'.format]
    moog_model_content.append(marcs_model['model structure'][['RHOX', 'T', 'Pg', 'Pe', 'KappaRoss']].to_string(formatters=ms_formatters, index=False, header=False))
    moog_model_content.append('    {:8.2E}'.format(marcs_model['vmicro']))

    # Element shift part
    if abun_change != None:
        abun_change_num = len(abun_change)
        moog_model_content.append('NATOMS      {}   {}'.format(abun_change_num, marcs_model['[M/H]']))
        ele_names = list(abun_change.keys())
        ele_names.sort()
        for element in ele_names:
            if element == 2:
                moog_model_content.append('      {:.2f}    {:.2f}'.format(element, xabu[element-1]+abun_change[element]))
            else:
                moog_model_content.append('      {:.2f}    {:.2f}'.format(element, xabu[element-1]+abun_change[element]+float(marcs_model['[M/H]'])))
    else:
        moog_model_content.append('NATOMS      0   {}'.format(marcs_model['[M/H]']))

    # Molecular line switches
    if molecules_include != None:
        molecules_num = len(molecules_include)
        moog_model_content.append('NMOL        {}'.format(molecules_num))
        molecules_str = ['{:11.1f}'.format(i) for i in molecules_include]
        molecules_str = ''.join(molecules_str)
        moog_model_content.append(molecules_str)
    else:
        moog_model_content.append('NMOL        0')
    moog_model_content = [i+'\n' for i in moog_model_content]
    
    with open(save_name, 'w') as file:
        file.writelines(moog_model_content)
        
    pass

def z2ao(z, out_type='grid'):
    '''
    Get the corresponding alpha (a) and O (o) abundance from z. 
    Applicable to MARCS Standard composition sphreical and plane-parallel models.
    Warning: sanity check for z is not performed here.
    
    Parameters
    ----------
    z : float
        Input overall metallicity.
    out_type : str, default 'grid'
        Output type, if 'grid' then only the z values in MARCS grid is allowed, and if 'other' all z values are allowed
    '''
    
    if out_type == 'grid':
        if z >= 0:
            return [0, 0]
        elif z == -0.25:
            return [0.1, 0.1]
        elif z == -0.50:
            return [0.2, 0.2]
        elif z == -0.75:
            return [0.3, 0.3]
        elif z <= -1:
            return [0.4, 0.4]
    else:
        fit_res = np.polyfit([-1, 0], [0.4, 0], 1)
        if z >= 0:
            return [0, 0]
        elif z <= -1:
            return [0.4, 0.4]
        else:
            return [np.polyval(fit_res, z), np.polyval(fit_res, z)]
    
def marcs_filename2para(file_name, out_type='list'):
    '''
    Convert marcs file name to grid point list.
    
    Parameters
    ----------
    file_name : str
        Name of the MARCS model file.
    out_type : str, default 'list'
        Output type, either 'list' or 'dict.' 
        
    Returns
    ----------
    marcs_point : list or dict
        The grid point values of MARCS model.
    '''
    
    marcs_name_single = file_name[:-4].split('_')
    marcs_point = {}
    marcs_point['chem'] = marcs_name_single[4]
    marcs_point['geo'] = marcs_name_single[0][0]
    marcs_point['teff'] = int(marcs_name_single[0][1:])
    _ = marcs_name_single.pop(4)
    _ = marcs_name_single.pop(0)
    for ele in marcs_name_single:
        marcs_point[ele[0]] = float(ele[1:])

    if out_type == 'dict':
        return marcs_point
    elif out_type == 'list':
        return list(marcs_point.values())
    
def para2marcs_filename(para_list):
    '''
    Convert parameters to MARCS file name.
    
    Parameters
    ----------
    para_list : list
        The list values of MARCS model. 
        
    Returns
    ----------
    file_name : str
        Name of the MARCS model file.
    '''
    return '{1}{2:.0f}_g{3:+4.1f}_m{4:3.1f}_t{5:02.0f}_{0}_z{6:+5.2f}_a{7:+5.2f}_c{8:+5.2f}_n{9:+5.2f}_o{10:+5.2f}_r{11:+5.2f}_s{12:+5.2f}.mod'.format(*para_list)

def interpolate_marcs_model(teff, logg, m_h, vmicro=2, mass=1, chem='st', geo='auto'):
    '''
    Interpolate the MARCS model.
    
    Parameters
    ----------
    teff : int
        The effective temperature of the model
    logg : float
        logg value of the model
    m_h : float
        [M/H] value (overall metallicity) of the model
    vmicro : float
        The microturbulance velocity of the model.
    mass : float
        The mass of the model. Only used in the spherical (geo='s').
    chem : str
        The chemical composition of the model. Either 'st', 'mc', 'hc', 'ap', 'ae', 'an', 'gs'.
    geo : str
        The geometry of the model. Either 's' for spherical, 'p' for plane-parallel or 'auto'. 
        If auto then the code will determine the geometry automatically according to stellar parameters.
    '''

    g = 10**logg
    # Determine the geo
    if geo == 'auto':
        if chem != 'st':
            raise ValueError("geo cannot be 'auto' when chem is not 'st'.")
        else:
            if logg >= 3.25:
                geo = 'p'
            else:
                geo = 's'

    if geo == 's':
        p = np.array([teff, g, m_h, vmicro, mass])
    else:
        p = np.array([teff, g, m_h, vmicro])
        
    marcs_grid = pd.read_csv(MOOG_file_path + '/pymoog_lf/model/marcs/{}/{}/grid_points.csv'.format(chem, geo))
    marcs_tri = pickle.load(open(MOOG_file_path + '/pymoog_lf/model/marcs/{}/{}/tri.pkl'.format(chem, geo), 'rb'))

    a, o = z2ao(p[2], out_type='other')

    tri_simplex = marcs_tri.find_simplex(p)
    if tri_simplex == -1:
            raise ValueError('The given stellar parameters are outside grid points, failed to interpolate.')
    else:
        tri_index = marcs_tri.simplices[tri_simplex]

    marcs_grid_sub = marcs_grid.loc[tri_index].reset_index(drop=True)

    if geo == 's':
        b = marcs_tri.transform[tri_simplex][:5].dot(np.transpose(p - marcs_tri.transform[tri_simplex][5]))
    elif geo == 'p':
        b = marcs_tri.transform[tri_simplex][:4].dot(np.transpose(p - marcs_tri.transform[tri_simplex][4]))
    b = np.concatenate([b, [1-sum(b)]])

    marcs_grid_use = marcs_grid_sub[b != 0].reset_index(drop=True)

    # Judge if the grid space is too large for interpolation
    teff_space_bad = np.ptp(marcs_grid_use['teff']) > 1500
    logg_space_bad = np.ptp(marcs_grid_use['logg']) > 0.5
    m_h_space_bad = np.ptp(marcs_grid_use['[M/H]']) > 0.5

    if np.any([teff_space_bad, logg_space_bad, m_h_space_bad]):
        raise ValueError('The separation between grid points is too large, failed to interpolate.')

    b = b[b != 0]

    if len(marcs_grid_use) == 1:
        # No interpolation
        marcs_model_interpolated = read_marcs_model(MOOG_file_path + '/pymoog_lf/model/marcs/{13}/{14}/{1}{2:.0f}_g{3:+4.1f}_m{4:3.1f}_t{5:02.0f}_{0:}_z{6:+5.2f}_a{7:+5.2f}_c{8:+5.2f}_n{9:+5.2f}_o{10:+5.2f}_r{11:+5.2f}_s{12:+5.2f}.mod'.format(*np.array(marcs_grid_use.loc[0, marcs_grid_use.columns[:-1]]), chem, geo))
    else:
        marcs_model_interpolated = {}
        # Interpolation
        marcs_model_interpolated['model name'] = para2marcs_filename([chem, geo, teff, logg, mass, vmicro, m_h, a, 0, 0, o, 0, 0])
        if geo == 's':
            marcs_model_interpolated['geometry'] = 'spherical'
        elif geo == 'p':
            marcs_model_interpolated['geometry'] = 'plane-parallel'
        for i in range(len(marcs_grid_use)):
            marcs_model_single = read_marcs_model(MOOG_file_path + '/pymoog_lf/model/marcs/{13}/{14}/{1}{2:.0f}_g{3:+4.1f}_m{4:3.1f}_t{5:02.0f}_{0:}_z{6:+5.2f}_a{7:+5.2f}_c{8:+5.2f}_n{9:+5.2f}_o{10:+5.2f}_r{11:+5.2f}_s{12:+5.2f}.mod'.format(*np.array(marcs_grid_use.loc[i, marcs_grid_use.columns[:-1]]), chem, geo))
            if i == 0:
                for name in ['teff', 'last iteration', 'flux', 'g', 'vmicro', 'mass', '[M/H]', '[alpha/Fe]', 'radius', 'luminosity', 
                             'conv:alpha', 'conv:nu', 'conv:y', 'conv:beta', 'X', 'Y', 'Z', '12C/13C',
                             'log abundance', 'model structure', 'logarithmic partial pressures']:
                    if name == 'last iteration':
                        marcs_model_interpolated[name] = marcs_model_single[name]
                    else:
                        marcs_model_interpolated[name] = marcs_model_single[name] * b[i]
                N_depth = marcs_model_single['number of depth points']
                marcs_model_interpolated['number of depth points'] = N_depth
            else:
                for name in ['teff','flux', 'g', 'vmicro', 'mass', '[M/H]', '[alpha/Fe]', 'radius', 'luminosity', 
                             'conv:alpha', 'conv:nu', 'conv:y', 'conv:beta', 'X', 'Y', 'Z', '12C/13C',
                             'log abundance', 'model structure', 'logarithmic partial pressures']:
                    marcs_model_interpolated[name] += marcs_model_single[name] * b[i]
                if marcs_model_single['number of depth points'] != N_depth:
                    raise ValueError('N-depth not consistent')

    return marcs_model_interpolated

def interpolate_model(teff, logg, m_h, vmicro=2, mass=1, abun_change=None, molecules_include=None, save_name=None, model_type='marcs', chem='st', geo='auto'):
    '''
    Interpolate Kurucz model or MARCS model and convert to moog format.
    '''

    if model_type == 'kurucz':
        interpolate_kurucz_model(teff, logg, m_h, vmicro=vmicro, abun_change=abun_change, molecules_include=molecules_include, to_path=save_name)
    elif model_type == 'marcs':
        marcs_model_interpolated = interpolate_marcs_model(teff, logg, m_h, vmicro=vmicro, mass=mass, chem=chem, geo=geo)
        marcs2moog(marcs_model_interpolated, save_name, abun_change=abun_change, molecules_include=molecules_include)

    pass