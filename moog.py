from . import private
import math
import pexpect
import os
import sys
import re
import subprocess
import numpy as np

# The wrapper of MOOG.
#  To DO:
#  Function to convert KURUCZ atmosphere file to MOOG format
#  Function to create batch.par file
#  Function to read the output of MOOG

def v2fn(value_list):
    '''
    The function to translate the value into name, according to the grids of Meszaros (2016). The values input have to be the same as the grid.

    Parameters
    ----------
    value_list : list, [m_h, c_m, alpha_m, teff, logg]
        List containing the input parameters.

    Returns
    ----------
    fn_list : list
        List of strings to be put into the path of atmosphere model.
    '''
    fn_list = []
    for value in value_list[:-2]:
        if value < 100:
            fn = '{:02d}'.format(int(math.ceil(abs(value)*10)))
        else:
            fn = '{:4}'.format(value)
        if value < 0:
            fn = 'm'+fn
        elif value < 100:
            fn = 'p'+fn
        fn_list.append(fn)
    for value in value_list[-2:]:
        if value < 100:
            fn = '{:02d}'.format(int(abs(value)*10))
        else:
            fn = '{:4}'.format(value)
        fn_list.append(fn)
    return fn_list

def KURUCZ_APOGEE_download(teff, logg, m_h, c_m=0, alpha_m=0, to_path='./'):
    '''
    Download the ATLAS9-APOGEE grid file.

    Parameters
    ----------
    teff : int
        The effective temperature of the model
    logg : float
        logg value of the model
    m_h : float
        [M/H] value (overall metallicity) of the model
    c_m : float
        [C/M] value of the model
    alpha_m : float
        [alpha/M] value of the model
    to_path : str
        The path to save the model file.

    Returns
    ----------
    to_name : str
        The path and name of saved model.
    dl_situtation : boolen
        Download situation, True if succeeded, False if falied.

    '''
    if to_path[-1] != '/':
        to_path = to_path + '/'
    vlist = [m_h, c_m, alpha_m, teff, logg]
    vlist = v2fn(vlist) + [to_path]
    # Download the model file.
    file_command = 'scp public@133.11.16.65:/share/CACHEDEV1_DATA/Public/Mezaros2012-APOGEE-grid/{0}_all_mod/m{0}c{1}o{2}/am{0}c{1}o{2}t{3}g{4}v20.mod {5}'.format(*vlist)
    p = pexpect.spawn(file_command)
    p.expect("{}@{}'s password:".format('public', '133.11.16.65'))
    result = p.sendline('Public')
    p.wait()

    # Judge if the file is downloaded
    to_name = '{5}am{0}c{1}o{2}t{3}g{4}v20.mod'.format(*vlist)
    dl_situtation = os.path.isfile('{5}am{0}c{1}o{2}t{3}g{4}v20.mod'.format(*vlist))

    return to_name, dl_situtation

def KURUCZ_APOGEE_convert(file_path, vmicro=2.0, abun_change=None):
    '''
    Convert the model file of ATLAS9-APOGEE in to MOOG format.

    Parameters
    ----------
    file_path : str
        The path of donloaded model file.
    v_micro : float, default 2.0
        Microturbulance velocity of the spectra.
    abun_change : list of pairs [int, float]
        Abundance change, have to be a list of pairs of atomic number and [X/Fe] values.

    Returns
    ----------
    k_model_path : str
        The path of converted file.
    cv_situation : boolen
        Convert situation, True if file exist, False if falied.
    '''

    #print(file_path)
    model_file = open(file_path)
    # Convert the model files into MOOG format. Only for the model in the website above.

    # Read and save the first two lines (except 'TITLE ') into header.
    header = model_file.readline()[:-2] + model_file.readline()[6:-2].strip()
    m_h = re.findall(r'\[(.*)\]', header)[0]


    # Read the abundance change as well as model lines.
    temp = model_file.readline() + model_file.readline() + model_file.readline()

    abun_list = ''
    temp = model_file.readline()
    while temp[:17] == ' ABUNDANCE CHANGE':
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
    temp = temp.split()
    model_lines = []
    model_linen = int(temp[2])
    for i in range(int(temp[2])):
        model_lines.append(model_file.readline().split()[:7])

    # Prepare the microtrubulance value.
    vmicro = '{}E00'.format(vmicro)

    # Write the model file.
    # header, abun09, model_lines and vmicro
    k_model_path = file_path.rsplit('/',1)[0] + '/k' + file_path.rsplit('/',1)[1][1:]
    k_model_file = open(k_model_path, 'w')

    # Header part
    k_model_file.writelines('KURUCZ\n')
    k_model_file.writelines(header + '\n')

    # Model part
    k_model_file.writelines('ntau=       ' + str(model_linen) + '\n')
    for i in model_lines:
        k_model_file.writelines(' ' + ' '.join(i) + '\n')

    # Microturbulant velocity part
    k_model_file.writelines('    ' + vmicro + '\n')

    # Element shift part
    if abun_change != None:
        abun_change_num = len(abun_change)
        k_model_file.writelines('NATOMS      {}   {}\n'.format(abun_change_num, m_h))
        for abun in abun_change:
            k_model_file.writelines('      {:.2f}    {:.2f}\n'.format(abun[0], xabu[abun[0]-1]+abun[1]+float(m_h)))
    else:
        k_model_file.writelines('NATOMS      0   {}\n'.format(m_h))
#         for i in abun09:
#             k_model_file.writelines('          {:2g} {:.3f}\n'.format(i[0], i[1]))

    # Molecular line switches (temporary closed)
    k_model_file.writelines('NMOL        0')
    k_model_file.close()

    cv_situation = os.path.isfile(k_model_path)
    return k_model_path, cv_situation

def create_para_file(k_model_path, linelist_path, start_wav=15167.0, end_wav=16767.0, del_wav=0.02, smooth='g', smooth_width=0.75, atmosphere=1, lines=1):
    '''
    Function for creating the parameter file of batch.par
    '''
    MOOG_para_file = open('batch.par', 'w')
    # Parameter list of MOOG: standard output file (1), summary output file (2), smoothed output file (3),
    #                         begin wavelength, end wavelength, wavelength step;
    #                         smoothing function, Gaussian FWHM, vsini, limb darkening coefficient,
    #                         Macrotrubulent FWHM, Lorentzian FWHM
    smooth_para = [smooth, smooth_width, 0.0, 0.0, 0.0, 0.0]
    #MOOG_para_file = open('batch.par', 'w')
    MOOG_contant = ["synth\n",
                    "atmosphere         {}\n".format(atmosphere),
                    "lines              {}\n".format(lines),
                    "standard_out       '{}'\n".format('MOOG.out1'),
                    "summary_out        '{}'\n".format('MOOG.out2'),
                    "smoothed_out       '{}'\n".format('MOOG.out3'),
                    "model_in           '{}'\n".format(k_model_path),
                    "lines_in           '{}'\n".format(linelist_path),
                    "terminal           'x11'\n",
                    "synlimits\n",
                    "  {}  {}  {}  2.50\n".format(start_wav, end_wav, del_wav),
                    "plot        3\n",
                    "plotpars    1\n",
                    "  0.0  0.0  0.0  0.0 \n",
                    "  0.0  0.0  0.0  0.0 \n",
                    "  '{}'  {}  {}  {}  {}  {}\n".format(*smooth_para)
                   ]
    MOOG_para_file.writelines(MOOG_contant)
    MOOG_para_file.close()

def moog_run(output=False, version='FEB2017'):
    '''
    Run MOOG and print the reuslt if required.

    Parameters
    ----------
    output : boolen, default False
        If set to True, then print the out put of MOOG.

    Returns
    ----------
    None. Three files MOOG.out1, MOOG.out2 and MOOG.out3 will be save in the current working path.
    '''
    if version == 'FEB2017':
        MOOG_run = subprocess.run(['/home/mingjie/software/MOOG/FEB2017_modified/MOOGSILENT'], stdout=subprocess.PIPE)
    elif version == 'NOV2019':
        MOOG_run = subprocess.run(['/home/mingjie/software/MOOG/NOV2019_modified/MOOGSILENT'], stdout=subprocess.PIPE)
    if output:
        MOOG_run = str(MOOG_run.stdout, encoding = "utf-8").split('\n')
        MOOG_output = []
        for i in MOOG_run:
            if len(i) > 12:
                ansi_escape = re.compile('\x1b\[...H')
                temp = ansi_escape.sub('', i)
                ansi_escape = re.compile('\x1b\[....H')
                temp = ansi_escape.sub('', temp)
                ansi_escape = re.compile('\x1b\[H')
                temp = ansi_escape.sub('', temp)
                ansi_escape = re.compile('\x1b\[2J')
                MOOG_output.append(ansi_escape.sub('', temp))
        for i in MOOG_output:
           print(i)

def read_spectra(type='smooth'):
    '''
    Read the output spectra of MOOG.

    Parameters
    ----------
    type : str, default 'smooth'
        Decide the type of spectra to be read. 'smooth' will read the smoothed spectra (MOOG.out3), and 'standard' will read the un-smoothed one (MOOG.out2).

    Returns
    ---------
    wav : a numpy array
        An array of wavelength
    flux : a numpy array
        An array of flux
    '''
    if type == 'standard':
        models_file = open('MOOG.out2')
        models = models_file.readline()
        models = models_file.readline()
        models = models_file.read().split()
        models = [float(i) for i in models]
        x = range(round((models[1]-models[0])/models[2])+1)
        models_wav = []
        for i in x:
            models_wav.append(models[0] + models[2]*i)
        model_flux = np.array(models[4:])
        return np.array(model_wav), np.array(model_flux)
    elif type == 'smooth':
        models_file = open('MOOG.out3')
        models = models_file.readline()
        models = models_file.readline()
        models = models_file.readlines()
        wavelength = []
        depth = []
        for i in models:
            temp = i.split()
            wavelength.append(float(temp[0]))
            depth.append(float(temp[1]))
        return np.array(wavelength), np.array(depth)


def read_vald_linelist():
    '''
    Read the VALD linelist in YJ band created by Taniguchi.
    '''

    vald_yj_linelist = private.pd.read_fwf(private.package_path + '/file/vald_yj',
            colspecs=[(0,11), (11,21), (21,31), (31,41), (41,51), (51,61)],
            names=['wavelength', 'id', 'EP', 'loggf', 'C6', 'D0'],
            skiprows=1)
    return vald_yj_linelist
