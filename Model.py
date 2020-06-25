from PyAstronomy import pyasl
import math
import re
import numpy as np
import os

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
    vlist = [m_h, c_m, alpha_m, teff, logg]
    vlist = v2fn(vlist) + [to_path]
    file_name = '{5}am{0}c{1}o{2}t{3}g{4}v20.mod'.format(*vlist)

    #The model is invoked here
    model = pyasl.getKuruczModel(teff,logg,m_h)

    outf = open(file_name,"w")

    for line in model:
        outf.write(line)
        outf.write("\n")

    outf.close()

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


if __name__=="__main__": 
    KURUCZ_APOGEE_download(4250, 4.5, 0.1,to_path='./demo/');
    KURUCZ_APOGEE_convert("./demo/amp01cp00op00t4250g45v20.mod")