from PyAstronomy import pyasl
import math
model = pyasl.getKuruczModel(4250, 4.5, 0.1)

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

    model = pyasl.getKuruczModel(teff,logg,m_h)

    outf = open(file_name,"w")

    for line in model:
        outf.write(line)
        outf.write("\n")

    outf.close()


KURUCZ_APOGEE_download(4250, 4.5, 0.1);
