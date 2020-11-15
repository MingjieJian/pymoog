#!/usr/bin/python

def extract_kappa_l(file_name, line_id='all'):
    # Extract the line opacities

    # Open the file
    file = open(file_name)
    out1 = file.readlines()

    # Extract the "LINE OPACITIES" part
    for i in range(len(out1)):
        if 'LINE OPACITIES' in out1[i]:
            start_index = i
        if 'SPECTRUM DEPTHS' in out1[i]:
            end_index = i
    line_opacity_content = out1[start_index:end_index]

    # Convert the line opacities into arrays
    line_opacity_content = line_opacity_content[2:]

    record_switch = 0
    line_opacity_array = []
    
    for line in line_opacity_content:
        if record_switch == 0:
            record_str = ''
        elif record_switch == 1:
            record_str = record_str + line
        if 'LINE' in line:
            line_wav = line.split(' ')[-1]
            record_switch = 1
        elif line == '\n':
            record_switch = 0
            line_opacity_array = line_opacity_array + [np.array(record_str.split(), dtype='float')]

    if line_id == 'all':
        return line_opacity_array
    elif line_id > 0:
        return line_opacity_array[line_id-1]            
    else:
        raise ValueError('line_id must larger than 0.')
        
def extract_kappa_c(file_name):
    # Read the continuum opacity

    # Open the file
    file = open(file_name)
    out1 = file.readlines()

    # Extract the "LINE OPACITIES" part
    for i in range(len(out1)):
        if 'log opacities due to H, He, e-' in out1[i]:
            start_index = i
        if 'log opacities due to "metal" edges' in out1[i]:
            end_index = i
    con_opacity_content = out1[start_index+1:end_index]

    file_temp = open('temp', 'w')
    file_temp.writelines(con_opacity_content)
    file_temp.close()
    kappa_c = pd.read_table('temp', sep=' +', engine='python')
    os.remove('temp')

    return kappa_c

def extract_atmosphere(file_name):
    # Open the file
    file = open(file_name)
    out1 = file.readlines()

    # Extract the "LINE OPACITIES" part
    for i in range(len(out1)):
        if 'INPUT ATMOSPHERE QUANTITIES' in out1[i]:
            start_index = i
        if 'INPUT ABUNDANCES: (log10 number' in out1[i]:
            end_index = i
    atmosphere_content = out1[start_index+1:end_index]
    for i in range(len(atmosphere_content)):
        atmosphere_content[i] = atmosphere_content[i].replace('D', 'E')
    file_temp = open('temp', 'w')
    file_temp.writelines(atmosphere_content)
    file_temp.close()

    atmosphere = pd.read_table('temp', sep=' +', engine='python')
    os.remove('temp')

    return atmosphere

def extract_kappa_ref(file_name):
    # Open the file
    file = open(file_name)
    out1 = file.readlines()

    # Extract the "LINE OPACITIES" part
    for i in range(len(out1)):
        if 'KAPREF ARRAY' in out1[i]:
            start_index = i
        if 'For these computations, some abu' in out1[i]:
            end_index = i
    line_opacity_content = out1[start_index+1:end_index]

    record_str = ''
    for line in line_opacity_content:
        record_str = record_str + line
    record_str = record_str.replace('D', 'E')
    kappa_ref = np.array(record_str.split(), dtype='float')

    return kappa_ref

def calculate_dz(z):
    '''
    A function to calculate dz.

    Parameters
    ----------
    z : np.array
        The geometrical depth of atmosphere grid.

    Returns
    -------
    dz : np.array
        An array of dz, with the first element as 0.
    '''

    if type(z) != 'numpy.ndarray':
        z = np.array(z)
    dz = [0]
    for i in range(len(z)-1):
        dz.append(z[i+1] - z[i])
    return np.array(dz)

def calculate_tau(tau_ref, kappa_c, kappa_l, kappa_ref):
    '''
    A function to calculate tau_c and tau_l.

    Parameters
    ----------
    z : numpy array
        The geometrical depth of atmosphere grid. The unit should be cm.

    kappa_df : pandas dataframe
        Output of the PANDORA "LINE (X/Y)" kappa table.

    Returns
    -------
    tau_c : numpy array
        An array of tau_c
    tau_l : numpy array
        An array of tau_l
    '''
    dtau_ref = calculate_dz(tau_ref)

    tau_c = []
    tau_l = []
    for i in range(len(dtau_ref)):
        tau_c_i = np.sum(kappa_c[:i+1] / kappa_ref[:i+1] * dtau_ref[:i+1])
        tau_c.append(tau_c_i)
        tau_l_i = np.sum(kappa_l[:i+1] / kappa_ref[:i+1] * dtau_ref[:i+1])
        tau_l.append(tau_l_i)
    tau_c = np.array(tau_c)
    tau_l = np.array(tau_l)

    return tau_c, tau_l

def cal_contri_func(file_name):
    kappa_l = extract_kappa_l(file_name, line_id=1)
    kappa_c = extract_kappa_c(file_name)
    kappa_ref = extract_kappa_ref(file_name)
    atmosphere = extract_atmosphere(file_name)

    tau_ref = atmosphere['tauref'].values
    tau_c, tau_l = calculate_tau(tau_ref, 10**kappa_c['kaplam'].values, kappa_l, kappa_ref)

    # Calculate source function from temperature
    S = blackbody_lambda(15207.53, kappa_c['T(i)']).value
    kappa_c = 10**kappa_c['kaplam'].values
    CF_I_c = S * np.exp(-tau_c) * atmosphere['tauref'] * np.log(10) * kappa_c / kappa_ref
    dlog_tau_ref = pandora.calculate_dz(np.log10(tau_ref))
    Ic_tau_ref = []
    for i in range(len(dlog_tau_ref)-1):
        Ic_tau_ref.append(np.sum(CF_I_c[i:] * dlog_tau_ref[i:]))
    Ic_tau_ref.append(0)
    Ic_tau_ref = np.array(Ic_tau_ref)
    CF_Dlp = Ic_tau_ref * np.exp(-tau_l) * tau_ref * np.log(10) * kappa_l / kappa_ref
    CF_Ilp = S * np.exp(-tau_l-tau_c) * np.log(10) * tau_ref * kappa_l / kappa_ref
    CF_Dl = CF_Dlp - CF_Ilp
    CF_Ilc = S * np.exp(-tau_l-tau_c) * np.log(10) * tau_ref * kappa_c / kappa_ref
    CF_Il = (1 + kappa_l/kappa_c) * S * np.exp(-tau_l-tau_c) * np.log(10) * tau_ref * kappa_c / kappa_ref

    CF_dict = {'tau_ref':tau_ref, 'CF_I_c':CF_I_c, 'CF_Dlp':CF_Dlp, 'CF_Ilp':CF_Ilp, 'CF_Dl':CF_Dl,
               'CF_Ilc':CF_Ilc, 'CF_Il':CF_Il}
    return CF_dict