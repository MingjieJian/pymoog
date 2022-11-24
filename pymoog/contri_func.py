#!/usr/bin/python
from pymoog import private
from pymoog import synth
from pymoog import line_data

# The functions to calculate the contribution functions.
# Note that the atmosphere and lines have to be set to 2 and 4, and only one line each time.
# s.rundir_path = '{}/.pymoog/rundir/'.format(private.os.environ['HOME'])

def kappa_l_single(line_opacity_content):
    '''
    Function to convert list of kappa_l output for one line into np.array. 
    '''
    
    line_index = int(line_opacity_content[1].split()[1])
    wavelength = float(line_opacity_content[1].split()[4])

    kappa = private.np.array(''.join(line_opacity_content[2:]).split(), dtype=float) 
    return line_index, wavelength, kappa

def extract_kappa_l(file_name, line_id='all'):
    # Extract the line opacities

    # Open the file
    file = open(file_name, 'r')
    out1 = file.readlines()

    # Extract the "LINE OPACITIES" part
    for i in range(len(out1)):
        if 'LINE OPACITIES' in out1[i]:
            start_index = i
        if 'SPECTRUM DEPTHS' in out1[i]:
            end_index = i
    line_opacity_content = out1[start_index:end_index]

    sep_index = []
    for i in range(len(line_opacity_content)):
        if line_opacity_content[i] == '\n':
            sep_index.append(i)

    kappa_l = {}
    for i in range(len(sep_index)-1):
        line_index, wavelength, kappa = kappa_l_single(line_opacity_content[sep_index[i]:sep_index[i+1]])
        kappa_l[line_index] = [wavelength, kappa]

    if line_id == 'all':
        return kappa_l
    elif type(line_id) == int:
        if line_id >= 1 and line_id <= len(kappa_l):
            return kappa_l[line_id]
        else:
            raise ValueError('line_id cannot be < 1 or >= the length of line list.')
    else:
        raise ValueError('line_id is not "all" or an integer.')

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
    kappa_c = private.pd.read_table('temp', sep=' +', engine='python')
    private.os.remove('temp')

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

    atmosphere = private.pd.read_table('temp', sep=' +', engine='python')
    private.os.remove('temp')

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
    kappa_ref = private.np.array(record_str.split(), dtype='float')

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
        z = private.np.array(z)
    dz = [0]
    for i in range(len(z)-1):
        dz.append(z[i+1] - z[i])
    return private.np.array(dz)

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
        tau_c_i = private.np.sum(kappa_c[:i+1] / kappa_ref[:i+1] * dtau_ref[:i+1])
        tau_c.append(tau_c_i)
        tau_l_i = private.np.sum(kappa_l[:i+1] / kappa_ref[:i+1] * dtau_ref[:i+1])
        tau_l.append(tau_l_i)
    tau_c = private.np.array(tau_c)
    tau_l = private.np.array(tau_l)

    return tau_c, tau_l

def cal_contri_func_file(file_name, line_id=1):
    wavelength, kappa_l = extract_kappa_l(file_name, line_id=line_id)
    kappa_c = extract_kappa_c(file_name)
    kappa_ref = extract_kappa_ref(file_name)
    atmosphere = extract_atmosphere(file_name)

    tau_ref = atmosphere['tauref'].values
    tau_c, tau_l = calculate_tau(tau_ref, 10**kappa_c['kaplam'].values, kappa_l, kappa_ref)

    # Calculate source function from temperature
    # S = private.BlackBody(wavelength, kappa_c['T(i)']).value
    S = private.BlackBody(temperature=kappa_c['T(i)']*private.u.K)(wavelength*private.u.AA).cgs.value
    kappa_c = 10**kappa_c['kaplam'].values
    CF_I_c = S * private.np.exp(-tau_c) * atmosphere['tauref'] * private.np.log(10) * kappa_c / kappa_ref
    dlog_tau_ref = calculate_dz(private.np.log10(tau_ref))
    Ic_tau_ref = []
    for i in range(len(dlog_tau_ref)-1):
        Ic_tau_ref.append(private.np.sum(CF_I_c[i:] * dlog_tau_ref[i:]))
    Ic_tau_ref.append(0)
    Ic_tau_ref = private.np.array(Ic_tau_ref)
    CF_Dlp = Ic_tau_ref * private.np.exp(-tau_l) * tau_ref * private.np.log(10) * kappa_l / kappa_ref
    CF_Ilp = S * private.np.exp(-tau_l-tau_c) * private.np.log(10) * tau_ref * kappa_l / kappa_ref
    CF_Dl = CF_Dlp - CF_Ilp
    CF_Ilc = S * private.np.exp(-tau_l-tau_c) * private.np.log(10) * tau_ref * kappa_c / kappa_ref
    CF_Il = (1 + kappa_l/kappa_c) * S * private.np.exp(-tau_l-tau_c) * private.np.log(10) * tau_ref * kappa_c / kappa_ref

    CF_dict = {'tau_ref':tau_ref, 'CF_Ic':CF_I_c.values, 'CF_Dlp':CF_Dlp, 'CF_Ilp':CF_Ilp, 'CF_Dl':CF_Dl,
               'CF_Ilc':CF_Ilc, 'CF_Il':CF_Il}
    return CF_dict, atmosphere

def cal_blending_ratio(teff, logg, fe_h, resolution, line_list, wav_range, weedout=True):
    
    wav_start = wav_range[0] - 7
    wav_end = wav_range[1] + 7
    s = synth.synth(teff, logg, fe_h, wav_start, wav_end, resolution, line_list=line_list, weedout=weedout)

    # Whole spectra
    s.prepare_file()
    s.run_moog()
    s.read_spectra(remove=False)
    s.weedout = False
    wav_all, flux_all = s.wav, s.flux

    # Target line excluded
    linelist_all = line_data.read_linelist(s.rundir_path + line_list)
    indices = (linelist_all['wavelength']>=wav_range[0]) & (linelist_all['wavelength']<=wav_range[1])
    linelist_out = linelist_all[indices].reset_index(drop=True)
    line_index_all = linelist_all[indices].index

    r_blend_depth = []
    r_blend_EW = []
    depth = []
    for line_index in line_index_all:

        linelist_exclude = linelist_all.drop(line_index).reset_index(drop=True)
        line_data.save_linelist(linelist_exclude, s.rundir_path + line_list)
        s.run_moog()
        s.read_spectra(remove=False)
        wav_exclude, flux_exclude = s.wav, s.flux

        # Calculate the EW and blending fraction
        linelist_target = linelist_all.loc[line_index:line_index].reset_index(drop=True)
        line_wavlength = linelist_target.loc[0, 'wavelength']
        EW = (private.np.sum(1-flux_all)*0.02 - private.np.sum(1-flux_exclude)*0.02) * 1000
        if flux_exclude[private.np.argmin(private.np.abs(wav_exclude-line_wavlength))] - flux_all[private.np.argmin(private.np.abs(wav_exclude-line_wavlength))] < 1e-5:
            r_blend_depth.append(1.0)
        else:
            r_blend_depth.append((1-flux_exclude[private.np.argmin(private.np.abs(wav_exclude-line_wavlength))]) / (1-flux_all[private.np.argmin(private.np.abs(wav_all-line_wavlength))]))
            
        EW_indices = private.np.abs(wav_exclude-line_wavlength) < line_wavlength / resolution *3
        r_blend_EW.append(private.np.sum(1-flux_exclude[EW_indices]) / private.np.sum(1-flux_all[EW_indices]))
        depth.append(1 - private.np.min(flux_all[private.np.abs(wav_all-line_wavlength) <= 0.03]))
    linelist_out['depth'] = depth
    linelist_out['r_blend_depth'] = r_blend_depth
    linelist_out['r_blend_EW'] = r_blend_EW

    return linelist_out

def plot_contri_func(teff, logg, fe_h, resolution, line_list, line_wav_input=None, line_id=None, target_line_df=None, smooth_para=None, plot=True, dpi=100, plot_Dlp=True, plot_Dl=True, plot_Ilp=True, plot_Ic=True):

    
    if target_line_df is None and (line_wav_input is None or line_id is None):
        raise ValueError('Please provide target_line_df or both line_wav_input and line_id.')
    
    if target_line_df is not None:
        wav_start = private.np.min(target_line_df['wavelength'])-7
        wav_end = private.np.max(target_line_df['wavelength'])+7
        s = synth.synth(teff, logg, fe_h, wav_start, wav_end, resolution, line_list=line_list)
    else:
        s = synth.synth(teff, logg, fe_h, line_wav_input-7, line_wav_input+7, resolution, line_list=line_list)

    # Whole spectra 
    if smooth_para is None:
        s.prepare_file()
    else:
        s.prepare_file(smooth_para=smooth_para)
    s.run_moog()
    s.read_spectra(remove=False)
    wav_all, flux_all = s.wav, s.flux

    # Target line excluded
    linelist_all = line_data.read_linelist(s.rundir_path + 'line.list')
    if target_line_df is not None:
        line_index_all = line_data.find_lines(target_line_df, linelist_all)
    else:
        smooth_width = line_wav_input / resolution
        indices = (private.np.abs(linelist_all['wavelength']-line_wav_input) <= smooth_width) & (linelist_all['id'] == line_id)
        line_index_all = linelist_all[indices].index

    if plot:
        private.plt.figure(figsize=(14, 5*len(line_index_all)), dpi=dpi)


    plot_index = 1
    CF_res = []
    for line_index in line_index_all:
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
        if smooth_para is None:
            s.prepare_file(atmosphere=2, lines=4)
        else:
            s.prepare_file(smooth_para=smooth_para, atmosphere=2, lines=4)
        line_data.save_linelist(linelist_target, s.rundir_path + 'line.list')
        s.run_moog()
        s.read_spectra(remove=False)
        wav_target, flux_target = s.wav, s.flux

        # Calculate the EW and blending fraction
        EW = (private.np.sum(1-flux_all)*0.02 - private.np.sum(1-flux_exclude)*0.02) * 1000
        try:
            depth = 1 - private.np.min(flux_all[private.np.abs(wav_all-line_wavlength) <= 0.1])
        except:
            # print(flux_all[private.np.abs(wav_all-line_wavlength) <= 0.1])
            depth = private.np.nan
        r_blend_depth = (1-flux_exclude[private.np.argmin(private.np.abs(wav_exclude-line_wavlength))]) / (1-flux_all[private.np.argmin(private.np.abs(wav_all-line_wavlength))])

        # Plot the line information
        CF_dict, atmosphere = cal_contri_func_file(s.rundir_path + 'MOOG.out1', line_id=1)
        s.read_spectra(remove=False)
        Ic_sum = private.np.sum(CF_dict['CF_Ic'])
        Dlp_sum = private.np.sum(CF_dict['CF_Dlp']) 

        if plot:
            ax = private.plt.subplot(len(line_index_all),2,plot_index)
            private.plt.plot(wav_all, flux_all, label='all lines included')
            private.plt.plot(wav_exclude, flux_exclude, label='target line excluded')
            private.plt.plot(wav_target, flux_target, label='target line only')
            private.plt.plot(line_wavlength, 1, alpha=0, label='EW={:.2f} m$\mathrm{{\AA}}$'.format(EW))
            private.plt.plot(line_wavlength, 1, alpha=0, label='$d$={:.2f}'.format(depth))
            private.plt.plot(line_wavlength, 1, alpha=0, label='$f_\mathrm{{blend,depth}}$={:.2f}'.format(r_blend_depth))
            private.plt.axvline(wav_all[private.np.argmin(private.np.abs(wav_all-line_wavlength))],zorder=0, ls='--', c='gray')
            private.plt.xlim(line_wavlength-5, line_wavlength+5)
            private.plt.legend()
            private.plt.xlabel(r'Wavelength ($\mathrm{\AA}$)')
            private.plt.ylabel('Normalized flux')
            if target_line_df is not None:
                line_id = linelist_all.loc[line_index, 'id']
            private.plt.title('Teff={:.0f}, logg={:.2f}, [Fe/H]={:.2f}, line_id={:.1f}, EP={:.2f}, loggf={:.3f}'.format(teff, logg, fe_h, line_id, line_EP, line_loggf))
            plot_index += 1

            ax = private.plt.subplot(len(line_index_all),2,plot_index)

            if plot_Dlp:
                private.plt.plot(private.np.log10(CF_dict['tau_ref']), CF_dict['CF_Dlp'] / max(CF_dict['CF_Ic']), label='CF: $D_l^p$')
            if plot_Ilp:
                private.plt.plot(private.np.log10(CF_dict['tau_ref']), CF_dict['CF_Ilp'] / max(CF_dict['CF_Ic']), label='CF: $I_l^p$')
            if plot_Dl:
                private.plt.plot(private.np.log10(CF_dict['tau_ref']), CF_dict['CF_Dl'] / max(CF_dict['CF_Ic']), label='CF: $D_l$')
            if Dlp_sum / Ic_sum < 0.1:
                private.plt.ylim(private.plt.ylim())
            if plot_Ic:
                private.plt.plot(private.np.log10(CF_dict['tau_ref']), CF_dict['CF_Ic'] / max(CF_dict['CF_Ic']), c='gray', label='CF: $I_c$')

            private.plt.xlabel(r'$\log{\tau_\mathrm{ref}}$')
            private.plt.ylabel('CF: max($I_c$) normalized as 1')
            private.plt.title('Contribution function')

            ax2 = private.plt.twinx()
            private.plt.plot(private.np.log10(atmosphere['tauref']), atmosphere['T'], ls='--', c='gray', zorder=0, label='T (K)')
            private.plt.ylabel('T (K)')
            lines, labels = ax.get_legend_handles_labels()
            lines2, labels2 = ax2.get_legend_handles_labels()
            ax2.legend(lines + lines2, labels + labels2)
            
            plot_index += 1
            private.plt.tight_layout()
        CF_res.append(CF_dict)
        
    return CF_res, atmosphere