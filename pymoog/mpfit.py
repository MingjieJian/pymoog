from . import private
from . import synth
from . import line_data
from . import weedout
import spectres

def cal_partial_f_v(teff, logg, m_h, wav_start, wav_end, fwhm_broad, vmicro_in, rv_in, wav_in, mode, abun_change=None, diff_v=0.01, line_list='vald_winered', prefix=''):
    '''
    Here the resolution in synth (20000) is a dummy parameter.
    mode must be either 'vbroad' or 'vmicro'.
    '''
    if mode not in ['vmicro', 'vbroad']:
        raise ValueError("mode must be either 'vbroad' or 'vmicro'.")
    
    diff_v_dict = {'vmicro':0, 'vbroad':0}
    diff_v_dict[mode] = diff_v
    
    if diff_v <= 0:
        raise ValueError("difff_v must > 0.")
    
    # convert rv to wavelength
    del_wav = private.np.abs(rv_in / 3e5 * private.np.mean(wav_in))

    s = synth.synth(teff, logg, m_h, wav_start-0.75-del_wav, wav_end+0.75+del_wav, 20000, line_list=line_list, weedout=True, prefix=prefix, vmicro=vmicro_in+diff_v_dict['vmicro'])
    s.prepare_file(smooth_para=['g', fwhm_broad+diff_v_dict['vbroad'], 0, 0, 0, 0],
                   abun_change=abun_change)
    s.run_moog()
    s.read_spectra()
    s.wav = s.wav * (1 + rv_in/3e5)
    flux_p = spectres.spectres(wav_in, s.wav, s.flux)

    s_ = synth.synth(teff, logg, m_h, wav_start-0.75-del_wav, wav_end+0.75+del_wav, 20000, line_list=line_list, weedout=True, prefix=prefix, vmicro=vmicro_in-diff_v_dict['vmicro'])
    s_.prepare_file(smooth_para=['g', fwhm_broad-diff_v_dict['vbroad'], 0, 0, 0, 0],
                    abun_change=abun_change)
    s_.run_moog()
    s_.read_spectra()
    s_.wav = s_.wav * (1 + rv_in/3e5)
    flux_m = spectres.spectres(wav_in, s_.wav, s_.flux)

    partial_flux = (flux_p - flux_m) / (2*diff_v)
    
    return partial_flux

def cal_partial_f_abun(teff, logg, m_h, wav_start, wav_end, fwhm_broad, vmicro_in, rv_in, wav_in, abun_change, diff_abun, diff_value=0.02, line_list='vald_winered', prefix=''):
    
    if type(diff_abun) != int:
        raise TypeError('Type of diff_abun have to be int.')
    
    if diff_abun not in abun_change.keys():
        raise ValueError('diff_abun have to be in the key of abun_change.')

    # convert rv to wavelength
    del_wav = private.np.abs(rv_in / 3e5 * private.np.mean(wav_in))
        
    abun_change[diff_abun] = abun_change[diff_abun] + diff_value
    s = synth.synth(teff, logg, m_h, wav_start-0.75-del_wav, wav_end+0.75+del_wav, 20000, line_list=line_list, weedout=True, prefix=prefix, vmicro=vmicro_in)
    s.prepare_file(smooth_para=['g', fwhm_broad, 0, 0, 0, 0],
                   abun_change=abun_change)
    s.run_moog()
    s.read_spectra()
    s.wav = s.wav * (1 + rv_in/3e5)
    flux_p = spectres.spectres(wav_in, s.wav, s.flux)

    abun_change[diff_abun] = abun_change[diff_abun] - 2*diff_value
    s_ = synth.synth(teff, logg, m_h, wav_start-0.75-del_wav, wav_end+0.75+del_wav, 20000, line_list=line_list, weedout=True, prefix=prefix, vmicro=vmicro_in)
    s_.prepare_file(smooth_para=['g', fwhm_broad, 0, 0, 0, 0],
                    abun_change=abun_change)
    s_.run_moog()
    s_.read_spectra() 
    s_.wav = s_.wav * (1 + rv_in/3e5)
    flux_m = spectres.spectres(wav_in, s_.wav, s_.flux)
    
    abun_change[diff_abun] = abun_change[diff_abun] + diff_value
    
    partial_flux = (flux_p - flux_m) / (2*diff_value)
    
    return partial_flux

def cal_partial_f_m_h(teff, logg, m_h, wav_start, wav_end, fwhm_broad, vmicro_in, rv_in, wav_in, abun_change=None, diff_m_h=0.02, line_list='vald_winered', prefix=''):

    # convert rv to wavelength
    del_wav = private.np.abs(rv_in / 3e5 * private.np.mean(wav_in))
    
    s = synth.synth(teff, logg, m_h+diff_m_h, wav_start-0.75-del_wav, wav_end+0.75+del_wav, 20000, line_list=line_list, weedout=True, prefix=prefix, vmicro=vmicro_in)
    s.prepare_file(smooth_para=['g', fwhm_broad, 0, 0, 0, 0],
                   abun_change=abun_change)
    s.run_moog()
    s.read_spectra()
    s.wav = s.wav * (1 + rv_in/3e5)
    flux_p = spectres.spectres(wav_in, s.wav, s.flux)

    s_ = synth.synth(teff, logg, m_h-diff_m_h, wav_start-0.75-del_wav, wav_end+0.75+del_wav, 20000, line_list=line_list, weedout=True, prefix=prefix, vmicro=vmicro_in)
    s_.prepare_file(smooth_para=['g', fwhm_broad, 0, 0, 0, 0],
                    abun_change=abun_change)
    s_.run_moog()
    s_.read_spectra() 
    s_.wav = s_.wav * (1 + rv_in/3e5)
    
    flux_m = spectres.spectres(wav_in, s_.wav, s_.flux)
        
    partial_flux = (flux_p - flux_m) / (2*diff_m_h)
    
    return partial_flux

def cal_partial_f_rv(teff, logg, m_h, wav_start, wav_end, fwhm_broad, vmicro_in, rv_in, wav_in, abun_change=None, diff_rv=0.1, line_list='vald_winered', prefix=''):
        
    # convert rv to wavelength
    del_wav = private.np.abs(rv_in / 3e5 * private.np.mean(wav_in))
    
    s = synth.synth(teff, logg, m_h, wav_start-0.75-del_wav, wav_end+0.75+del_wav, 20000, line_list=line_list, weedout=True, prefix=prefix, vmicro=vmicro_in)
    s.prepare_file(smooth_para=['g', fwhm_broad, 0, 0, 0, 0],
                   abun_change=abun_change)
    s.run_moog()
    s.read_spectra()
    s.wav = s.wav * (1 + (rv_in + diff_rv)/3e5)
    flux_p = spectres.spectres(wav_in, s.wav, s.flux)

    s_ = synth.synth(teff, logg, m_h, wav_start-0.75-del_wav, wav_end+0.75+del_wav, 20000, line_list=line_list, weedout=True, prefix=prefix, vmicro=vmicro_in)
    s_.prepare_file(smooth_para=['g', fwhm_broad, 0, 0, 0, 0],
                    abun_change=abun_change)
    s_.run_moog()
    s_.read_spectra() 
    s_.wav = s_.wav * (1 + (rv_in - diff_rv)/3e5)
    
    flux_m = spectres.spectres(wav_in, s_.wav, s_.flux)
         
    partial_flux = (flux_p - flux_m) / (2*diff_rv)
    
    return partial_flux

def mpfit_main(wav_in, flux_in, vmicro_in, fwhm_broad, rv_in, m_h, abun_change_in, fit_paras_list, teff, logg, 
               line_list='vald_3000_24000', niter_max=30, iter_printout=False, fitting_boundary=None, boundary_mode='back', boundary_printout=False, prefix=''):
    '''
    Fit the stellar parameters using the process from MPFIT.
    '''
    
    N = len(wav_in)
    del_x = private.np.ones([len(fit_paras_list),1])

    # Record all the parameters in each run to para_record 
    para_record = {'vmicro':[vmicro_in], 'vbroad':[fwhm_broad], 'm_h':[m_h], 'rv':[rv_in], 'C0':[0]}
    for fit_para in fit_paras_list:
        if 'abun' in fit_para:
             para_record[fit_para] = [abun_change_in[fit_para]]

    # Initilizate the fitting boundary
    fitting_boundary_use = {
        'vbroad':[0.05, 5],
        'm_h':[-5, 1],
        'rv':[-15, 15],
    }
    
    for fit_para in fit_paras_list:
        if 'abun' in fit_para:
             fitting_boundary_use[fit_para] = [-3, 3]
    if fitting_boundary is not None:
        for fit_para in fitting_boundary.keys():
            if fit_para not in fit_paras_list:
                raise ValueError('Key values in fitting_boundary must also present in fit_paras_list.')
            fitting_boundary_use[fit_para] = fitting_boundary[fit_para]

    niter = 0
    while private.np.any(private.np.abs(del_x) > 0.005):

        # Stop if niter >= niter_max
        if niter >= niter_max:
            print('Reached maximum iteration number, iteration stopped and the mean result of the iteration after 20 is set as final result.')
            for k in fit_paras_list:
                if k == 'vbroad':
                    fwhm_broad = private.np.mean(para_record[k][20:])
                elif 'abun' in k: 
                    abun_change_in[k] = private.np.mean(para_record[k][20:])
                elif k == 'm_h':
                    m_h = private.np.mean(para_record[k][20:])
                elif k == 'rv':
                    rv_in = private.np.mean(para_record[k][20:])
            break
        
        # Calculate F_0
        s = synth.synth(teff, logg, m_h, wav_in[0]-0.75, wav_in[-1]+0.75, 20000, line_list=line_list, weedout=True, prefix=prefix, vmicro=vmicro_in)
        s.prepare_file(smooth_para=['g', fwhm_broad, 0, 0, 0, 0], abun_change=abun_change_in)
        s.run_moog()
        s.read_spectra()
        s.wav = s.wav * (1 + rv_in/3e5)

        flux_0 = spectres.spectres(wav_in, s.wav, s.flux)

        # Calculate partial C_0
        C_0 = 1/N * (private.np.sum(flux_in) - private.np.sum(flux_0))
        para_record['C0'].append(C_0)

        # Calculate partial F_0
        # Divide these into different functions.
        #   To do: loggf, C6
        #   Done: Vbroad, Vmicro, m_h, x_fe, rv

        partial_flux_dict = {}
        for k in fit_paras_list:
            if k == 'vbroad':
                partial_flux_dict['vbroad'] = cal_partial_f_v(teff, logg, m_h, wav_in[0], wav_in[-1], fwhm_broad, vmicro_in, rv_in, wav_in, 
                                                              'vbroad', abun_change=abun_change_in, line_list=line_list, prefix=prefix)
            elif k == 'vmicro':
                partial_flux_dict['vmicro'] = cal_partial_f_v(teff, logg, m_h, wav_in[0], wav_in[-1], fwhm_broad, vmicro_in, rv_in, wav_in, 
                                                              'vmicro', abun_change=abun_change_in, line_list=line_list, prefix=prefix)
            elif 'abun' in k: 
                partial_flux_dict[k] = cal_partial_f_abun(teff, logg, m_h, wav_in[0], wav_in[-1], fwhm_broad, vmicro_in, rv_in, wav_in, 
                                                          abun_change_in, int(k.split('_')[1]), line_list=line_list, prefix=prefix)
            elif k == 'm_h':
                partial_flux_dict['m_h'] = cal_partial_f_m_h(teff, logg, m_h, wav_in[0], wav_in[-1], fwhm_broad, vmicro_in, rv_in, wav_in, abun_change=abun_change_in, line_list=line_list, prefix=prefix)
            elif k == 'rv':
                partial_flux_dict['rv'] = cal_partial_f_rv(teff, logg, m_h, wav_in[0], wav_in[-1], fwhm_broad, vmicro_in, rv_in, wav_in, abun_change=abun_change_in, line_list=line_list, prefix=prefix)

            if private.np.max(private.np.abs(partial_flux_dict[k])) <= 0.005 and k != 'rv':
                raise ValueError('The sensitivity of {} on current synthetic spectra is smaller than 0.005, thus the result may not converge. Consider removing this fitting parameter.'.format(k))

        # Set up the size of LHS and RHS
        LHS = private.np.zeros([len(fit_paras_list),len(fit_paras_list)])
        RHS = private.np.zeros([len(fit_paras_list),1])

        i = 0    
        for k_out in fit_paras_list:
            j = 0
            RHS[i,0] = private.np.sum((flux_in - flux_0 - C_0) * partial_flux_dict[k_out])
            for k_in in fit_paras_list:
                LHS[i, j] = private.np.sum(partial_flux_dict[k_in] * partial_flux_dict[k_out]) - 1/N * private.np.sum(partial_flux_dict[k_in]) * private.np.sum(partial_flux_dict[k_out])
                j += 1
            i += 1

        del_x = private.np.linalg.inv(LHS).dot(RHS)
        if niter > 20:
            del_x = 0.7 * del_x

        # Renew the parameters; set to boundary value of the parameter exceed the boundary
        i = 0
        
        # If boundary_mode is 'stop', then stop the whole process if any of the fitting parameter exceed the boundary.
        if boundary_mode == 'stop':
            vbroad_exceed = fwhm_broad < fitting_boundary_use['vbroad'][0] or fwhm_broad > fitting_boundary_use['vbroad'][1]
            abund_exceed = [abun_change_in[int(k.split('_')[1])] < fitting_boundary_use[k][0] or abun_change_in[int(k.split('_')[1])] > fitting_boundary_use[k][1] for k in fit_paras_list if 'abun' in k]
            m_h_exceed = m_h < fitting_boundary_use['m_h'][0] or m_h > fitting_boundary_use['m_h'][1]
            rv_exceed = rv_in < fitting_boundary_use['rv'][0] or rv_in > fitting_boundary_use['rv'][1]
            if vbroad_exceed or abund_exceed or m_h_exceed or rv_exceed:
                print('{} exceeded boundary, MPFIT stopped.'.format(list(private.compress(['vbroad', 'abund', 'm_h', 'rv'], [vbroad_exceed, abund_exceed, m_h_exceed, rv_exceed]))))
                return vmicro_in, fwhm_broad, rv_in, m_h, abun_change_in, C_0, niter, para_record
            
        for k in fit_paras_list:
            if k == 'vbroad':
                fwhm_broad += del_x[i,0]
                if fwhm_broad < fitting_boundary_use['vbroad'][0]:
                    if boundary_printout:
                        print('vbroad < vbroad min boundary. Adjust vbroad to {}'.format(fitting_boundary_use['vbroad'][0]))
                    fwhm_broad = fitting_boundary_use['vbroad'][0]
                elif fwhm_broad > fitting_boundary_use['vbroad'][1]:
                    if boundary_printout:
                        print('vbroad > vbroad min boundary. Adjust vbroad to {}'.format(fitting_boundary_use['vbroad'][1]))
                    fwhm_broad = fitting_boundary_use['vbroad'][1]
                para_record[k].append(fwhm_broad)
            elif k == 'vmicro':
                vmicro_in += del_x[i,0]
                para_record[k].append(vmicro_in)
            elif 'abun' in k: 
                abun_change_in[int(k.split('_')[1])] += del_x[i,0]
                if abun_change_in[int(k.split('_')[1])] < fitting_boundary_use[k][0]:
                    if boundary_printout:
                        print('{0} < {0} baundary. Adjust {0} to {1}'.format(k, fitting_boundary_use[k][0]))
                    abun_change_in[int(k.split('_')[1])] = fitting_boundary_use[k][0]
                elif abun_change_in[int(k.split('_')[1])] > fitting_boundary_use[k][1]:
                    print('{0} < {0} baundary. Adjust {0} to {1}'.format(k, fitting_boundary_use[k][1]))
                    if boundary_printout:
                        abun_change_in[int(k.split('_')[1])] = fitting_boundary_use[k][1]
                para_record[k].append(abun_change_in[int(k.split('_')[1])])
            elif k == 'm_h':
                m_h += del_x[i,0]
                if m_h < fitting_boundary_use['m_h'][0]:
                    if boundary_printout:
                        print('m_h < m_h min boundary. Adjust m_h to {}'.format(fitting_boundary_use['m_h'][0]))
                    m_h = fitting_boundary_use['m_h'][0]
                elif m_h > fitting_boundary_use['m_h'][1]:
                    if boundary_printout:
                        print('m_h > m_h min boundary. Adjust m_h to {}'.format(fitting_boundary_use['m_h'][1]))
                    m_h = fitting_boundary_use['m_h'][1]
                para_record[k].append(m_h)
            elif k == 'rv':
                rv_in += del_x[i, 0]
                if rv_in < fitting_boundary_use['rv'][0]:
                    if boundary_printout:
                        print('rv_in < rv_in min boundary. Adjust rv_in to {}'.format(fitting_boundary_use['rv'][0]))
                    rv_in = fitting_boundary_use['rv'][0]
                elif rv_in > fitting_boundary_use['rv'][1]:
                    if boundary_printout:
                        print('rv_in > rv_in min boundary. Adjust rv_in to {}'.format(fitting_boundary_use['rv'][1]))
                    rv_in = fitting_boundary_use['rv'][1]
                para_record[k].append(rv_in)
            i += 1

        if iter_printout:
            print('Iteration # {}. MAX del_x is: '.format(niter, private.np.max(private.np.abs(del_x))))
            print('vmi:{:.2f}, vb:{:.2f}, rv:{:.2f}, m_h:{:.2f}, abun:{}, , C0:{:.2f}, niter:{}'.format(vmicro_in, fwhm_broad, rv_in, m_h, abun_change_in, C_0, niter))
            
        niter += 1
        
    return vmicro_in, fwhm_broad, rv_in, m_h, abun_change_in, C_0, niter, para_record

def cal_depth(line_lis_in, teff, logg, m_h, vmicro, vbroad, abun_change=None, tqdm_disable=True, prefix=''):
    
    depth_list = []
    for i in private.tqdm(range(len(line_lis_in)), disable=tqdm_disable, leave=False):
#     for i in tqdm(range(10)):
            
        # Extract line wavelength
        line_wav = line_lis_in.iloc[i]['wavelength']
        # Save the linelist
        line_data.save_linelist(line_lis_in.iloc[i:i+1], 'use.list')
        
        s = synth.synth(teff, logg, m_h, line_wav-2, line_wav+2, 20000, line_list='use.list', weedout=False, prefix=prefix, vmicro=vmicro)
        s.prepare_file(smooth_para=['g', vbroad, 0, 0, 0, 0], abun_change=abun_change)
        s.run_moog()
        s.read_spectra()
        
        depth_list.append(1 - private.np.min(s.flux))
        
    line_lis_in['depth'] = depth_list
        
    return line_lis_in

def cal_depth_blending_ratio(line_lis_in, teff, logg, m_h, vmicro, vbroad, abun_change=None, del_wav=2, tqdm_disable=True, prefix=''):
    
    depth_blending_ratio_list = []
    for i in private.tqdm(range(len(line_lis_in)), disable=tqdm_disable, leave=False):
    
        line_wav = line_lis_in.iloc[i]['wavelength']
        s = synth.synth(teff, logg, m_h, line_wav-del_wav, line_wav+del_wav, 50000, 
                               line_list='vald_3000_24000', weedout=False, prefix=prefix, vmicro=vmicro)
        if i == 0:
            s.prepare_file(smooth_para=['g', vbroad, 0, 0, 0, 0], abun_change=abun_change)
            private.copyfile(s.rundir_path + 'model.mod', './model.mod')
        else:
            s.prepare_file(model_file='model.mod', smooth_para=['g', vbroad, 0, 0, 0, 0], abun_change=abun_change)
        s.run_moog()
        s.read_spectra(remove=False)
        wav_all, flux_all = s.wav, s.flux
        
        line_list_syn = line_data.read_linelist(s.rundir_path + 'line.list')
        line_index = line_data.find_lines(line_lis_in.iloc[i:i+1], line_list_syn)
        
        line_list_syn_exclude = line_list_syn.drop(line_index)
        line_data.save_linelist(line_list_syn_exclude, s.rundir_path + 'line.list')
        s.run_moog()
        s.read_spectra()
        wav_exclude, flux_exclude = s.wav, s.flux
        
        
        r_blend_depth = (1-flux_exclude[private.np.argmin(private.np.abs(wav_exclude-line_wav))]) / (1-flux_all[private.np.argmin(private.np.abs(wav_all-line_wav))])

        depth_blending_ratio_list.append(r_blend_depth)

    line_lis_in['f_d_blend'] = depth_blending_ratio_list
    
    return line_lis_in

def cal_X_index(line_list, teff):
    '''
    Calculate the X index defined in Kondo+2019.
    '''
    
    X_index = line_list['loggf'] - line_list['EP'] * (5040 /  (0.86 * teff))
    line_list['X_index'] = X_index
    
    return line_list

def find_isolated_lines(teff, logg, m_h, wav_range_list, ele_index, resolution, obs_spec=None, line_list_in='vald_3000_24000', x_index_thres=-6, depth_obs_thres=0.05):

    '''
    Find the isolated lines of one element within a wavelength range.

    Parameters
    ----------
    teff : float
        The effective temperature of the model
    logg : float
        logg value of the model
    m_h : float
        [M/H] value (overall metallicity) of the model
    wav_range_list : list, [[wav_start_1, wav_end_1], [wav_start_2, wav_end_2], ...]
        List of wavelength range to search
    ele_index : int
        The atomic number of element for the lines to be searched.
    resolution : float
        Resolution of the synthetic spectra; this will passed to MOOG and convolute with initial spectra.
    obs_spec : list of [wav, flux], optional
        The observed spectra.
    line_list_in : str or pd.DataFrame
        The name of the linelist file. If not specified will use built-in VALD linelist (vald_3000_24000).
    x_index_thres : float, default -6
        The threshold of lines to be included (with x-index < x_index_thres)
    depth_obs_thres : float, default 0.05
        The threshold of lines to be included (with obsered depth >= depth_obs_thres)

    Returns
    ----------
    line_list_ele_all : pandas.DataFrame
        The isolated lines list.
    '''

    # Only the Fe lines with depth > depth_lower_limit and blending ratio < blending_limit will be selected.
    depth_lower_limit = 0.05
    blending_limit = 0.2

    j = 0
    for wav_range in wav_range_list:
    # for wav_range in tqdm(wav_range_list[0:1]):

        wav_start, wav_end = wav_range[0], wav_range[1]

        if isinstance(line_list_in, str):
            # Linelist file is not specified, use internal line list;
            line_list = line_data.read_linelist(line_list_in)
        elif isinstance(line_list_in, private.pd.DataFrame):
            line_list = line_list_in
        else:
            raise TypeError('Type of input linelist have to be either str or pandas.DataFrame.')
        
        indices = (line_list['id'] == ele_index) & (line_list['wavelength'] >= wav_start) & (line_list['wavelength'] <= wav_end)
        line_list_ele = line_list[indices].reset_index(drop=True)
    #     pymoog.line_data.save_linelist(line_list_ele, '{}/vald_fe_use.list'.format(folder_name))

        w = weedout.weedout(teff, logg, m_h, 
                                   private.np.min(line_list_ele['wavelength']), np.max(line_list_ele['wavelength']), 
                                   line_list=line_list_ele, kappa_ratio=0.1)
        w.prepare_file()
        w.run_moog()
        w.read_linelist()

        line_list_ele = w.keep_list

        line_list_ele = cal_depth(line_list_ele, teff, logg, m_h, 2, private.np.mean(wav_range)/resolution, tqdm_disable=False)
        line_list_ele = line_list_ele[line_list_ele['depth'] > depth_lower_limit].reset_index(drop=True)

        line_list_ele = cal_depth_blending_ratio(line_list_ele, teff, logg, m_h, 2, private.np.mean(wav_range)/resolution, tqdm_disable=False)
        line_list_ele = line_list_ele[line_list_ele['f_d_blend'] < blending_limit].reset_index(drop=True)

        if j == 0:
            line_list_ele_all = line_list_ele
        else:
            line_list_ele_all = private.pd.concat([line_list_ele_all, line_list_ele])
        j += 1

    line_list_ele_all = line_list_ele_all.reset_index(drop=True)

    # TODO: Calculate the tel_blending ratio
    t = []
    obs_depth = []

    obs_spec = private.pd.read_csv('solar_spec.txt', sep=' ', names=['wav', 'obs_spec'])
    obs_spec['obs_spec'] /= private.np.median(obs_spec['obs_spec'])

    wavs = line_list_ele_all['wavelength'].values
    for wav in wavs:
        obs_indices = private.np.abs(obs_spec['wav']-wav) <= 0.2
        obs_depth.append(1 - private.np.min(obs_spec.loc[obs_indices, 'obs_spec']))
        t.append(wav)

    line_list_ele_all['depth_obs'] = obs_depth

    line_list_ele_all = cal_X_index(line_list_ele_all, teff)
    line_list_ele_all = line_list_ele_all[line_list_ele_all['X_index'] < x_index_thres].reset_index(drop=True)

    line_list_ele_all = line_list_ele_all[line_list_ele_all['depth_obs'] > depth_obs_thres].reset_index(drop=True)
    
    return line_list_ele_all


def determine_vmicro_m_h(spline_res, x_index_dict):
    '''
    Determine the microturbulence and [M/H] from the MPFIT result.
    '''
    slope_list = []
    
    vmicro_list = private.np.linspace(0.5, 2.5, 20)
    for vmicro_index in range(len(vmicro_list)):

        m_h_list = []
        x_index_list = []

        for i in range(len(spline_res)):

            m_h = spline_res[i][1](vmicro_list[vmicro_index])
            x_index = x_index_dict[i][1]

            m_h_list.append(m_h)
            x_index_list.append(x_index)

        x_index_list, m_h_list = private.np.array(x_index_list), private.np.array(m_h_list)    

        fit_res = private.np.polyfit(x_index_list[~private.np.isnan(m_h_list)], m_h_list[~private.np.isnan(m_h_list)], 1)
        slope = fit_res[0]
        slope_list.append(slope)

    slope_list = private.np.array(slope_list)
    
    slope_spline = private.scipy.interpolate.UnivariateSpline(vmicro_list, slope_list, k=3, s=10**-4)
        
    # Calculte vmicro
    vmicro_final = slope_spline.roots()[0]

    # Calculate [M/H] using the determined vmicro.
    m_h_list = []
    for i in range(len(spline_res)):
        m_h = spline_res[i][1](vmicro_final)
        m_h_list.append(m_h)

    m_h_final = private.np.mean(m_h_list)
    
    return vmicro_final, m_h_final, slope_spline