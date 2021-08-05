from . import private
from . import synth
import spectres

def cal_partial_f_v(teff, logg, m_h, wav_start, wav_end, vbroad_in, vmicro_in, rv_in, wav_in, mode, abun_change=None, diff_v=0.01, line_list='vald_winered'):
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

    s = synth.synth(teff, logg, m_h, wav_start-0.5, wav_end+0.5, 20000, line_list=line_list, weedout=True)
    s.prepare_file(vmicro=vmicro_in+diff_v_dict['vmicro'], smooth_para=['g', vbroad_in+diff_v_dict['vbroad'], 0, 0, 0, 0],
                   abun_change=abun_change)
    s.run_moog()
    s.read_spectra()
    s.wav = s.wav * (1 + rv_in/3e5)
    flux_p = spectres.spectres(wav_in, s.wav, s.flux)

    s_ = synth.synth(teff, logg, m_h, wav_start-0.5, wav_end+0.5, 20000, line_list=line_list, weedout=True)
    s_.prepare_file(vmicro=vmicro_in-diff_v_dict['vmicro'], smooth_para=['g', vbroad_in-diff_v_dict['vbroad'], 0, 0, 0, 0],
                    abun_change=abun_change)
    s_.run_moog()
    s_.read_spectra()
    s_.wav = s_.wav * (1 + rv_in/3e5)
    flux_m = spectres.spectres(wav_in, s_.wav, s_.flux)

    partial_flux = (flux_p - flux_m) / (2*diff_v)
    
    return partial_flux

def cal_partial_f_abun(teff, logg, m_h, wav_start, wav_end, vbroad_in, vmicro_in, rv_in, wav_in, abun_change, diff_abun, diff_value=0.02, line_list='vald_winered'):
    
    if type(diff_abun) != int:
        raise TypeError('Type of diff_abun have to be int.')
    
    if diff_abun not in abun_change.keys():
        raise ValueError('diff_abun have to be in the key of abun_change.')
    
    abun_change[diff_abun] = abun_change[diff_abun] + diff_value
    s = synth.synth(teff, logg, m_h, wav_start-0.5, wav_end+0.5, 20000, line_list=line_list, weedout=True)
    s.prepare_file(vmicro=vmicro_in, smooth_para=['g', vbroad_in, 0, 0, 0, 0],
                   abun_change=abun_change)
    s.run_moog()
    s.read_spectra()
    s.wav = s.wav * (1 + rv_in/3e5)
    flux_p = spectres.spectres(wav_in, s.wav, s.flux)

    abun_change[diff_abun] = abun_change[diff_abun] - 2*diff_value
    s_ = synth.synth(teff, logg, m_h, wav_start-0.5, wav_end+0.5, 20000, line_list=line_list, weedout=True)
    s_.prepare_file(vmicro=vmicro_in, smooth_para=['g', vbroad_in, 0, 0, 0, 0],
                    abun_change=abun_change)
    s_.run_moog()
    s_.read_spectra() 
    s_.wav = s_.wav * (1 + rv_in/3e5)
    flux_m = spectres.spectres(wav_in, s_.wav, s_.flux)
    
    abun_change[diff_abun] = abun_change[diff_abun] + diff_value
    
    partial_flux = (flux_p - flux_m) / (2*diff_value)
    
    return partial_flux

def cal_partial_f_m_h(teff, logg, m_h, wav_start, wav_end, vbroad_in, vmicro_in, rv_in, wav_in, abun_change=None, diff_m_h=0.02, line_list='vald_winered'):
    
    s = synth.synth(teff, logg, m_h+diff_m_h, wav_start-0.5, wav_end+0.5, 20000, line_list=line_list, weedout=True)
    s.prepare_file(vmicro=vmicro_in, smooth_para=['g', vbroad_in, 0, 0, 0, 0],
                   abun_change=abun_change)
    s.run_moog()
    s.read_spectra()
    s.wav = s.wav * (1 + rv_in/3e5)
    flux_p = spectres.spectres(wav_in, s.wav, s.flux)

    s_ = synth.synth(teff, logg, m_h-diff_m_h, wav_start-0.5, wav_end+0.5, 20000, line_list=line_list, weedout=True)
    s_.prepare_file(vmicro=vmicro_in, smooth_para=['g', vbroad_in, 0, 0, 0, 0],
                    abun_change=abun_change)
    s_.run_moog()
    s_.read_spectra() 
    s_.wav = s_.wav * (1 + rv_in/3e5)
    flux_m = spectres.spectres(wav_in, s_.wav, s_.flux)
        
    partial_flux = (flux_p - flux_m) / (2*diff_m_h)
    
    return partial_flux

def cal_partial_f_rv(teff, logg, m_h, wav_start, wav_end, vbroad_in, vmicro_in, rv_in, wav_in, abun_change=None, diff_rv=0.1, line_list='vald_winered'):
        
    # convert rv to wavelength
    del_wav = rv_in / 3e5 * private.np.mean(wav_in)
    
    s = synth.synth(teff, logg, m_h, wav_start-0.5-del_wav, wav_end+0.5+del_wav, 20000, line_list=line_list, weedout=True)
    s.prepare_file(vmicro=vmicro_in, smooth_para=['g', vbroad_in, 0, 0, 0, 0],
                   abun_change=abun_change)
    s.run_moog()
    s.read_spectra()
    s.wav = s.wav * (1 + (rv_in + diff_rv)/3e5)
    flux_p = spectres.spectres(wav_in, s.wav, s.flux)

    s_ = synth.synth(teff, logg, m_h, wav_start-0.5-del_wav, wav_end+0.5+del_wav, 20000, line_list=line_list, weedout=True)
    s_.prepare_file(vmicro=vmicro_in, smooth_para=['g', vbroad_in, 0, 0, 0, 0],
                    abun_change=abun_change)
    s_.run_moog()
    s_.read_spectra() 
    s_.wav = s_.wav * (1 + (rv_in - diff_rv)/3e5)
    flux_m = spectres.spectres(wav_in, s_.wav, s_.flux)
        
    partial_flux = (flux_p - flux_m) / (2*diff_rv)
    
    return partial_flux

def mpfit_main(wav_in, flux_in, vmicro_in, vbroad_in, rv_in, m_h, abun_change_in, fit_paras_list, teff, logg, 
               line_list='vald_winered', niter_max=15, printout=False):
    '''
    Fit the stellar parameters using the process from MPFIT.
    '''
    
    N = len(wav_in)
    del_x = private.np.ones([len(fit_paras_list),1])

    niter = 0
    while private.np.any(private.np.abs(del_x) > 0.005):

        # Stop if niter >= niter_max
        if niter >= niter_max:
            raise ValueError('Reached maximum iteration number, stopped.')
        
        # Calculate F_0
        s = synth.synth(teff, logg, m_h, wav_in[0]-0.5, wav_in[-1]+0.5, 20000, line_list=line_list, weedout=True)
        s.prepare_file(vmicro=vmicro_in, smooth_para=['g', vbroad_in, 0, 0, 0, 0], abun_change=abun_change_in)
        s.run_moog()
        s.read_spectra()
        s.wav = s.wav * (1 + rv_in/3e5)

        flux_0 = spectres.spectres(wav_in, s.wav, s.flux)

        # Calculate partial C_0
        C_0 = 1/N * (private.np.sum(flux_in) - private.np.sum(flux_0))

        # Calculate partial F_0
        # Divide these into different functions.
        #   To do: loggf, C6
        #   Done: Vbroad, Vmicro, m_h, x_fe, rv

        partial_flux_dict = {}
        for k in fit_paras_list:
            if k == 'vbroad':
                partial_flux_dict['vbroad'] = cal_partial_f_v(teff, logg, m_h, wav_in[0], wav_in[-1], vbroad_in, vmicro_in, rv_in, wav_in, 
                                                              'vbroad', abun_change=abun_change_in, line_list=line_list)
            elif k == 'vmicro':
                partial_flux_dict['vmicro'] = cal_partial_f_v(teff, logg, m_h, wav_in[0], wav_in[-1], vbroad_in, vmicro_in, rv_in, wav_in, 
                                                              'vmicro', abun_change=abun_change_in, line_list=line_list)
            elif 'abun' in k: 
                partial_flux_dict[k] = cal_partial_f_abun(teff, logg, m_h, wav_in[0], wav_in[-1], vbroad_in, vmicro_in, rv_in, wav_in, 
                                                          abun_change_in, int(k.split('_')[1]), line_list=line_list)
            elif k == 'm_h':
                partial_flux_dict['m_h'] = cal_partial_f_m_h(teff, logg, m_h, wav_in[0], wav_in[-1], vbroad_in, vmicro_in, rv_in, wav_in, abun_change=abun_change_in, line_list=line_list)
            elif k == 'rv':
                partial_flux_dict['rv'] = cal_partial_f_rv(teff, logg, m_h, wav_in[0], wav_in[-1], vbroad_in, vmicro_in, rv_in, wav_in, abun_change=abun_change_in, line_list=line_list)

            if private.np.max(private.np.abs(partial_flux_dict[k])) <= 0.01 and k != 'rv':
                raise ValueError('The sensitivity of {} on current synthetic spectra is smaller than 0.01, thus the result may not converge. Consider removing this fitting parameter.'.format(k))

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

        # Renew the parameters
        i = 0
        for k in fit_paras_list:
            if k == 'vbroad':
                vbroad_in += del_x[i,0]
                if vbroad_in < 0:
                    print('adjust vbroad to 0.05')
                    vbroad_in = 0.05
            elif k == 'vmicro':
                vmicro_in += del_x[i,0]
            elif 'abun' in k: 
                abun_change_in[int(k.split('_')[1])] += del_x[i,0]
            elif k == 'm_h':
                m_h += del_x[i,0]
            elif k == 'rv':
                rv_in += del_x[i, 0]
            i += 1

        niter += 1
        if printout:
            print('vmi:{:.2f}, vb:{:.2f}, rv:{:.2f}, m_h:{:.2f}, abun:{}, , C0:{:.2f}, niter:{}'.format(vmicro_in, vbroad_in, rv_in, m_h, abun_change_in, C_0, niter))
            print('MAX del_x is: ', private.np.max(private.np.abs(del_x)))
    return vmicro_in, vbroad_in, rv_in, m_h, abun_change_in, C_0, niter