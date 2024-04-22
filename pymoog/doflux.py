#!/usr/bin/python
from . import moog_structure
from . import private

class doflux(moog_structure.moog_structure):
    def __init__(self, teff, logg, m_h, start_wav, end_wav, resolution, vmicro=2, mass=1, prefix='', vmicro_mode='flexible'):
        '''
        Initiate a doflux Instance and read the parameters.
        
        Parameters
        ----------
        teff : float
            The effective temperature of the model
        logg : float
            logg value of the model
        m_h : float
            [M/H] value (overall metallicity) of the model
        start_wav : float
            The start wavelength of synthetic spectra
        end_wav : float
            The end wavelength of synthetic spectra
        resolution : float
            Resolution of the synthetic spectra; this will passed to MOOG and convolute with initial spectra.
        vmicro : float, default 2
            The microturbulance velocity of the model. 
        mass : float, default 1
            The stellar mass of the input model. Only used when the model type is MARCS spherical.
        del_wav : float, default 0.02
            The wavelength step of the synthetic spectra. 
        weedout : bool or float, default False
            The switch for running weedout driver before synth. If False then weedout is not run; if True the weedout is run with kappa_ratio=0.01, and if a float (> 0 and < 1) is given then weedout is run with the kappa_ratio set as the number.
        prefix : str, default ''.
            The prefix to be added to the name of rundir. Convenient when you want to find a specified rundir if there are many.
        vmicro_mode : str, default 'fix'
            The mode of the vmicro in calculation. If 'fixed', will use the same vmicro in model interpolation and synthesis; if 'flexible', then will use the cloest vmicro in model interpolation if the given vmicro is outside the grid. 
        '''
        super(doflux, self).__init__('doflux', prefix=prefix)
        self.teff = teff
        self.logg = logg
        self.m_h = m_h
        self.vmicro = vmicro
        self.mass = mass
        self.start_wav = start_wav
        self.end_wav = end_wav
        self.resolution = resolution
        self.prefix = prefix
        self.vmicro_mode = vmicro_mode        

        if start_wav >= end_wav:
            raise ValueError('start_wav has to be smaller than end_wav.')
        if end_wav - start_wav >= 2000:
            raise ValueError('MOOG may provide incorrect spectra when the synthetic length is longer than 2000A. Please split the task into tasks with length <2000 and combine them later on.')        

    def read_spectra(self, remove=True):
        '''
        Read the output spectra of MOOG.

        Parameters
        ----------
        spec_type : str, default 'smooth'
            Decide the type of spectra to be read. 'smooth' will read the smoothed spectra (MOOG.out3), and 'standard' will read the un-smoothed one (MOOG.out2).
        remove : bool, default True
            Whether remove the working folder after this function.

        Returns
        ---------
        wav : a numpy array
            An array of wavelength
        flux : a numpy array
            An array of flux
        '''
        flux_df = private.pd.read_csv(self.rundir_path+'MOOG.out2', sep=' +', names=['wave', 'flux', 'waveinv', 'fluxlog'], engine='python')
        indices = (flux_df['wave'] >= self.start_wav) & (flux_df['wave'] <= self.end_wav) 
        flux_df = flux_df[indices]
        self.wav, self.flux_cont, self.waveinv, self.logflux_cont = flux_df['wave'].values, flux_df['flux'].values, flux_df['waveinv'].values, flux_df['fluxlog'].values
            
        # Crop the spectra to fit the synthetic wavelength.
            
        if remove:
            self.remove_rundir()