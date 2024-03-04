#!/usr/bin/python
import os
from . import moog_structure
from . import line_data
from . import synth

class weedout(moog_structure.moog_structure):
    def __init__(self, teff, logg, m_h, start_wav, end_wav, vmicro=2, mass=1, kappa_ratio=0.01, line_list='vald_3000_24000', keeplines='keep.list', tosslines='toss.list', prefix='', vmicro_mode='flexible'):
        '''
        Initiate a weedout Instance and read the parameters.
        
        Parameters
        ----------
        teff : int
            The effective temperature of the model
        logg : float
            logg value of the model
        m_h : float
            [M/H] value (overall metallicity) of the model
        start_wav : float
            The start wavelength of line list
        end_wav : float
            The end wavelength of line list
        vmicro : float, default 2
            The microturbulance velocity of the model. 
        mass : float, default 1
            The stellar mass of the input model. Only used when the model type is MARCS spherical.
        kappa_ratio : float, default 0.01
            Minimum line/continuum opacity ratio.
        line_list : str, default 'ges'
            The name of the linelist file. If not specified will use built-in VALD linelist.
        keeplines : str, default 'keep.list'
            The name of the linelist for the lines being kept (in the pymoog working path).
        tosslines : str, default 'toss.list'
            The name of the linelist for the lines tossed (in the pymoog working path).
        prefix : str, default ''.
            The prefix to be added to the name of rundir. Convenient when you want to find a specified rundir if there are many.
        '''
        super(weedout, self).__init__('weedout', prefix=prefix)
        self.teff = teff
        self.logg = logg
        self.m_h = m_h
        self.vmicro = vmicro
        self.mass = mass
        self.start_wav = start_wav
        self.end_wav = end_wav
        self.kappa_ratio = kappa_ratio
        self.line_list_in = line_list
        self.keeplines = keeplines
        self.tosslines = tosslines
        self.vmicro_mode = vmicro_mode

    def read_linelist(self, tosslines=False, remove=True):
        '''
        Read the keep (and toss) linelist of weedout driver.

        Parameters
        ----------
        tosslines : bool, default False
            If True then also output the tossed linelist.
        remove : bool, default True
            Whether remove the working folder after this function.

        Returns
        ---------
        self.keep_list : a panda DataFrame
            A DataFrame for the lines being kept.
        self.toss_list : a panda DataFrame only appear when tosslines is True
            A DataFrame for the lines tossed.
        '''
        keep_list = line_data.read_linelist(self.rundir_path + self.keeplines)
        if tosslines:
            toss_list = line_data.read_linelist(self.rundir_path + self.tosslines)
            self.keep_list, self.toss_list =  keep_list, toss_list
        else:
            self.keep_list = keep_list
            
        if remove:
            self.remove_rundir()
        
    def compare(self, resolution, output=False):
        '''
        Compare the synthetic spectra before and after weedout.
 
        Parameters
        ----------
        resolution : float
            Resolution of the synthetic spectra; this will passed to MOOG and convolute with initial spectra.

        Returns
        ---------
        self.wav_all : a numpy array
            An array of wavelength for the spectra using all the lines input.
        self.flux_all : a numpy array
            An array of flux for the spectra using all the lines input.
        self.wav_keep : a numpy array
            An array of wavelength for the spectra using the lines being kept.
        self.flux_keep : a numpy array
            An array of flux for the spectra using the lines being kept.
        '''
        
        # run synth to get the spectra
        s_all = synth.synth(self.teff, self.logg, self.m_h, self.start_wav, self.end_wav, resolution, line_list=self.rundir_path + 'line.list')
        s_all.prepare_file()
        s_all.run_moog(output=output)
        s_all.read_spectra()
        
        s_keep = synth.synth(self.teff, self.logg, self.m_h, self.start_wav, self.end_wav, resolution, line_list=self.rundir_path + 'keep.list')
        s_keep.prepare_file()
        s_keep.run_moog(output=output)
        s_keep.read_spectra()
        
        self.wav_all, self.flux_all, self.wav_keep, self.flux_keep = s_all.wav, s_all.flux, s_keep.wav, s_keep.flux
        self.remove_rundir()