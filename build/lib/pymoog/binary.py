#!/usr/bin/python
import numpy as np
from . import moog_structure
from . import private

'''example batch.par
binary
bin_raw_out        'bin_raw.out'
bin_smo_out        'bin_smo.out'
atmosphere         1
lines              1
molecules          1
deltaradvel        300
lumratio           1
terminal           'x11'
plot               3
plotpars    1
  0.0  0.0  0.0  0.0 
  0.0  0.0  0.0  0.0 
  g  0.212  0.000  0.000  0.000  0.000
RUN                1
standard_out       'MOOG1.out1'
summary_out        'MOOG1.out2'
smoothed_out       'MOOG1.out3'
model_in           'model1.mod'
lines_in           'line.list'
synlimits
  10539.56  10620.44  0.02  1.00
RUN                 2
standard_out       'MOOG2.out1'
summary_out        'MOOG2.out2'
smoothed_out       'MOOG2.out3'
model_in           'model2.mod'
lines_in           'line.list'
synlimits
  10539.56  10620.44  0.02  1.00
'''

class binary(moog_structure.moog_structure):
    def __init__(self, teff, logg, m_h, start_wav, end_wav, resolution, vmicro=2, mass=1, line_list='vald_3000_24000', weedout=False, prefix='', vmicro_mode='flexible'):
        '''
        Initiate a synth Instance and read the parameters.
        
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
        line_list : str or pd.DataFrame
            The name of the linelist file. If not specified will use built-in VALD linelist (vald_3000_24000).
        weedout : bool or float, default False
            The switch for running weedout driver before synth. If False then weedout is not run; if True the weedout is run with kappa_ratio=0.01, and if a float (> 0 and < 1) is given then weedout is run with the kappa_ratio set as the number.
        prefix : str, default ''.
            The prefix to be added to the name of rundir. Convenient when you want to find a specified rundir if there are many.
        vmicro_mode : str, default 'fix'
            The mode of the vmicro in calculation. If 'fixed', will use the same vmicro in model interpolation and synthesis; if 'flexible', then will use the cloest vmicro in model interpolation if the given vmicro is outside the grid. 
        '''
        super(binary, self).__init__('binary', prefix=prefix)
        self.teff = teff
        self.logg = logg
        self.m_h = m_h
        self.vmicro = vmicro
        self.mass = mass
        self.start_wav = start_wav
        self.end_wav = end_wav
        self.resolution = resolution
        self.line_list = line_list
        self.weedout = weedout
        self.prefix = prefix
        self.vmicro_mode = vmicro_mode

        if start_wav >= end_wav:
            raise ValueError('start_wav has to be smaller than end_wav.')
        if end_wav - start_wav >= 2000:
            raise ValueError('MOOG may provide incorrect spectra when the synthetic length is longer than 2000A. Please split the task into tasks with length <2000 and combine them later on.')        
        
        # Weedout the line list 
        if self.weedout == True:
            if self.weedout == True:
                w = weedout.weedout(self.teff, self.logg, self.m_h, self.start_wav, self.end_wav, line_list=self.rundir_path+self.line_list, prefix=self.prefix)
            else:
                w = weedout.weedout(self.teff, self.logg, self.m_h, self.start_wav, self.end_wav, kappa_ratio=self.weedout, line_list=self.rundir_path+self.line_list, prefix=self.prefix)
            w.prepare_file()
            w.run_moog()
            w.read_linelist()
            self.line_list = w.keep_list

    def read_spectra(self, spec_type='smooth', remove=True):
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
        if spec_type == 'standard':
            models_file = open(self.rundir_path+'MOOG.out2')
            models = models_file.readline()
            models = models_file.readline()
            models = models_file.read().split()
            models = [float(i) for i in models]
            x = range(round((models[1]-models[0])/models[2])+1)
            model_wav = []
            for i in x:
                model_wav.append(models[0] + models[2]*i)
            model_flux = np.array(models[4:])
            self.wav, self.flux = np.array(model_wav), np.array(model_flux)
        elif spec_type == 'smooth':
            models_file = open(self.rundir_path+'MOOG.out3')
            models = models_file.readline()
            models = models_file.readline()
            models = models_file.readlines()
            wavelength = []
            depth = []
            for i in models:
                temp = i.split()
                wavelength.append(float(temp[0]))
                depth.append(float(temp[1]))
            self.wav, self.flux =  np.array(wavelength), np.array(depth)
            
        # Crop the spectra to fit the synthetic wavelength.
        indices = (self.wav >= self.start_wav) & (self.wav <= self.end_wav) 
        self.wav = self.wav[indices]
        self.flux = self.flux[indices]    
            
        if remove:
            self.remove_rundir()
            