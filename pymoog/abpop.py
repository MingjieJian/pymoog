import numpy as np
from . import moog_structure
from . import private, model, line_data
import subprocess

'''example batch.par
synpop
modprefix MODEL
synprefix mod
title 47Tuc models
abundances 5
     11 6 7 13 12
isotopes 2
     607.01214 607.01314
models
 1 77.85 3 5.00 6.50 8.50 4.0 6.0 1.0 30.0
 2 70.45 2 5.00 6.50 8.50 4.0 6.0 1.0 30.0
 3 61.74 3 5.00 6.50 8.50 4.0 6.0 1.0 30.0
'''

class synpop(moog_structure.moog_structure):
    def __init__(self, stellar_paras_list, model_RM, start_wav, end_wav, resolution, vmicro=2, mass=1, line_list='vald_3000_24000', weedout=False, prefix='', vmicro_mode='flexible', model_abundances={}, model_isotopes={}):
        '''
        Initiate a synpop instance and read the parameters.
        
        Parameters
        ----------
        stellar_paras_list : list, [[teff1, logg1, m_h1], [teff2, logg2, m_h2], ...]
            List of the stellar parameters of the binaries, max length 99.
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
        super(synpop, self).__init__('synpop', prefix=prefix)
        N_box = len(stellar_paras_list)
        if not hasattr(vmicro, '__len__'):
            vmicro = [vmicro] * N_box
        if not hasattr(mass, '__len__'):
            mass = [mass] * N_box

        self.stellar_paras = stellar_paras_list
        self.vmicro = vmicro
        self.mass = mass
        self.N_box = N_box

        self.start_wav = start_wav
        self.end_wav = end_wav
        self.resolution = resolution
        self.line_list_in = line_list
        self.weedout = weedout
        self.prefix = prefix
        self.vmicro_mode = vmicro_mode

        self.model_RM = model_RM
        self.model_abundances = model_abundances
        self.model_isotopes = model_isotopes

        if start_wav >= end_wav:
            raise ValueError('start_wav has to be smaller than end_wav.')
        if end_wav - start_wav >= 2000:
            raise ValueError('MOOG may provide incorrect spectra when the synthetic length is longer than 2000A. Please split the task into tasks with length <2000 and combine them later on.')        
        
        # Weedout the line list 
        if self.weedout != False:
            self.line_list = []
            for i in range(2):
                if self.weedout == True:
                    w = weedout.weedout(*self.stellar_paras[i], self.start_wav, self.end_wav, line_list=self.rundir_path+self.line_list, prefix=self.prefix, vmicro=self.vmicro[i], mass=self.mass[i])
                else:
                    w = weedout.weedout(*self.stellar_paras[i], self.start_wav, self.end_wav, kappa_ratio=self.weedout, line_list=self.rundir_path+self.line_list, prefix=self.prefix, vmicro=self.vmicro[i], mass=self.mass[i])
                w.prepare_file()
                w.run_moog()
                w.read_linelist()
                self.line_list.append(w.keep_list)

    def prepare_file(self, model_file=None, model_format='moog', loggf_cut=None, abun_change=None, molecules_include=None, model_type='marcs', model_chem='st', model_geo='auto', **args):
        '''
        Prepare the model, linelist and control files for synpop driver.
        Can either provide stellar parameters and wavelengths or provide file names.
        If fine name(s) provided, the files will be copied to working directory for calculation.  
        
        Parameters
        ----------
        model_file : str, optional
            The name of the model file. If not specified, the code will use internal model.
        model_format : str, optional
            The type of the INPUT model file. Default is "moog" (then no conversion of format will be done); can be "moog", "kurucz-atlas9", "kurucz-atlas12" or "marcs". Should left as it is when not providing the input model file. 
        loggf_cut : float, optional
            The cut in loggf; if specified will only include the lines with loggf >= loggf_cut.
        abun_change : dict of pairs {int:float, ...}
            Abundance change, have to be a dict of pairs of atomic number and [X/Fe] values.
        atmosphere : int, default 1
            The atmosphere value described in MOOG documention, section III.
        lines : int, default 1
            The lines value described in MOOG documention, section III.
        molecules : int, default 1
            The molecules value described in MOOG documention, section III.
        molecules_include : list, default None
            Molecules to be included to molecular calculation. Follows the MOOG notation.
        smooth_para : None or list, default None
            The smoothing parameter list of the synthetic spectra.
        model_type : str, default marcs
            The type of internal atmosphere model. Must be kurucz or marcs.
        model_chem : str, default st
            The chemical composition of marcs model. Only valid when model_type is marcs. 
        model_geo : str, default auto
            The geometry of MARCS model, either 's' for spherical, 'p' for plane-parallel or 'auto'.
        '''

        # Use defaule smooth parameter if not specified. 
        if 'smooth_para' not in args.keys():
            args['smooth_para'] = ['g', 0.0, 0.0, 0.0, 0.0, 0.0]

        if type(abun_change) == dict or abun_change == None:
            self.abun_change = [abun_change, abun_change]

        # Create model file.
        self.model_file = []
        self.vmicro_model = []
        for i in range(self.N_box):
            if model_file == None:
                # Model file is not specified, will use builtin model according to stellar parameters.
                self.model1 = model.interpolate_model(*self.stellar_paras[i], vmicro=self.vmicro[i], vmicro_mode=self.vmicro_mode, mass=self.mass[i], abun_change=self.abun_change[i], molecules_include=molecules_include, save_name=self.rundir_path + 'model{}'.format(i+1), model_type=model_type, chem=model_chem, geo=model_geo)
                self.model_file.append('model{}'.format(i+1))
                self.vmicro_model.append(self.model1['vmicro_model'])
            else:
                # Model file is specified; record model file name and copy to working directory.
                if len(model_file) != 2 or type(model_file) != list:
                    raise ValueError('When model file is specified, model_file must be length 2 list.')
                if model_format == 'moog':
                    subprocess.run(['cp', model_file[i], self.rundir_path], encoding='UTF-8', stdout=subprocess.PIPE)
                    self.model_file.append(model_file.split('/')[-1])
                elif model_format[:6] == 'kurucz':
                    model.kurucz2moog(model_path=model_file[i], abun_change=self.abun_change[i], model_format=model_format[7:], molecules_include=molecules_include, converted_model_path=self.rundir_path + 'model{}'.format(i+1))
                    self.model_file.append('model{}'.format(i+1))
                elif model_format == 'marcs':
                    marcs_model = model.read_marcs_model(model_file[i])
                    model.marcs2moog(marcs_model, self.rundir_path + 'model{}'.format(i+1), abun_change=self.abun_change[i], molecules_include=molecules_include)
                    self.model_file.append('model{}'.format(i+1)) 
                else:
                    raise ValueError("The input model_type is not supported. Have to be either 'moog', 'kurucz' or 'marcs.")

        # Create line list.
        if 'del_wav' not in args.keys():
            args['del_wav'] = 0.02
            self.del_wav = args['del_wav']
        if 'del_wav_opac' not in args.keys():
            args['del_wav_opac'] = 1.0
            self.del_wav_opac = args['del_wav_opac']
        if args['del_wav'] < 0.001:
            raise ValueError('del_wav cannot be smaller than 0.001; the calculation and I/O precision is not enough.')
        smooth_width = np.mean([self.start_wav / self.resolution, self.end_wav / self.resolution])
        smooth_width_num = int(np.ceil(smooth_width / args['del_wav']))

        if isinstance(self.line_list_in, str):
            if self.line_list_in[-5:] != '.list':
                # Linelist file is not specified, use internal line list;
                line_list = line_data.read_linelist(self.line_list_in, loggf_cut=loggf_cut)
                line_data.save_linelist(line_list, self.rundir_path + 'line.list', wav_start=self.start_wav-smooth_width_num*2*args['del_wav'], wav_end=self.end_wav+smooth_width_num*2*args['del_wav'])
                self.line_list_name = 'line.list'
                self.line_list = line_list
            elif self.line_list_in[-5:] == '.list':
                # Linelist file is specified; record linelist file name and copy to working directory.
                subprocess.run(['cp', self.line_list_in, self.rundir_path], encoding='UTF-8', stdout=subprocess.PIPE)
                self.line_list_name = self.line_list.split('/')[-1]
                args['lines_in'] = self.line_list_name
                self.line_list = None
        elif isinstance(self.line_list_in, private.pd.DataFrame):
            line_data.save_linelist(self.line_list_in, self.rundir_path + 'line.list', wav_start=self.start_wav-smooth_width_num*2*args['del_wav'], wav_end=self.end_wav+smooth_width_num*2*args['del_wav'])
            self.line_list = 'line.list'
        else:
            raise TypeError('Type of input linelist have to be either str or pandas.DataFrame.')
                
        # Create parameter file.
        self.create_para_file(args=args)

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
            models_file = open(self.rundir_path+'bin_raw.out')
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
            