#!/usr/bin/python
from multiprocessing.sharedctypes import Value
import subprocess
import numpy as np
import re
from . import line_data
from . import model
from . import rundir_num
from . import weedout
from . import private, internal
import time

MOOG_path = '{}/.pymoog/moog_nosm/moog_nosm_NOV2019/'.format(private.os.environ['HOME'])
# self.rundir_path = r.rundir_path
MOOG_file_path = '{}/.pymoog/files/'.format(private.os.environ['HOME'])

class synth(rundir_num.rundir_num):
    def __init__(self, teff, logg, m_h, start_wav, end_wav, resolution, vmicro=2, mass=1, del_wav=0.02, line_list='vald_3000_24000', weedout=False, prefix=''):
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
        '''
        super(synth, self).__init__('{}/.pymoog/'.format(private.os.environ['HOME']), 'synth', prefix=prefix)
        self.teff = teff
        self.logg = logg
        self.m_h = m_h
        self.vmicro = vmicro
        self.mass = mass
        self.start_wav = start_wav
        self.end_wav = end_wav
        self.resolution = resolution
        self.del_wav = del_wav
        self.line_list = line_list
        self.weedout = weedout
        self.prefix = prefix

        # Perform some sanity check
        if del_wav < 0.001:
            raise ValueError('del_wav cannot be smaller than 0.001; the calculation and I/O precision is not enough.')

        if start_wav >= end_wav:
            raise ValueError('start_wav has to be smaller than end_wav.')
        if end_wav - start_wav >= 2000:
            raise ValueError('MOOG may provide incorrect spectra when the synthetic length is longer than 2000A. Please split the task into tasks with length <2000 and combine them later on.')

    def prepare_file(self, model_file=None, model_format='moog', loggf_cut=None, abun_change=None, atmosphere=1, lines=1, molecules=1, molecules_include=None, smooth_para=None, model_type='marcs', model_chem='st', model_geo='auto'):
        '''
        Prepare the model, linelist and control files for MOOG.
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
        if smooth_para is None:
            smooth_para = ['g', 0.0, 0.0, 0.0, 0.0, 0.0]
        
        # Create model file.
        if model_file == None:
            # Model file is not specified, will use Kurucz model according to stellar parameters.
            model.interpolate_model(self.teff, self.logg, self.m_h, vmicro=self.vmicro, mass=self.mass, abun_change=abun_change, molecules_include=molecules_include, save_name=self.rundir_path + 'model.mod', model_type=model_type, chem=model_chem, geo=model_geo)
            self.model_file = 'model.mod'
        else:
            # Model file is specified; record model file name and copy to working directory.
            if model_format == 'moog':
                subprocess.run(['cp', model_file, self.rundir_path], encoding='UTF-8', stdout=subprocess.PIPE)
                self.model_file = model_file.split('/')[-1]
            elif model_format[:6] == 'kurucz':
                model.kurucz2moog(model_path=model_file, abun_change=abun_change, model_format=model_format[7:], molecules_include=molecules_include, converted_model_path=self.rundir_path + 'model.mod')
                self.model_file = 'model.mod'
            elif model_format == 'marcs':
                marcs_model = model.read_marcs_model(model_file)
                model.marcs2moog(marcs_model, self.rundir_path + 'model.mod', abun_change=abun_change, molecules_include=molecules_include)
                self.model_file = 'model.mod'
            else:
                raise ValueError("The input model_type is not supported. Have to be either 'moog', 'kurucz' or 'marcs.")

        # Create line list.
        if isinstance(self.line_list, str):
            if self.line_list[-5:] != '.list':
                # Linelist file is not specified, use internal line list;
                line_list = line_data.read_linelist(self.line_list, loggf_cut=loggf_cut)
                line_data.save_linelist(line_list, self.rundir_path + 'line.list', wav_start=self.start_wav, wav_end=self.end_wav)
                self.line_list = 'line.list'
            elif self.line_list[-5:] == '.list':
                # Linelist file is specified; record linelist file name and copy to working directory.
                subprocess.run(['cp', self.line_list, self.rundir_path], encoding='UTF-8', stdout=subprocess.PIPE)
                self.line_list = self.line_list.split('/')[-1]
        elif isinstance(self.line_list, private.pd.DataFrame):
            line_data.save_linelist(self.line_list, self.rundir_path + 'line.list', wav_start=self.start_wav, wav_end=self.end_wav)
            self.line_list = 'line.list'
        else:
            raise TypeError('Type of input linelist have to be either str or pandas.DataFrame.')
            
        # Weedout the line list 
        if self.weedout == True:
            if self.weedout == True:
                w = weedout.weedout(self.teff, self.logg, self.m_h, self.start_wav, self.end_wav, line_list=self.rundir_path+self.line_list, prefix=self.prefix)
            else:
                w = weedout.weedout(self.teff, self.logg, self.m_h, self.start_wav, self.end_wav, kappa_ratio=self.weedout, line_list=self.rundir_path+self.line_list, prefix=self.prefix)
            w.prepare_file()
            w.run_moog()
            w.read_linelist()
            line_data.save_linelist(w.keep_list, self.rundir_path + self.line_list)
            self.keep_list = w.keep_list
                
        # Create parameter file.
        self.create_para_file(atmosphere=atmosphere, lines=lines, molecules=molecules, del_wav=self.del_wav, smooth_para=smooth_para)
        
    def create_para_file(self, del_wav=0.02, del_wav_opac=1.0, smooth_para=['g', 0.0, 0.0, 0.0, 0.0, 0.0], atmosphere=1, lines=1, molecules=2):
        '''
        Function for creating the parameter file of batch.par
        
        Parameters
        ----------
        del_wav : float, default 0.02
            The sampling distance of synthetic spectra. 
        del_wav_opac : float, default 1.0
            The delta wavelength from a spectrum point to consider opacity contributions from neighboring transitions.
        smooth_para : list, default ['g', 0.0, 0.0, 0.0, 0.0, 0.0]
            Smoothing parameters in the third line of plotpars, as decribed in MOOG documention. Note that if the second value j=0, then it will be changed to the width corresponding to the resolution. 
        atmosphere : int, default 1
            The atmosphere value described in MOOG documention, section III.
        lines : int, default 1
            The lines value described in MOOG documention, section III.
        molecules : int, default 1
            The molecules value described in MOOG documention, section III.
        '''
        MOOG_para_file = open(self.rundir_path + '/batch.par', 'w')
        # Parameter list of MOOG: standard output file (1), summary output file (2), smoothed output file (3),
        #                         begin wavelength, end wavelength, wavelength step;
        #                         smoothing function, Gaussian FWHM, vsini, limb darkening coefficient,
        #                         Macrotrubulent FWHM, Lorentzian FWHM
        smooth_width = np.mean([self.start_wav / self.resolution, self.end_wav / self.resolution])

        if smooth_para[1] == 0:
            smooth_para[1] = smooth_width

        # The fitting range is enlarge by smooth_para[1]*2 A, to cover a full range of wavelength without being cut by the smoothing function.
        MOOG_contant = ["synth\n",
                        "standard_out       '{}'\n".format('MOOG.out1'),
                        "summary_out        '{}'\n".format('MOOG.out2'),
                        "smoothed_out       '{}'\n".format('MOOG.out3'),
                        "model_in           '{}'\n".format(self.model_file),
                        "lines_in           '{}'\n".format(self.line_list),
                        "atmosphere         {}\n".format(atmosphere),
                        "lines              {}\n".format(lines),
                        "molecules          {}\n".format(molecules),
                        "terminal           'x11'\n",
                        "synlimits\n",
                        "  {:.2f}  {:.2f}  {}  {:.2f}\n".format(self.start_wav - smooth_para[1]*2, self.end_wav + smooth_para[1]*2, del_wav, del_wav_opac),
                        "plot        3\n",
                        "plotpars    1\n",
                        "  0.0  0.0  0.0  0.0 \n",
                        "  0.0  0.0  0.0  0.0 \n",
                        "  {}  {:.3f}  {:.3f}  {:.3f}  {:.3f}  {:.3f}\n".format(*smooth_para)
                    ]
        MOOG_para_file.writelines(MOOG_contant)
        MOOG_para_file.close()
    
    def run_moog(self, output=False):
        '''
        Run MOOG and print the reuslt if required.

        Parameters
        ----------
        output : boolen, default False
            If set to True, then print the out put of MOOG.

        Returns
        ----------
        None. Three files MOOG.out1, MOOG.out2 and MOOG.out3 will be save in the pymoog working path.
        '''

        MOOG_run = subprocess.run([MOOG_path + '/MOOGSILENT'], stdout=subprocess.PIPE,
                                  cwd=self.rundir_path)
        
        MOOG_run = str(MOOG_run.stdout, encoding = "utf-8").split('\n')
        MOOG_output = []
        for i in MOOG_run:
            if len(i) > 12:
                ansi_escape = re.compile(r'\x1b\[...H')
                temp = ansi_escape.sub('', i)
                ansi_escape = re.compile(r'\x1b\[....H')
                temp = ansi_escape.sub('', temp)
                ansi_escape = re.compile(r'\x1b\[H')
                temp = ansi_escape.sub('', temp)
                ansi_escape = re.compile(r'\x1b\[2J')
                MOOG_output.append(ansi_escape.sub('', temp))
                
        if output:
            for i in MOOG_output:
                print(i)
        
        if 'ERROR' in ''.join(MOOG_run):
            raise ValueError('There is error during the running of MOOG.')
        
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
            self.remove()
            
    def read_model(self, remove=True):
        '''
        Read the output model of MOOG. This model have tauref calculated from MOOG.

        Parameters
        ----------
        remove : bool, default True
            Whether remove the working folder after this function.

        Returns
        ---------
        model_df : pandas DataFrame
            An DataFrame containing the model
        '''
        
        with open(self.rundir_path+'MOOG.out1') as file:
            content = file.readlines()
        i_list = []
        for i in range(len(content)):
            if 'INPUT ATMOSPHERE QUANTITIES' in content[i] or 'INPUT ABUNDANCES:' in content[i]:
                i_list.append(i)

        i_list[0] += 1
        i_list[1] = len(content) - i_list[1]

        self.model = private.pd.read_csv(self.rundir_path+'MOOG.out1', skiprows=i_list[0], skipfooter=i_list[1], sep=' +', engine='python')
        
        for column in ['tauref', 'Pgas', 'Ne', 'Vturb']:
            self.model[column] = self.model[column].map(private.D2E)
            
        self.model = self.model.astype(np.float64)
        
        if remove:
            self.remove()