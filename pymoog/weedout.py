#!/usr/bin/python
import os
import subprocess
import numpy as np
import re
from . import line_data
from . import model
from . import synth
from . import rundir_num
from . import private

MOOG_path = '{}/.pymoog/moog_nosm/moog_nosm_NOV2019/'.format(os.environ['HOME'])
# self.rundir_path = '{}/.pymoog/rundir/'.format(os.environ['HOME'])
MOOG_file_path = '{}/.pymoog/files/'.format(os.environ['HOME'])

class weedout(rundir_num.rundir_num):
    def __init__(self, teff, logg, m_h, start_wav, end_wav, kappa_ratio=0.01, line_list='ges', keeplines='keep.list', tosslines='toss.list', prefix=''):
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
        kappa_ratio : float, default 0.01
            Minimum line/continuum opacity ratio.
        line_list : str, default 'ges'
            The name of the linelist file. If not specified will use built-in VALD linelist.
        keeplines : str, default 'keep.list'
            The name of the linelist for the lines being kept (in the pymoog working path).
        tosslines : str, default 'toss.list'
            The name of the linelist for the lines tossed (in the pymoog working path).
        '''
        super(weedout, self).__init__('{}/.pymoog/'.format(private.os.environ['HOME']), 'weedout', prefix=prefix)
        self.teff = teff
        self.logg = logg
        self.m_h = m_h
        self.start_wav = start_wav
        self.end_wav = end_wav
        self.kappa_ratio = kappa_ratio
        self.line_list = line_list
        self.keeplines = keeplines
        self.tosslines = tosslines
        
    def prepare_file(self, model_file=None, model_type='moog', loggf_cut=None, abun_change=None, molecules=None, atmosphere=1, lines=1):
        '''
        Prepare the model, linelist and control files for MOOG.
        Can either provide stellar parameters and wavelengths or provide file names.
        If fine name(s) provided, the files will be copied to working directory for calculation.  
        
        Parameters
        ----------
        model_file : str, optional
            The name of the model file. If not specified will use internal Kurucz model.
             
        model_type : str, optional
            The type of the model file. Default is "moog" (then no conversion of format will be done); can be "moog", "kurucz-atlas9" and "kurucz-atlas12". 
        
        logf_cut : float, optional
            The cut in loggf; if specified will only include the lines with loggf >= loggf_cut.
            
        abun_change : dict of pairs {int:float, ...}
            Abundance change, have to be a dict of pairs of atomic number and [X/Fe] values.
        '''
        
        if model_file == None:
            # Model file is not specified, will download Kurucz model according to stellar parameters.
            model.interpolate_model(self.teff, self.logg, self.m_h, abun_change=abun_change, molecules=molecules, to_path=self.rundir_path + 'model.mod')
            self.model_file = 'model.mod'
        else:
            # Model file is specified; record model file name and copy to working directory.
            if model_type == 'moog':
                subprocess.run(['cp', model_file, self.rundir_path], encoding='UTF-8', stdout=subprocess.PIPE)
                self.model_file = model_file.split('/')[-1]
            elif model_type[:6] == 'kurucz':
                model.KURUCZ_convert(model_path=model_file, abun_change=abun_change, model_type=model_type[7:], molecules=molecules, converted_model_path=self.rundir_path + 'model.mod')
                self.model_file = 'model.mod'

        if self.line_list[-5:] != '.list':
            line_list = line_data.read_linelist(self.line_list, loggf_cut=loggf_cut, mode='npy')
            line_data.save_linelist(line_list, self.rundir_path + 'line.list', wav_start=self.start_wav, wav_end=self.end_wav)
            self.line_list = 'line.list'
        elif self.line_list[-5:] == '.list':
            # Linelist file is specified; record linelist file name and copy to working directory.
            line_list = line_data.read_linelist(self.line_list, loggf_cut=loggf_cut)
            line_data.save_linelist(line_list, self.rundir_path + 'line.list', wav_start=self.start_wav, wav_end=self.end_wav)
            self.line_list = 'line.list'
            # subprocess.run(['cp', self.line_list, self.rundir_path], encoding='UTF-8', stdout=subprocess.PIPE)
            self.line_list = self.line_list.split('/')[-1]
            
        # Create parameter file.
        self.create_para_file(atmosphere=atmosphere, lines=lines)    
        
    def create_para_file(self, atmosphere=1, lines=1, molecules=2):
        '''
        Function for creating the parameter file of batch.par
        
        Parameters
        ----------
        del_wav : float, optional
            The sampling distance of weedout spectra. Default 0.02.
        smooth : str, optional
            Line profile to be used to smooth the weedout spectra, same as decribed in MOOG documention. Default Gaussian. 
        '''
        MOOG_para_file = open(self.rundir_path + '/batch.par', 'w')
        # Parameter list of MOOG: standard output file (1), summary output file (2), smoothed output file (3),
        #                         begin wavelength, end wavelength, wavelength step;
        #                         smoothing function, Gaussian FWHM, vsini, limb darkening coefficient,
        #                         Macrotrubulent FWHM, Lorentzian FWHM

        MOOG_contant = ["weedout\n",
                        "standard_out       '{}'\n".format('MOOG.out1'),
                        "model_in           '{}'\n".format(self.model_file),
                        "lines_in           '{}'\n".format(self.line_list),
                        "keeplines_out      '{}'\n".format(self.keeplines),
                        "tosslines_out      '{}'\n".format(self.tosslines),
                        "atmosphere         {}\n".format(atmosphere),
                        "lines              {}\n".format(lines),
                        "molecules          {}\n".format(molecules),
                        "terminal           'x11'\n",
                    ]
        MOOG_para_file.writelines(MOOG_contant)
        MOOG_para_file.close()
    
    def run_moog(self, output=False, unlock=False):
        '''
        Run MOOG and print the result if required.

        Parameters
        ----------
        output : boolen, default False
            If set to True, then print the out put of MOOG.

        Returns
        ----------
        None. Three files MOOG.out1, MOOG.out2 and MOOG.out3 will be save in the pymoog working path.
        '''
        
        MOOG_run = subprocess.run([MOOG_path + '/MOOGSILENT'], stdout=subprocess.PIPE, input=bytes('{}'.format(self.kappa_ratio), 'utf-8'), cwd=self.rundir_path)

        
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
        
        # Move line.list as all.list
        subprocess.run(['mv', self.rundir_path + 'line.list', self.rundir_path + 'all.list'])
            
        if output:    
            for i in MOOG_output:
                print(i)
        
        if 'ERROR' in ''.join(MOOG_run):
            raise ValueError('There is error during the running of MOOG.')

    def read_linelist(self, tosslines=False, remove=True):
        '''
        Read the keep (and toss) linelist of weedout driver.

        Parameters
        ----------
        tosslines : bool, default False
            If True then also output the tossed linelist.

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
            self.remove()
        
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
        s_all = synth.synth(self.teff, self.logg, self.m_h, self.start_wav, self.end_wav, resolution, line_list=self.rundir_path + 'all.list')
        s_all.prepare_file()
        s_all.run_moog(output=output)
        s_all.read_spectra()
        
        s_keep = synth.synth(self.teff, self.logg, self.m_h, self.start_wav, self.end_wav, resolution, line_list=self.rundir_path + 'keep.list')
        s_keep.prepare_file()
        s_keep.run_moog(output=output)
        s_keep.read_spectra()
        
        self.wav_all, self.flux_all, self.wav_keep, self.flux_keep =  s_all.wav, s_all.flux, s_keep.wav, s_keep.flux
        