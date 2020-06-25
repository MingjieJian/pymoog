import os
import subprocess
import numpy as np
import re

MOOG_path = '{}/.pymoog/moog_nosm/moog_nosm_FEB2017/'.format(os.environ['HOME'])
MOOGrun_path = '{}/.pymoog/rundir/'.format(os.environ['HOME'])

class synth:
    def __init__(self, teff, logg, m_h, start_wav, end_wav, resolution, del_wav=0.02, smooth='g'):
        '''
        Initiate a SpecSyntheiszer and read the parameters.
        '''
        self.teff = teff
        self.logg = logg
        self.m_h = m_h
        self.start_wav = start_wav
        self.end_wav = end_wav
        self.resolution = resolution
        
    def create_para_file(self, k_model_path, linelist_path, start_wav=15167.0, end_wav=16767.0, del_wav=0.02, smooth='g', smooth_width=0.75, atmosphere=1, lines=1):
        '''
        Function for creating the parameter file of batch.par
        '''
        MOOG_para_file = open(MOOGrun_path + '/batch.par', 'w')
        # Parameter list of MOOG: standard output file (1), summary output file (2), smoothed output file (3),
        #                         begin wavelength, end wavelength, wavelength step;
        #                         smoothing function, Gaussian FWHM, vsini, limb darkening coefficient,
        #                         Macrotrubulent FWHM, Lorentzian FWHM
        smooth_para = [smooth, smooth_width, 0.0, 0.0, 0.0, 0.0]
        #MOOG_para_file = open('batch.par', 'w')
        MOOG_contant = ["synth\n",
                        "atmosphere         {}\n".format(atmosphere),
                        "lines              {}\n".format(lines),
                        "standard_out       '{}'\n".format('MOOG.out1'),
                        "summary_out        '{}'\n".format('MOOG.out2'),
                        "smoothed_out       '{}'\n".format('MOOG.out3'),
                        "model_in           '{}'\n".format(k_model_path),
                        "lines_in           '{}'\n".format(linelist_path),
                        "terminal           'x11'\n",
                        "synlimits\n",
                        "  {}  {}  {}  5.0\n".format(start_wav, end_wav, del_wav),
                        "plot        3\n",
                        "plotpars    1\n",
                        "  0.0  0.0  0.0  0.0 \n",
                        "  0.0  0.0  0.0  0.0 \n",
                        "  '{}'  {}  {}  {}  {}  {}\n".format(*smooth_para)
                    ]
        MOOG_para_file.writelines(MOOG_contant)
        MOOG_para_file.close()
        
    def prepare_file(self, model_file=None, line_list=None):
        '''
        Prepare the model, linelist and control files for MOOG.
        Can either provide stellar parameters and wavelengths or provide file names.
        If fine name(s) provided, the files will be copied to working directory for calculation.  
        '''
        if model_file == None:
            # Model file is not specified, will download Kurucz model according to stellar parameters.
            pass
        else:
            # Model file is specified; record model file name and copy to working directory.
            subprocess.run(['cp', model_file, MOOGrun_path], encoding='UTF-8', stdout=subprocess.PIPE)
            self.model_file = model_file.split('/')[-1]

            
        if line_list == None:
            # Linelist file is not specified, will use built-in VALD linelist according to wavelength specification.
            
            self.line_list = 'vald_sub'
        else:
            # Linelist file is specified; record linelist file name and copy to working directory.
            subprocess.run(['cp', line_list, MOOGrun_path], encoding='UTF-8', stdout=subprocess.PIPE)
            self.line_list = line_list.split('/')[-1]
            
        # Create parameter file.
        self.create_para_file(self.model_file, self.line_list, start_wav=self.start_wav, end_wav=self.end_wav, )


    
    def moog_run(self, output=False):
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
        # if version == 'FEB2017':
        MOOG_run = subprocess.run([MOOG_path + '/MOOGSILENT'], stdout=subprocess.PIPE,
                                  cwd=MOOGrun_path)
        # elif version == 'NOV2019':
        #     MOOG_run = subprocess.run(['/home/mingjie/software/MOOG/NOV2019_modified/MOOGSILENT'], stdout=subprocess.PIPE)
        if output:
            MOOG_run = str(MOOG_run.stdout, encoding = "utf-8").split('\n')
            MOOG_output = []
            for i in MOOG_run:
                if len(i) > 12:
                    ansi_escape = re.compile('\x1b\[...H')
                    temp = ansi_escape.sub('', i)
                    ansi_escape = re.compile('\x1b\[....H')
                    temp = ansi_escape.sub('', temp)
                    ansi_escape = re.compile('\x1b\[H')
                    temp = ansi_escape.sub('', temp)
                    ansi_escape = re.compile('\x1b\[2J')
                    MOOG_output.append(ansi_escape.sub('', temp))
            for i in MOOG_output:
                print(i)

    def read_spectra(self, type='smooth'):
        '''
        Read the output spectra of MOOG.

        Parameters
        ----------
        type : str, default 'smooth'
            Decide the type of spectra to be read. 'smooth' will read the smoothed spectra (MOOG.out3), and 'standard' will read the un-smoothed one (MOOG.out2).

        Returns
        ---------
        wav : a numpy array
            An array of wavelength
        flux : a numpy array
            An array of flux
        '''
        if type == 'standard':
            models_file = open(MOOGrun_path+'MOOG.out2')
            models = models_file.readline()
            models = models_file.readline()
            models = models_file.read().split()
            models = [float(i) for i in models]
            x = range(round((models[1]-models[0])/models[2])+1)
            model_wav = []
            for i in x:
                model_wav.append(models[0] + models[2]*i)
            model_flux = np.array(models[4:])
            self.wav, self.flux =  np.array(model_wav), np.array(model_flux)
        elif type == 'smooth':
            models_file = open(MOOGrun_path+'MOOG.out3')
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