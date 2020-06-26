#!/usr/bin/python
import os
from PyAstronomy import pyasl
import subprocess
import numpy as np
import re
import line_data

MOOG_path = '{}/.pymoog/moog_nosm/moog_nosm_FEB2017/'.format(os.environ['HOME'])
MOOGrun_path = '{}/.pymoog/rundir/'.format(os.environ['HOME'])

class synth:
    def __init__(self, teff, logg, m_h, start_wav, end_wav, resolution, del_wav=0.02, smooth='g'):
        '''
        Initiate a synth Instance and read the parameters.
        Parameters
        ----------
        teff : int
            The effective temperature of the model
        logg : float
            logg value of the model
        m_h : float
            [M/H] value (overall metallicity) of the model
        start_wav : float
            The start wavelength of synthetic spenctra
        end_wav : float
            The end wavelength of synthetic spenctra
        resolution : float
            Resolution of the synthetic spectra; this will passed to MOOG and convolute with initial spectra.
        '''
        self.teff = teff
        self.logg = logg
        self.m_h = m_h
        self.start_wav = start_wav
        self.end_wav = end_wav
        self.resolution = resolution
      
    def KURUCZ_download(self, save_name=None):
        '''
        Download the Kurucz ATLAS9 model using pyasl.

        Parameters
        ----------
        teff : int
            The effective temperature of the model
        logg : float
            logg value of the model
        m_h : float
            [M/H] value (overall metallicity) of the model
        save_name : str
            The path to save the model file (including the name).

        Returns
        ----------
        self.model_path : str
            The path and name of saved model.

        '''

        #The model is invoked here
        model = pyasl.getKuruczModel(self.teff, self.logg, self.m_h)

        if save_name == None:
            output_file_name = MOOGrun_path + '/model.mod'
        else:
            output_file_name = save_name
        output_file = open(output_file_name,"w")

        for line in model:
            output_file.write(line + "\n")

        output_file.close()

        self.model_path = output_file_name


    def KURUCZ_convert(self, model_path=None, vmicro=2.0, abun_change=None, converted_model_path=None):
        '''
        Convert the model file of ATLAS9-APOGEE in to MOOG format.

        Parameters
        ----------
        model_path : str, optional
            The path of donloaded model file.
        v_micro : float, default 2.0
            Microturbulance velocity of the spectra.
        abun_change : list of pairs [int, float]
            Abundance change, have to be a list of pairs of atomic number and [X/Fe] values.
        converted_model_path : str, optional
            The name of converted model. Default will be saved into MOOG working folder.
        '''

        if model_path == None:
            model_file = open(self.model_path)
        else:
            model_file = open(model_path)
        # Convert the model files into MOOG format.

        # Read and save the first two lines (except 'TITLE ') into header.
        header = model_file.readline()[:-2] + model_file.readline()[6:-2].strip()
        m_h = re.findall(r'\[(.*)\]', header)[0]

        # Read the abundance change as well as model lines.
        temp = model_file.readline() + model_file.readline() + model_file.readline()

        abun_list = ''
        temp = model_file.readline()
        while temp[:17] == ' ABUNDANCE CHANGE':
            abun_list = abun_list + temp[17:]
            temp = model_file.readline()
        abun = np.array(abun_list.split(), dtype='f').reshape(int(len(abun_list.split())/2), 2)

        # Load the abundances from Asplund 2009 (which MOOG used; hard-coded 20180531)
        xabu = [12.00,10.93, 1.05, 1.38, 2.70, 8.43, 7.83, 8.69, 4.56, 7.93,
                6.24, 7.60, 6.45, 7.51, 5.41, 7.12, 5.50, 6.40, 5.03, 6.34,
                3.15, 4.95, 3.93, 5.64, 5.43, 7.50, 4.99, 6.22, 4.19, 4.56,
                3.04, 3.65, 2.30, 3.34, 2.54, 3.25, 2.52, 2.87, 2.21, 2.58,
                1.46, 1.88,-5.00, 1.75, 0.91, 1.57, 0.94, 1.71, 0.80, 2.04,
                1.01, 2.18, 1.55, 2.24, 1.08, 2.18, 1.10, 1.58, 0.72, 1.42,
                -5.00, 0.96, 0.52, 1.07, 0.30, 1.10, 0.48, 0.92, 0.10, 0.84,
                0.10, 0.85,-0.12, 0.85, 0.26, 1.40, 1.38, 1.62, 0.92, 1.17,
                0.90, 1.75, 0.65,-5.00,-5.00,-5.00,-5.00,-5.00,-5.00, 0.02,
                -5.00,-0.54,-5.00,-5.00,-5.00]

        # Read the model lines
        temp = temp.split()
        model_lines = []
        model_linen = int(temp[2])
        for i in range(int(temp[2])):
            model_lines.append(model_file.readline().split()[:7])

        # Prepare the microtrubulance value.
        vmicro = '{}E00'.format(vmicro)

        # Write the model file.
        # header, abun09, model_lines and vmicro
        if converted_model_path == None:
            c_model_path = MOOGrun_path + '/model.mod'
        else:
            c_model_path = converted_model_path
        c_model_file = open(c_model_path, 'w')

        # Header part
        c_model_file.writelines('KURUCZ\n')
        c_model_file.writelines(header + '\n')

        # Model part
        c_model_file.writelines('ntau=       ' + str(model_linen) + '\n')
        for i in model_lines:
            c_model_file.writelines(' ' + ' '.join(i) + '\n')

        # Microturbulant velocity part
        c_model_file.writelines('    ' + vmicro + '\n')

        # Element shift part
        if abun_change != None:
            abun_change_num = len(abun_change)
            c_model_file.writelines('NATOMS      {}   {}\n'.format(abun_change_num, m_h))
            for abun in abun_change:
                c_model_file.writelines('      {:.2f}    {:.2f}\n'.format(abun[0], xabu[abun[0]-1]+abun[1]+float(m_h)))
        else:
            c_model_file.writelines('NATOMS      0   {}\n'.format(m_h))
    #         for i in abun09:
    #             k_model_file.writelines('          {:2g} {:.3f}\n'.format(i[0], i[1]))

        # Molecular line switches (temporary closed)
        c_model_file.writelines('NMOL        0')
        c_model_file.close()

        # cv_situation = os.path.isfile(c_model_path)
        # return c_model_file, cv_situation  
        
    def prepare_file(self, model_file=None, line_list=None, loggf_cut=None):
        '''
        Prepare the model, linelist and control files for MOOG.
        Can either provide stellar parameters and wavelengths or provide file names.
        If fine name(s) provided, the files will be copied to working directory for calculation.  
        
        Parameters
        ----------
        model_file : str, optional
            The name of the model file. If not specified will download Kurucz model. 
        line_list : str
            The name of the linelist file. If not specified will use built-in VALD linelist.
        logf_cut : float, optional
            The cut in loggf; if specified will only include the lines with loggf >= loggf_cut.
        '''
        if model_file == None:
            # Model file is not specified, will download Kurucz model according to stellar parameters.
            self.KURUCZ_download()
            self.KURUCZ_convert()
            self.model_file = 'model.mod'
        else:
            # Model file is specified; record model file name and copy to working directory.
            subprocess.run(['cp', model_file, MOOGrun_path], encoding='UTF-8', stdout=subprocess.PIPE)
            self.model_file = model_file.split('/')[-1]

            
        if line_list == None:
            # Linelist file is not specified, will use built-in VALD linelist according to wavelength specification.
            vald = line_data.read_linelist('files/linelist/vald', loggf_cut=loggf_cut)
            line_data.save_linelist(vald, MOOGrun_path + 'vald_sub', wav_start=self.start_wav, wav_end=self.end_wav)
            self.line_list = 'vald_sub'
        else:
            # Linelist file is specified; record linelist file name and copy to working directory.
            subprocess.run(['cp', line_list, MOOGrun_path], encoding='UTF-8', stdout=subprocess.PIPE)
            self.line_list = line_list.split('/')[-1]
            
        # Create parameter file.
        self.create_para_file()    
        
    def create_para_file(self, del_wav=0.02, smooth='g', atmosphere=1, lines=1):
        '''
        Function for creating the parameter file of batch.par
        
        Parameters
        ----------
        del_wav : float, optional
            The sampling distance of synthetic spectra. Default 0.02.
        smooth : str, optional
            Line profile to be used to smooth the synthetic spectra, same as decribed in MOOG documention. Default Gaussian. 
        '''
        MOOG_para_file = open(MOOGrun_path + '/batch.par', 'w')
        # Parameter list of MOOG: standard output file (1), summary output file (2), smoothed output file (3),
        #                         begin wavelength, end wavelength, wavelength step;
        #                         smoothing function, Gaussian FWHM, vsini, limb darkening coefficient,
        #                         Macrotrubulent FWHM, Lorentzian FWHM
        smooth_width = np.mean([self.start_wav / self.resolution, self.end_wav / self.resolution])
        smooth_para = [smooth, smooth_width, 0.0, 0.0, 0.0, 0.0]
        #MOOG_para_file = open('batch.par', 'w')
        MOOG_contant = ["synth\n",
                        "atmosphere         {}\n".format(atmosphere),
                        "lines              {}\n".format(lines),
                        "standard_out       '{}'\n".format('MOOG.out1'),
                        "summary_out        '{}'\n".format('MOOG.out2'),
                        "smoothed_out       '{}'\n".format('MOOG.out3'),
                        "model_in           '{}'\n".format(self.model_file),
                        "lines_in           '{}'\n".format(self.line_list),
                        "terminal           'x11'\n",
                        "synlimits\n",
                        "  {}  {}  {}  4.0\n".format(self.start_wav, self.end_wav, del_wav),
                        "plot        3\n",
                        "plotpars    1\n",
                        "  0.0  0.0  0.0  0.0 \n",
                        "  0.0  0.0  0.0  0.0 \n",
                        "  '{}'  {:.3f}  {}  {}  {}  {}\n".format(*smooth_para)
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
                    ansi_escape = re.compile(r'\x1b\[...H')
                    temp = ansi_escape.sub('', i)
                    ansi_escape = re.compile(r'\x1b\[....H')
                    temp = ansi_escape.sub('', temp)
                    ansi_escape = re.compile(r'\x1b\[H')
                    temp = ansi_escape.sub('', temp)
                    ansi_escape = re.compile(r'\x1b\[2J')
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