#!/usr/bin/python
from . import private
from . import line_data
from . import model
from . import rundir_num

MOOG_path = '{}/.pymoog/moog_nosm/moog_nosm_NOV2019/'.format(private.os.environ['HOME'])
# self.rundir_path = '{}/.pymoog/rundir/'.format(private.os.environ['HOME'])
MOOG_file_path = '{}/.pymoog/files/'.format(private.os.environ['HOME'])

class abfind(rundir_num.rundir_num):
    def __init__(self, teff, logg, m_h, line_list='ges', prefix=''):
        '''
        Initiate a abfind Instance and read the parameters.
        
        Parameters
        ----------
        teff : int
            The effective temperature of the model
        logg : float
            logg value of the model
        m_h : float
            [M/H] value (overall metallicity) of the model
        line_list : str, default 'ges'
            The name of the linelist file. If not specified will use built-in VALD linelist.
        '''
        super(abfind, self).__init__('{}/.pymoog/'.format(private.os.environ['HOME']), 'abfind', prefix=prefix)
        self.teff = teff
        self.logg = logg
        self.m_h = m_h
        self.line_list = line_list
        
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
                private.subprocess.run(['cp', model_file, self.rundir_path], encoding='UTF-8', stdout=private.subprocess.PIPE)
                self.model_file = model_file.split('/')[-1]
            elif model_type[:6] == 'kurucz':
                model.KURUCZ_convert(model_path=model_file, abun_change=abun_change, model_type=model_type[7:], molecules=molecules, converted_model_path=self.rundir_path + 'model.mod')
                self.model_file = 'model.mod'

        # Linelist file is specified; record linelist file name and copy to working directory.
        private.subprocess.run(['cp', self.line_list, self.rundir_path], encoding='UTF-8', stdout=private.subprocess.PIPE)
        self.line_list = self.line_list.split('/')[-1]
            
        # Create parameter file.
        self.create_para_file(atmosphere=atmosphere, lines=lines)    
        
    def create_para_file(self, atmosphere=1, lines=1, molecules=1):
        '''
        Function for creating the parameter file of batch.par for abfind.
        
        Parameters
        ----------
        
        '''
        MOOG_para_file = open(self.rundir_path + '/batch.par', 'w')
        # Parameter list of MOOG: standard output file (1), summary output file (2), smoothed output file (3),
        #                         begin wavelength, end wavelength, wavelength step;
        #                         smoothing function, Gaussian FWHM, vsini, limb darkening coefficient,
        #                         Macrotrubulent FWHM, Lorentzian FWHM
        #MOOG_para_file = open('batch.par', 'w')
        MOOG_contant = ["abfind\n",
                        "standard_out       '{}'\n".format('MOOG.out1'),
                        "summary_out        '{}'\n".format('MOOG.out2'),
                        "model_in           '{}'\n".format(self.model_file),
                        "lines_in           '{}'\n".format(self.line_list),
                        "atmosphere         {}\n".format(atmosphere),
                        "lines              {}\n".format(lines),
                        "molecules          {}\n".format(molecules),
                        "terminal           'x11'\n",
                    ]
        MOOG_para_file.writelines(MOOG_contant)
        MOOG_para_file.close()
    
    def run_moog(self, output=False, unlock=False):
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
        
        MOOG_run = private.subprocess.run([MOOG_path + '/MOOGSILENT'], stdout=private.subprocess.PIPE, cwd=self.rundir_path)

        # if unlock:
        #     self.unlock()
        
        MOOG_run = str(MOOG_run.stdout, encoding = "utf-8").split('\n')
        MOOG_output = []
        for i in MOOG_run:
            if len(i) > 12:
                ansi_escape = private.re.compile(r'\x1b\[...H')
                temp = ansi_escape.sub('', i)
                ansi_escape = private.re.compile(r'\x1b\[....H')
                temp = ansi_escape.sub('', temp)
                ansi_escape = private.re.compile(r'\x1b\[H')
                temp = ansi_escape.sub('', temp)
                ansi_escape = private.re.compile(r'\x1b\[2J')
                MOOG_output.append(ansi_escape.sub('', temp))
                
        if output:
            for i in MOOG_output:
                print(i)
        
        if 'ERROR' in ''.join(MOOG_run):
            raise ValueError('There is error during the running of MOOG.')

    def read_output(self, unlock=True):
        '''
        Read the output of abfind.

        Parameters
        ----------
        None.

        Returns
        ---------
        abfind_dict : dict
            A dictionary which the keys are the element index and values are DataFrames of abfind result.
        '''
        file = open(self.rundir_path + 'MOOG.out2', 'r')
        abfind_content = file.readlines()
        para = private.np.array(private.re.findall('[0-9]+.[0-9]+', abfind_content[2]), dtype=float)

        sep_index = []
        for i in range(len(abfind_content)):
            if 'Abundance Results' in abfind_content[i]:
                sep_index.append(i)

        abfind_dict = {}
        for i in range(len(sep_index)-1):
            abfind_single = abfind_content[sep_index[i]:sep_index[i+1]]
            for j in range(len(abfind_single)):
                if 'average abundance' in abfind_single[j]:
                    end_index = j

            abfind_s_df = private.pd.DataFrame(private.np.array([ele.split() for ele in abfind_single[2:end_index]], dtype=float), columns=abfind_single[1].split())

            ele_index = abfind_s_df.loc[0, 'ID']
            abfind_dict[ele_index] = abfind_s_df
            
        # if unlock:
        #     self.unlock()
        
        return abfind_dict