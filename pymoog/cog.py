#!/usr/bin/python
from . import private
from . import line_data
from . import model

MOOG_path = '{}/.pymoog/moog_nosm/moog_nosm_NOV2019/'.format(private.os.environ['HOME'])
MOOG_run_path = '{}/.pymoog/rundir/'.format(private.os.environ['HOME'])
MOOG_file_path = '{}/.pymoog/files/'.format(private.os.environ['HOME'])

class cog:
    def __init__(self, teff, logg, m_h, line_list, cog_low=-7.5, cog_up=-3.5, cog_step=0.05, lp_step=0):
        '''
        Initiate a cog Instance and read the parameters. 
        line_list muse be provided with suffix of '.list', and it can only contain one line.
        
        Parameters
        ----------
        teff : int
            The effective temperature of the model
        logg : float
            logg value of the model
        m_h : float
            [M/H] value (overall metallicity) of the model
        line_list : str
            The name of the linelist file.
        cog_low : float, default -7.5
            The log(W/lambda) lower limit for the cog. If larger than -5 then it will be set as -5 in MOOG.
        cog_up : float, default -3.5
            The log(W/lambda) upper limit for the cog. When larger than -3 step jump of cog will present.
        cog_step : float, default 0.05
            The log(W/lambda) step size for the cog.
        lp_step : float, default 0
            An explicit line profile wavelength step size, if desired.
        '''
        self.teff = teff
        self.logg = logg
        self.m_h = m_h
        self.line_list = line_list
        self.cog_low = cog_low
        self.cog_up = cog_up
        self.cog_step = cog_step
        self.lp_step = lp_step
        
    def prepare_file(self, model_file=None, model_type='moog', abun_change=None, molecules=None, atmosphere=1, lines=1):
        '''
        Prepare the model, linelist and control files for cog.
        Can either provide stellar parameters or provide file names.
        If fine name(s) provided, the files will be copied to working directory for calculation.  
        
        Parameters
        ----------
        model_file : str, optional
            The name of the model file. If not specified will use internal Kurucz model.
             
        model_type : str, optional
            The type of the model file. Default is "moog" (then no conversion of format will be done); can be "moog", "kurucz-atlas9" and "kurucz-atlas12". 
            
        abun_change : dict of pairs {int:float, ...}
            Abundance change, have to be a dict of pairs of atomic number and [X/Fe] values.
        '''
        private.subprocess.run(['rm', MOOG_run_path + 'batch.par'])
        private.subprocess.run(['rm', MOOG_run_path + 'model.mod'])
        private.subprocess.run(['rm', MOOG_run_path + 'line.list'])
        private.os.system('rm ' + MOOG_run_path + 'MOOG.out*')
        
        if model_file == None:
            # Model file is not specified, will download Kurucz model according to stellar parameters.
            model.interpolate_model(self.teff, self.logg, self.m_h, abun_change=abun_change, molecules=molecules)
            self.model_file = 'model.mod'
        else:
            # Model file is specified; record model file name and copy to working directory.
            if model_type == 'moog':
                private.subprocess.run(['cp', model_file, MOOG_run_path], encoding='UTF-8', stdout=private.subprocess.PIPE)
                self.model_file = model_file.split('/')[-1]
            elif model_type[:6] == 'kurucz':
                model.KURUCZ_convert(model_path=model_file, abun_change=abun_change, model_type=model_type[7:], molecules=molecules)
                self.model_file = 'model.mod'

        # Linelist file is specified; record linelist file name and copy to working directory.
        private.subprocess.run(['cp', self.line_list, MOOG_run_path], encoding='UTF-8', stdout=private.subprocess.PIPE)
        self.line_list = self.line_list.split('/')[-1]
            
        # Create parameter file.
        self.create_para_file(atmosphere=atmosphere, lines=lines)    
        
    def create_para_file(self, atmosphere=1, lines=1, molecules=1):
        '''
        Function for creating the parameter file of batch.par for cog.
        
        Parameters
        ----------
        
        Returns
        ----------
        None. A control file batch.par will be save in the pymoog working path.
        '''
        MOOG_para_file = open(MOOG_run_path + '/batch.par', 'w')
        # Parameter list of MOOG: standard output file (1), summary output file (2), smoothed output file (3),
        #                         begin wavelength, end wavelength, wavelength step;
        #                         smoothing function, Gaussian FWHM, vsini, limb darkening coefficient,
        #                         Macrotrubulent FWHM, Lorentzian FWHM
        #MOOG_para_file = open('batch.par', 'w')
        MOOG_contant = ["cog\n",
                        "standard_out       '{}'\n".format('MOOG.out1'),
                        "summary_out        '{}'\n".format('MOOG.out2'),
                        "model_in           '{}'\n".format(self.model_file),
                        "lines_in           '{}'\n".format(self.line_list),
                        "atmosphere         {}\n".format(atmosphere),
                        "lines              {}\n".format(lines),
                        "molecules          {}\n".format(molecules),
                        "terminal           'x11'\n",
                        "coglimits\n",
                        "  {}  {}  {}  {}  0\n".format(self.cog_low, self.cog_up, self.cog_step, self.lp_step)
                    ]
        MOOG_para_file.writelines(MOOG_contant)
        MOOG_para_file.close()
    
    def run_moog(self, output=False):
        '''
        Run MOOG and print the result if required.

        Parameters
        ----------
        output : boolen, default False
            If set to True, then print the out put of MOOG.

        Returns
        ----------
        None. Two files MOOG.out1 and MOOG.out2 will be save in the pymoog working path.
        '''
        
        MOOG_run = private.subprocess.run([MOOG_path + '/MOOGSILENT'], stdout=private.subprocess.PIPE, input=bytes('n', 'utf-8'), cwd=MOOG_run_path)

        
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

    def read_output(self):
        '''
        Read the output of cog.

        Parameters
        ----------

        Returns
        ---------
        self.loggf : a numpy array
            An array of loggf for cog.
        self.logrw : a numpy array
            An array of log(W/lambda) for cog.
        '''
        file = open(MOOG_run_path + 'MOOG.out2', 'r')
        cog_content = file.readlines()

        cog_single_content = [ele.replace(',', '').split() for ele in cog_content[6:]]
        cog_single_content = [item for sublist in cog_single_content for item in sublist]
        cog_array = private.np.array(cog_single_content, dtype=float).reshape(-1,2)
        loggf = cog_array[:,0]
        logrw = cog_array[:,1]
            
        self.loggf = loggf
        self.logrw = logrw