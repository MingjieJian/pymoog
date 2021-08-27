#!/usr/bin/python
from . import private
from . import line_data
from . import model
from . import rundir_num

MOOG_path = '{}/.pymoog/moog_nosm/moog_nosm_NOV2019/'.format(private.os.environ['HOME'])
MOOG_file_path = '{}/.pymoog/files/'.format(private.os.environ['HOME'])

class cog(rundir_num.rundir_num):
    def __init__(self, teff, logg, m_h, line_list, vmicro=2., cog_low=-7.5, cog_up=-3.5, cog_step=0.05, lp_step=0, prefix=''):
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
        line_list : pandas.DataFrame
            The dataframe containting the line information.
        vmicro : float, default 2
            The microtrubulance velocity of the synthesized spectra (this is different from the v_micro of the atmosphere model which is always 2)
        cog_low : float, default -7.5
            The log(W/lambda) lower limit for the cog. If larger than -6.5 then it will be set as -6.5 in MOOG.
        cog_up : float, default -3.5
            The log(W/lambda) upper limit for the cog. When larger than -3 step jump of cog will present.
        cog_step : float, default 0.05
            The log(W/lambda) step size for the cog.
        lp_step : float, default 0
            An explicit line profile wavelength step size, if desired.
        '''
        super(cog, self).__init__('{}/.pymoog/'.format(private.os.environ['HOME']), 'cog', prefix=prefix)
        self.teff = teff
        self.logg = logg
        self.m_h = m_h
        self.line_list = line_list
        self.cog_low = cog_low
        self.cog_up = cog_up
        self.cog_step = cog_step
        self.lp_step = lp_step
        self.vmicro = vmicro
        
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
        
        if model_file == None:
            # Model file is not specified, will download Kurucz model according to stellar parameters.
            model.interpolate_model(self.teff, self.logg, self.m_h, abun_change=abun_change, molecules=molecules, to_path=self.rundir_path + 'model.mod', vmicro=self.vmicro)
            self.model_file = 'model.mod'
        else:
            # Model file is specified; record model file name and copy to working directory.
            if model_type == 'moog':
                private.subprocess.run(['cp', model_file, self.rundir_path], encoding='UTF-8', stdout=private.subprocess.PIPE)
                self.model_file = model_file.split('/')[-1]
            elif model_type[:6] == 'kurucz':
                model.KURUCZ_convert(model_path=model_file, abun_change=abun_change, model_type=model_type[7:], molecules=molecules, converted_model_path=self.rundir_path + 'model.mod')
                self.model_file = 'model.mod'

        # Linelist file must be specified; record linelist file name and copy to working directory.
        # self.line_list = line_data.read_linelist(self.line_list)
        # line_data.save_linelist(self.line_list, self.rundir_path + '/line.list')
        private.subprocess.run(['cp', self.line_list, self.rundir_path+'/line.list'], encoding='UTF-8', stdout=private.subprocess.PIPE)
        self.line_list_name = 'line.list'
        self.line_df = line_data.read_linelist(self.line_list)
            
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
        MOOG_para_file = open(self.rundir_path + '/batch.par', 'w')
        # Parameter list of MOOG: standard output file (1), summary output file (2), smoothed output file (3),
        #                         begin wavelength, end wavelength, wavelength step;
        #                         smoothing function, Gaussian FWHM, vsini, limb darkening coefficient,
        #                         Macrotrubulent FWHM, Lorentzian FWHM
        #MOOG_para_file = open('batch.par', 'w')
        MOOG_contant = ["cog\n",
                        "standard_out       '{}'\n".format('MOOG.out1'),
                        "summary_out        '{}'\n".format('MOOG.out2'),
                        "model_in           '{}'\n".format(self.model_file),
                        "lines_in           '{}'\n".format(self.line_list_name),
                        "atmosphere         {}\n".format(atmosphere),
                        "lines              {}\n".format(lines),
                        "molecules          {}\n".format(molecules),
                        "terminal           'x11'\n",
                        "coglimits\n",
                        "  {}  {}  {}  {}  0\n".format(self.cog_low, self.cog_up, self.cog_step, self.lp_step)
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
        None. Two files MOOG.out1 and MOOG.out2 will be save in the pymoog working path.
        '''
        
        MOOG_run = private.subprocess.run([MOOG_path + '/MOOGSILENT'], stdout=private.subprocess.PIPE, input=bytes('n', 'utf-8'), cwd=self.rundir_path)

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
        file = open(self.rundir_path + 'MOOG.out2', 'r')
        cog_content = file.readlines()

        cog_single_content = [ele.replace(',', '').split() for ele in cog_content[6:]]
        cog_single_content = [item for sublist in cog_single_content for item in sublist]
        cog_array = private.np.array(cog_single_content, dtype=float).reshape(-1,2)
        loggf = cog_array[:,0]
        logrw = cog_array[:,1]
            
        self.loggf = loggf
        self.logrw = logrw
        
        # if unlock:
        #     self.unlock()
            
    def get_line_status(self, input_loggf=private.np.nan, plot='none'):
        
        '''
        Derive spectral line status (Linear/Saturated/Damped) of the line input to cog using input loggf value.
        
        Parameters
        ----------
        input_loggf : float, default np.nan
            The input loggf for determining the line status. If not specified, will use the loggf of the line in cog.
        plot : str, default 'none'
            Plot status. If 'none', nothing will be plot; if 'cog', then cog will be plot; if 'all', then cog, the first and second derivative will be plot.    

        Returns
        ---------
        self.line_status : list
            List containing the line status: ['<Linear/Saturated/Damped>', <loggf separation between linear/saturated region>, <loggf separation between saturated/damped region>]
        '''
        
        # Check whether input_loggf is specified. If not then use line loggf
        if private.np.isnan(input_loggf):
            input_loggf = self.line_df['loggf'].values[0]
        
        logrw_poly = private.np.polyfit(self.loggf,self.logrw, 20)
        if private.np.any(private.np.abs(self.logrw - private.np.polyval(logrw_poly, self.loggf)) > 0.05):
            private.warnings.warn('Fitting to the cog have an discrepancy larger than 0.05 in some point, please check the cog in detail if necessary')
        logrw_1d = private.np.polyval(private.np.polyder(logrw_poly), self.loggf)
        logrw_2d = private.np.polyval(private.np.polyder(private.np.polyder(logrw_poly)), self.loggf)
        # Determine the length of logrw_2d to be set to 0.
        zero_length = int(private.np.ceil(2/151*len(self.loggf)))
        if zero_length < 2:
            zero_length = 2
        logrw_2d[:zero_length] = 0
        logrw_2d[-zero_length:] = 0

        sat_loggf_sep = self.loggf[private.np.argmin(logrw_2d)]
        dam_loggf_sep = self.loggf[private.np.argmax(logrw_2d)]
        sat_logrw_val = self.logrw[private.np.argmin(logrw_2d)]
        dam_logrw_val = self.logrw[private.np.argmax(logrw_2d)]

        # Check the values of seperations
        if (sat_loggf_sep > dam_loggf_sep and sat_logrw_val > dam_logrw_val) or logrw_1d[private.np.argmax(logrw_2d)] > 0.8:
            # Damping part not present
            dam_loggf_sep = private.np.nan
            dam_logrw_val = private.np.nan

        if logrw_1d[private.np.argmin(logrw_2d)] > 0.8:
            # Saturation part not present
            sat_loggf_sep = private.np.nan
            sat_logrw_val = private.np.nan

        # Derive the status of the line
        if input_loggf >= max(self.loggf) or input_loggf <= min(self.loggf):
            line_status = 'Unknown'
        else:
            if private.np.isnan(sat_loggf_sep) and private.np.isnan(dam_loggf_sep):
                if private.np.polyval(private.np.polyder(logrw_poly), input_loggf) >= 0.8:
                    line_status = 'Linear'
                else:
                    line_status = 'Unknown'
            elif private.np.isnan(dam_loggf_sep):
                if input_loggf <= sat_loggf_sep:
                    line_status = 'Linear'
                else:
                    line_status = 'Saturated'
            else:
                if input_loggf <= sat_loggf_sep:
                    line_status = 'Linear'
                elif input_loggf <= dam_loggf_sep:
                    line_status = 'Saturated'
                else:
                    line_status = 'Damped'

        if plot == 'cog':
            private.plt.figure(figsize=(6,3), dpi=100)
            private.plt.plot(self.loggf, self.logrw)
            private.plt.ylabel('$\log{(W/\lambda)}$')
            private.plt.xlabel("$\log{(gf)}$")
            
            private.plt.axvline(sat_loggf_sep, ls='--', c='C3', label='linear-stutared separation')
            private.plt.axvline(dam_loggf_sep, ls='--', c='C4', label='stutared-damped separation')
            private.plt.axvline(input_loggf, c='C5', label='Input $\log{(gf)}$')
            private.plt.legend(fontsize=7, loc=4)
            private.plt.title('Line status: {}'.format(line_status))
        elif plot == 'all':
            private.plt.figure(figsize=(6,3*3), dpi=100)

            ax1 = private.plt.subplot(311)
            private.plt.plot(self.loggf, self.logrw)
            private.plt.ylabel('$\log{(W/\lambda)}$')

            ax2 = private.plt.subplot(312)
            private.plt.plot(self.loggf, logrw_1d)
            private.plt.ylabel("$[\log{(W/\lambda)}]'$")

            ax3 = private.plt.subplot(313)
            private.plt.plot(self.loggf, logrw_2d)
            private.plt.ylabel("$[\log{(W/\lambda)}]''$")
            private.plt.xlabel("$\log{(gf)}$")

            for ax in [ax1, ax2, ax3]:
                ax.axvline(sat_loggf_sep, ls='--', c='C3', label='linear-stutared separation')
                ax.axvline(dam_loggf_sep, ls='--', c='C4', label='stutared-damped separation')
                ax.axvline(input_loggf, c='C5', label='Input $\log{(gf)}$')

            ax1.set_title('Line status: {}'.format(line_status))
            ax1.legend(fontsize=7, loc=4)
            
        self.line_status = [line_status, sat_loggf_sep, dam_loggf_sep]