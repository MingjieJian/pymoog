#!/usr/bin/python
from . import private
from . import line_data
from . import model
from . import rundir_num

MOOG_path = '{}/.pymoog/moog_nosm/moog_nosm_NOV2019/'.format(private.os.environ['HOME'])
# self.rundir_path = '{}/.pymoog/rundir/'.format(private.os.environ['HOME'])
MOOG_file_path = '{}/.pymoog/files/'.format(private.os.environ['HOME'])

class blends(rundir_num.rundir_num):
    def __init__(self, teff, logg, m_h, start_wav, end_wav, EW, ele, line_list='ges', prefix=''):
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
        start_wav : float
            The start wavelength of the blended feature.
        end_wav : float
            The end wavelength of the blended feature.
        EW : float
            The measured equivalent width.
        ele : float
            The element index of the dominant line in the feature, e.g., Fe I -> 26.0.
        line_list : str, default 'ges'
            The name of the linelist file. If not specified will use built-in VALD linelist.
        '''
        super(blends, self).__init__('{}/.pymoog/'.format(private.os.environ['HOME']), 'blends', prefix=prefix)
        self.teff = teff
        self.logg = logg
        self.m_h = m_h
        self.start_wav = start_wav
        self.end_wav = end_wav
        self.EW = EW
        self.ele = ele
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

        if self.line_list[-5:] != '.list':
            # Linelist file is not specified, use internal line list;
            line_list = line_data.read_linelist(self.line_list, loggf_cut=loggf_cut, mode='npy')
            
            # Input EW into the linelist
            line_list = line_list[(line_list['wavelength'] >= self.start_wav) & (line_list['wavelength'] <= self.end_wav)].reset_index(drop=True)
            line_list.loc[1:, 'wavelength'] = -line_list.loc[1:, 'wavelength']
            line_list['EW'] = private.np.nan
            line_list.loc[0, 'EW'] = self.EW
            
            line_data.save_linelist(line_list, self.rundir_path + 'line.list', negative=True)
            self.line_list = 'line.list'
        elif self.line_list[-5:] == '.list':
            # Linelist file is specified; record linelist file name and copy to working directory.
            private.subprocess.run(['cp', self.line_list, self.rundir_path], encoding='UTF-8', stdout=private.subprocess.PIPE)
            self.line_list = self.line_list.split('/')[-1]
            # Input EW into the linelist
            line_list = line_data.read_linelist(self.rundir_path + self.line_list)
            line_list.loc[1:, 'wavelength'] = -line_list.loc[1:, 'wavelength']
            line_list['EW'] = private.np.nan
            line_list.loc[0, 'EW'] = self.EW
            line_data.save_linelist(line_list, self.rundir_path + 'line.list', negative=True)
            
        # Create parameter file.
        self.create_para_file(self.ele, atmosphere=atmosphere, lines=lines)    
        
    def create_para_file(self, ele, atmosphere=1, lines=1, molecules=1, edge_width=0.5, step=0.005):
        '''
        Function for creating the parameter file of batch.par for blends.
        
        Parameters
        ----------
        
        '''
        MOOG_para_file = open(self.rundir_path + '/batch.par', 'w')
        # Parameter list of MOOG: standard output file (1), summary output file (2), smoothed output file (3),
        #                         begin wavelength, end wavelength, wavelength step;
        #                         smoothing function, Gaussian FWHM, vsini, limb darkening coefficient,
        #                         Macrotrubulent FWHM, Lorentzian FWHM
        #MOOG_para_file = open('batch.par', 'w')
        MOOG_contant = ["blends\n",
                        "standard_out       '{}'\n".format('MOOG.out1'),
                        "summary_out        '{}'\n".format('MOOG.out2'),
                        "model_in           '{}'\n".format(self.model_file),
                        "lines_in           '{}'\n".format(self.line_list),
                        "atmosphere         {}\n".format(atmosphere),
                        "lines              {}\n".format(lines),
                        "molecules          {}\n".format(molecules),
                        "terminal           'x11'\n",
                        "blenlimits\n",
                        "    {}  {}  {:.1f}".format(edge_width, step, ele)
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
                
        # if unlock:
        #     self.unlock()
                
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
        self.blends_s_df : pandas DataFrame
            A pandas DataFrame containting one-line result of blends.
        '''
        file = open(self.rundir_path + 'MOOG.out2', 'r')
        blends_content = file.readlines()
        para = private.np.array(private.re.findall('[0-9]+.[0-9]+', blends_content[1]), dtype=float)

        sep_index = []
        for i in range(len(blends_content)):
            if 'wavelength' in blends_content[i]:
                begin_index = i
            elif 'average abundance' in blends_content[i]:
                end_index = i
                print(begin_index, end_index)
                break

        blends_s_df = private.pd.DataFrame(private.np.array([ele.split() for ele in blends_content[begin_index+1:end_index]], dtype=float), columns=blends_content[begin_index].split())
        # Exclude the lines with no measurement.
        blends_s_df = blends_s_df[blends_s_df['abund'] != 999.99].reset_index(drop=True)

        self.blends_s_df = blends_s_df

        # if unlock:
        #     self.unlock()