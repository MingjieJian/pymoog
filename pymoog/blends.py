#!/usr/bin/python
from . import private
from . import line_data
from . import model
from . import rundir_num

MOOG_path = '{}/.pymoog/moog_nosm/moog_nosm_NOV2019/'.format(private.os.environ['HOME'])
# self.rundir_path = '{}/.pymoog/rundir/'.format(private.os.environ['HOME'])
MOOG_file_path = '{}/.pymoog/files/'.format(private.os.environ['HOME'])

class blends(rundir_num.rundir_num):
    def __init__(self, teff, logg, m_h, start_wav, end_wav, EW, ele, vmicro=2, mass=1, line_list='ges', prefix=''):
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
        vmicro : float, default 2
            The microturbulance velocity of the model. 
        mass : float, default 1
            The stellar mass of the input model. Only used when the model type is MARCS spherical.
        line_list : str, default 'ges'
            The name of the linelist file. If not specified will use built-in VALD linelist.
        prefix : str, default ''.
            The prefix to be added to the name of rundir. Convenient when you want to find a specified rundir if there are many.
        '''
        super(blends, self).__init__('{}/.pymoog/'.format(private.os.environ['HOME']), 'blends', prefix=prefix)
        self.teff = teff
        self.logg = logg
        self.m_h = m_h
        self.vmicro = vmicro
        self.mass = mass
        self.start_wav = start_wav
        self.end_wav = end_wav
        self.EW = EW
        self.ele = ele
        self.line_list = line_list
        
    def prepare_file(self, model_file=None, model_format='moog', loggf_cut=None, abun_change=None, atmosphere=1, lines=1, molecules=1, molecules_include=None, model_type='marcs', model_chem='st', model_geo='auto'):
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
        
        if model_file == None:
            # Model file is not specified, will download Kurucz model according to stellar parameters.
            model.interpolate_model(self.teff, self.logg, self.m_h, vmicro=self.vmicro, mass=self.mass, abun_change=abun_change, molecules_include=molecules_include, save_name=self.rundir_path + 'model.mod', model_type=model_type, chem=model_chem, geo=model_geo)
            self.model_file = 'model.mod'
        else:
            # Model file is specified; record model file name and copy to working directory.
            if model_format == 'moog':
                private.subprocess.run(['cp', model_file, self.rundir_path], encoding='UTF-8', stdout=private.subprocess.PIPE)
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
        ele : float
            Element of the target line.
        atmosphere : int, default 1
            The atmosphere value described in MOOG documention, section III.
        lines : int, default 1
            The lines value described in MOOG documention, section III.
        molecules : int, default 1
            The molecules value described in MOOG documention, section III.
        edge_width : float, defalut 0.5
            The width to be included in the fitting around the central wavelength.
        step : float
            The wavelength step size of the syntheses.
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

    def read_output(self, remove=True):
        '''
        Read the output of abfind.

        Parameters
        ----------
        remove : bool, default True
            Whether remove the working folder after this function.

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
                # print(begin_index, end_index)
                break

        blends_s_df = private.pd.DataFrame(private.np.array([ele.split() for ele in blends_content[begin_index+1:end_index]], dtype=float), columns=blends_content[begin_index].split())
        # Exclude the lines with no measurement.
        blends_s_df = blends_s_df[blends_s_df['abund'] != 999.99].reset_index(drop=True)

        self.blends_s_df = blends_s_df

        if remove:
            self.remove()