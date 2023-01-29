#!/usr/bin/python
from . import private
from . import line_data
from . import model
from . import rundir_num
import subprocess

MOOG_path = '{}/.pymoog/moog_nosm/moog_nosm_NOV2019/'.format(private.os.environ['HOME'])
# self.rundir_path = '{}/.pymoog/rundir/'.format(private.os.environ['HOME'])
MOOG_file_path = '{}/.pymoog/files/'.format(private.os.environ['HOME'])

class abfind(rundir_num.rundir_num):
    def __init__(self, teff, logg, m_h, vmicro=2, mass=1, line_list='ges', prefix=''):
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
        prefix : str, default ''.
            The prefix to be added to the name of rundir. Convenient when you want to find a specified rundir if there are many.
        '''
        super(abfind, self).__init__('{}/.pymoog/'.format(private.os.environ['HOME']), 'abfind', prefix=prefix)
        self.teff = teff
        self.logg = logg
        self.m_h = m_h
        self.vmicro = vmicro
        self.mass = mass
        self.line_list = line_list
        
    def prepare_file(self, model_file=None, model_format='moog', abun_change=None, atmosphere=1, lines=1, molecules=1, molecules_include=None, model_type='marcs', model_chem='st', model_geo='auto'):
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


        # Create line list.
        if isinstance(self.line_list, str):
            # Linelist file have to be specified; record linelist file name and copy to working directory.
            subprocess.run(['cp', self.line_list, self.rundir_path], encoding='UTF-8', stdout=subprocess.PIPE)
            self.line_list = self.line_list.split('/')[-1]
        elif isinstance(self.line_list, private.pd.DataFrame):
            line_data.save_linelist(self.line_list.sort_values('id'), self.rundir_path + 'line.list')
            self.line_list = 'line.list'
        else:
            raise TypeError('Type of input linelist have to be either str or pandas.DataFrame.')
            
        # Create parameter file.
        self.create_para_file(atmosphere=atmosphere, lines=lines)    
        
    def create_para_file(self, atmosphere=1, lines=1, molecules=1):
        '''
        Function for creating the parameter file of batch.par for abfind.
        
        Parameters
        ----------
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

    def read_output(self, remove=True):
        '''
        Read the output of abfind.

        Parameters
        ----------
        remove : bool, default True
            Whether remove the working folder after this function.

        Returns
        ---------
        abfind_dict : dict
            A dictionary which the keys are the element index and values are DataFrames of abfind result.
        '''
        # file = open(self.rundir_path + 'MOOG.out2', 'r')

        with open(self.rundir_path + 'MOOG.out2', 'r') as file:
            abfind_content = file.readlines()
        para = private.np.array(private.re.findall('[0-9]+.[0-9]+', abfind_content[2]), dtype=float)

        sep_index = []
        for i in range(len(abfind_content)):
            if 'Abundance Results' in abfind_content[i]:
                sep_index.append(i)
        sep_index.append(len(abfind_content)+1)
                
        abfind_dict = {}
        for i in range(len(sep_index)-1):
            abfind_s_df = private.pd.read_csv(self.rundir_path + 'MOOG.out2', 
                                    skiprows=sep_index[i]+1, skipfooter=len(abfind_content)-1-(sep_index[i+1]-6),
                                    engine='python', sep=' +')
            ele_index = abfind_s_df.loc[0, 'ID']
            abfind_dict[ele_index] = abfind_s_df
            
        if remove:
            self.remove()
        
        self.abfind_res = abfind_dict
