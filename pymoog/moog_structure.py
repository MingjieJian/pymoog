#!/usr/bin/python
from . import private
import subprocess
from . import model, line_data
import re
import numpy as np

MOOG_path = '{}/.pymoog/moog_nosm/moog_nosm_NOV2019/'.format(private.os.environ['HOME'])
MOOG_rundir_path = '{}/.pymoog/'.format(private.os.environ['HOME'])
MOOG_file_path = '{}/.pymoog/files/'.format(private.os.environ['HOME'])

valid_batch_pars = {
    'abfind': ['standard_out', 'summary_out', 'model_in', 'lines_in', 'atmosphere', 'lines', 'molecules', 'terminal'],
    'binary': ['bin_raw_out', 'bin_smo_out', 'atmosphere', 'lines', 'molecules', 'deltaradvel', 'lumratio', 'terminal', 'plot', 'plotpars', 'binary:RUN1', 'synlimits', 'binary:RUN2', 'synlimits'],
    'blends': ['standard_out', 'summary_out', 'model_in', 'lines_in', 'atmosphere', 'lines', 'molecules', 'terminal', 'blenlimits'],
    'cog': ['standard_out', 'summary_out', 'model_in', 'lines_in', 'atmosphere', 'lines', 'molecules', 'terminal', 'coglimits'],
    'doflux': ['standard_out', 'summary_out', 'model_in', 'lines_in', 'atmosphere', 'lines', 'molecules', 'terminal', 'synlimits'],
    'synth': ['standard_out', 'summary_out', 'smoothed_out', 'model_in', 'lines_in', 'atmosphere', 'lines', 'molecules', 'terminal', 'isotopes', 'synlimits', 'plot', 'plotpars'],
    'weedout':['standard_out', 'model_in', 'lines_in', 'keeplines_out', 'tosslines_out', 'atmosphere', 'lines', 'molecules', 'terminal'],
    'synpop': ['standard_out', 'summary_out', 'smoothed_out', 'table_in', 'table_out', 'lines_in', 'popsyn_out', 'atmosphere', 'lines', 'molecules', 'isotopes', 'terminal', 'synlimits', 'plot', 'plotpars']
}

batch_pars_default = {
    'standard_out':'MOOG.out1', 'summary_out':'MOOG.out2', 'smoothed_out':'MOOG.out3', 
    'keeplines_out':'keep.list', 'tosslines_out':'toss.list',
    'bin_raw_out':'bin_raw.out', 'bin_smo_out':'bin_smo_out',
    'model_in':'model.mod', 'lines_in':'line.list', 
    'atmosphere':1, 'lines':1, 'molecules':1, 
    'lumratio':1,
    'plot':3,
    'terminal':'x11',
    'isotopes':'     0     0', 
    'bin_raw_out':'bin_raw.out',
    'bin_smo_out':'bin_smo.out',
    'binary:RUN1':["RUN                1",
                   "standard_out       'MOOG1.out1'",
                   "summary_out        'MOOG1.out2'",
                   "smoothed_out       'MOOG1.out3'",
                   "model_in           'model1.mod'",
                   "lines_in           'line.list'"],
    'binary:RUN2':["RUN                2",
                   "standard_out       'MOOG2.out1'",
                   "summary_out        'MOOG2.out2'",
                   "smoothed_out       'MOOG2.out3'",
                   "model_in           'model2.mod'",
                   "lines_in           'line.list'"],
    'popsyn_out':'rawpop'
}

batch_pars_format = {
    'standard_out':'str', 'summary_out':'str', 'smoothed_out':'str', 
    'keeplines_out':'str', 'tosslines_out':'str',
    'model_in':'str', 'lines_in':'str', 
    'atmosphere':'int', 'lines':'int', 'molecules':'int', 
    'plot':'int',
    'terminal':'str',
    'synlimits':'synlimits',
    'plotpars':'plotpars',
    'blenlimits':'blenlimits',
    'coglimits':'coglimits',
    'isotopes':'isotopes',
    'bin_raw_out':'str', 'bin_smo_out':'str',
    'deltaradvel':'float', 'lumratio':'float',
    'binary:RUN1':'binary:RUN1', 
    'binary:RUN2':'binary:RUN2'
}

para_format = {
            'str':"{:19s}'{}'\n",
            'int':"{:19s}{}\n",
            'float':"{:19s}{:7.2f}\n",
            'synlimits':"{}\n  {:.2f}  {:.2f}  {}  {:.2f}\n",
            'plotpars':"{}    1\n  0.0  0.0  0.0  0.0 \n  0.0  0.0  0.0  0.0 \n  {}  {:.3f}  {:.3f}  {:.3f}  {:.3f}  {:.3f}\n",
            'isotopes':'{} {}\n', 
            'blenlimits':"{}\n    {}  {}  {:.1f}",
            'coglimits':"{}\n  {}  {}  {}  {}  0\n",
            'binary:RUN1':'{}\n'*6,
            'binary:RUN2':'{}\n'*6
        }

class moog_structure(object):

    def __init__(self, run_type, prefix=''):
        '''
        Initiate moog_sturcture. This is the class which set the basic functions of pymoog, including creating run_dir, prepare  
        '''

        if prefix == '':
            self.rundir_path = '{}/{}-{}-{}/'.format(MOOG_rundir_path, run_type, private.datetime.now().strftime("%H:%M:%S.%f"), ''.join(private.secrets.choice(private.string.ascii_uppercase + private.string.ascii_lowercase) for i in range(9)))
        else:
            self.rundir_path = '{}/{}-{}-{}-{}/'.format(MOOG_rundir_path, prefix, run_type, private.datetime.now().strftime("%H:%M:%S.%f"), ''.join(private.secrets.choice(private.string.ascii_uppercase + private.string.ascii_lowercase) for i in range(9)))

        private.subprocess.run(['mkdir', '-p', self.rundir_path])
        
        self.run_type = run_type

    def create_model_file(self, model_file=None, model_format='moog', abun_change=None, molecules_include=None, model_type='marcs', model_chem='st', model_geo='auto'):
        # Create model file.
        if model_file == None:
            # Model file is not specified, will use builtin model according to stellar parameters.
            self.model = model.interpolate_model(self.teff, self.logg, self.m_h, vmicro=self.vmicro, vmicro_mode=self.vmicro_mode, mass=self.mass, abun_change=abun_change, molecules_include=molecules_include, save_name=self.rundir_path + 'model.mod', model_type=model_type, chem=model_chem, geo=model_geo)
            self.model_file = 'model.mod'
            self.vmicro_model = self.model['vmicro_model']
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

    def create_line_list(self, **args):
        # General line list routine for other drivers.
        # For drivers synthy, binary

        if 'del_wav' not in args.keys():
            args['del_wav'] = 0.02
            self.del_wav = args['del_wav']
        if 'del_wav_opac' not in args.keys():
            args['del_wav_opac'] = 1.0
            self.del_wav_opac = args['del_wav']
        if args['del_wav'] < 0.001:
            raise ValueError('del_wav cannot be smaller than 0.001; the calculation and I/O precision is not enough.')
        smooth_width = np.mean([self.start_wav / self.resolution, self.end_wav / self.resolution])
        smooth_width_num = int(np.ceil(smooth_width / args['del_wav']))

        if isinstance(self.line_list_in, str):
            if self.line_list_in[-5:] != '.list':
                # Linelist file is not specified, use internal line list;
                line_list = line_data.read_linelist(self.line_list_in, loggf_cut=args['loggf_cut'])
                line_data.save_linelist(line_list, self.rundir_path + 'line.list', wav_start=self.start_wav-smooth_width_num*2*args['del_wav'], wav_end=self.end_wav+smooth_width_num*2*args['del_wav'])
                self.line_list_name = 'line.list'
                self.line_list = line_list
            elif self.line_list[-5:] == '.list':
                # Linelist file is specified; record linelist file name and copy to working directory.
                subprocess.run(['cp', self.line_list_in, self.rundir_path], encoding='UTF-8', stdout=subprocess.PIPE)
                self.line_list_name = self.line_list_in.split('/')[-1]
                args['lines_in'] = self.line_list_name
        elif isinstance(self.line_list_in, private.pd.DataFrame):
            line_data.save_linelist(self.line_list_in, self.rundir_path + 'line.list', wav_start=self.start_wav-smooth_width_num*2*args['del_wav'], wav_end=self.end_wav+smooth_width_num*2*args['del_wav'])
            self.line_list_name = 'line.list'
            self.line_list = self.line_list_in
        else:
            raise TypeError('Type of input linelist have to be either str or pandas.DataFrame.')

    def prepare_file(self, model_file=None, model_format='moog', loggf_cut=None, abun_change=None, molecules_include=None, model_type='marcs', model_chem='st', model_geo='auto', **args):
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
        isotopes : dict of pairs {int:float, ...}
            Isotope ratios, have to be a dict of pairs of isotope number and ratio values.
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
        if ('smooth_para' not in args.keys() or args['smooth_para'] is None) and self.run_type in ['synth', 'binary', 'synpop']:
            args['smooth_para'] = ['g', 0.0, 0.0, 0.0, 0.0, 0.0]

        # Create model file.
        if model_file == None:
            # Model file is not specified, will use builtin model according to stellar parameters.
            self.model = model.interpolate_model(self.teff, self.logg, self.m_h, vmicro=self.vmicro, vmicro_mode=self.vmicro_mode, mass=self.mass, abun_change=abun_change, molecules_include=molecules_include, save_name=self.rundir_path + 'model.mod', model_type=model_type, chem=model_chem, geo=model_geo)
            self.model_file = 'model.mod'
            self.vmicro_model = self.model['vmicro_model']
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
        if self.run_type == 'abfind':
            if isinstance(self.line_list_in, str):
                # Linelist file have to be specified; record linelist file name and copy to working directory.
                subprocess.run(['cp', self.line_list_in, self.rundir_path], encoding='UTF-8', stdout=subprocess.PIPE)
                self.line_list_name = self.line_list_in.split('/')[-1]
            elif isinstance(self.line_list_in, private.pd.DataFrame):
                line_data.save_linelist(self.line_list_in.sort_values('id'), self.rundir_path + 'line.list')
                self.line_list_name = 'line.list'
                self.line_list = self.line_list_in
            else:
                raise TypeError('Type of input linelist have to be either str or pandas.DataFrame.')
        elif self.run_type == 'blends':
            if self.line_list_in[-5:] != '.list':
                # Linelist file is not specified, use internal line list;
                line_list = line_data.read_linelist(self.line_list_in, loggf_cut=loggf_cut, mode='npy')
                
                # Input EW into the linelist
                line_list = line_list[(line_list['wavelength'] >= self.start_wav) & (line_list['wavelength'] <= self.end_wav)].reset_index(drop=True)
                line_list.loc[1:, 'wavelength'] = -line_list.loc[1:, 'wavelength']
                line_list['EW'] = private.np.nan
                line_list.loc[0, 'EW'] = self.EW
                
                line_data.save_linelist(line_list, self.rundir_path + 'line.list', negative=True)
                self.line_list_name = 'line.list'
                self.line_list = line_list
            elif self.line_list[-5:] == '.list':
                # Linelist file is specified; record linelist file name and copy to working directory.
                private.subprocess.run(['cp', self.line_list_in, self.rundir_path], encoding='UTF-8', stdout=private.subprocess.PIPE)
                self.line_list_name = self.line_list_in.split('/')[-1]
                # Input EW into the linelist
                line_list = line_data.read_linelist(self.rundir_path + self.line_list_in)
                line_list.loc[1:, 'wavelength'] = -line_list.loc[1:, 'wavelength']
                line_list['EW'] = private.np.nan
                line_list.loc[0, 'EW'] = self.EW
                line_data.save_linelist(line_list, self.rundir_path + 'line.list', negative=True)
                self.line_list_name = 'line.list'
                self.line_list = line_list
        elif self.run_type == 'cog':
            # Linelist file must be specified; record linelist file name and copy to working directory.
            if isinstance(self.line_list_in, str):
                private.subprocess.run(['cp', self.line_list_in, self.rundir_path+'/line.list'], encoding='UTF-8', stdout=private.subprocess.PIPE)
                self.line_list_name = 'line.list'
                self.line_list = line_data.read_linelist(self.line_list_in)
            elif isinstance(self.line_list_in, private.pd.DataFrame):
                line_data.save_linelist(self.line_list_in, self.rundir_path + 'line.list')
                self.line_list_name = 'line.list'
                self.line_list = self.line_list_in
            else:
                raise TypeError('Type of input linelist have to be either str or pandas.DataFrame.')
        elif self.run_type == 'doflux':
            pass
        else:
            # General line list routine for other drivers.
            if self.run_type in ['synth', 'doflux', 'binary']:
                if 'del_wav' not in args.keys():
                    args['del_wav'] = 0.02
                    self.del_wav = args['del_wav']
                if 'del_wav_opac' not in args.keys():
                    args['del_wav_opac'] = 1.0
                    self.del_wav_opac = args['del_wav_opac']
                if args['del_wav'] < 0.001:
                    raise ValueError('del_wav cannot be smaller than 0.001; the calculation and I/O precision is not enough.')
                smooth_width = np.mean([self.start_wav / self.resolution, self.end_wav / self.resolution])
                smooth_width_num = int(np.ceil(smooth_width / args['del_wav']))
            else:
                smooth_width = 0 
            if isinstance(self.line_list_in, str):
                if self.line_list_in[-5:] != '.list':
                    # Linelist file is not specified, use internal line list;
                    self.line_list = line_data.read_linelist(self.line_list_in, loggf_cut=loggf_cut)
                    line_data.save_linelist(self.line_list, self.rundir_path + 'line.list', wav_start=self.start_wav-smooth_width_num*2*args['del_wav'], wav_end=self.end_wav+smooth_width_num*2*args['del_wav'])
                    self.line_list_name = 'line.list'
                elif self.line_list_in[-5:] == '.list':
                    # Linelist file is specified; record linelist file name and copy to working directory.
                    subprocess.run(['cp', self.line_list_in, self.rundir_path], encoding='UTF-8', stdout=subprocess.PIPE)
                    self.line_list_name = self.line_list.split('/')[-1]
                    args['lines_in'] = self.line_list_name
                    self.line_list = None
            elif isinstance(self.line_list_in, private.pd.DataFrame):
                line_data.save_linelist(self.line_list_in, self.rundir_path + 'line.list', wav_start=self.start_wav-smooth_width_num*2*args['del_wav'], wav_end=self.end_wav+smooth_width_num*2*args['del_wav'])
                self.line_list_name = 'line.list'
            else:
                raise TypeError('Type of input linelist have to be either str or pandas.DataFrame.')
                
        # Create parameter file.
        if self.run_type == 'blends':
            args['blenlimits'] = [self.edge_width, self.step, self.ele]
        if self.run_type == 'cog':
            args['coglimits'] = [self.cog_low, self.cog_up, self.cog_step, self.lp_step]
        self.create_para_file(args=args)

        # Misc for some drivers
        if self.run_type == 'synth' and self.doflux_cont:
            self.doflux_cont()
        
    def create_para_file(self, args=None):
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

        if self.run_type in ['synth', 'binary']:
            if 'del_wav' not in args.keys():
                args['del_wav'] = 0.02
            if 'del_wav_opac' not in args.keys():
                args['del_wav_opac'] = 1.0
            smooth_width = np.mean([self.start_wav / self.resolution, self.end_wav / self.resolution])
            smooth_width_num = int(np.ceil(smooth_width / args['del_wav']))
            if args['smooth_para'][1] == 0:
                args['smooth_para'][1] = smooth_width
            
            args['plotpars'] = args['smooth_para']
            
            args['synlimits'] = [self.start_wav - smooth_width_num*2*args['del_wav'], self.end_wav + smooth_width_num*2*args['del_wav'], args['del_wav'], args['del_wav_opac']]
        elif self.run_type == 'doflux':
            if 'del_wav' not in args.keys():
                args['del_wav'] = 0.02
            if 'ddel_wav_opacel_wav' not in args.keys():
                args['del_wav_opac'] = 1.0
            smooth_width = np.mean([self.start_wav / self.resolution, self.end_wav / self.resolution])
            smooth_width_num = int(np.ceil(smooth_width / args['del_wav']))
            
            args['synlimits'] = [self.start_wav - smooth_width_num*2*args['del_wav'], self.end_wav + smooth_width_num*2*args['del_wav'], args['del_wav'], args['del_wav_opac']]

        MOOG_para_file = open(self.rundir_path + '/batch.par', 'w')

        MOOG_contant = ["{}\n".format(self.run_type)]

        # Convert the isotope to the list
        if 'isotopes' in args.keys():
            assert np.all([not(hasattr(ele, '__len__')) for ele in args['isotopes'].values()])
            isotopes_sub_content = []
            isotopes_sub_content += ['     {:2.0f}     1'.format(len(args['isotopes']))]
            isotopes_sub_content += [f'  {ele:9.5f}    {args["isotopes"][ele]}' if ele > 100 else f'  {ele:9.3f}    {args["isotopes"][ele]}' for ele in args['isotopes'].keys()]
            args['isotopes'] = '\n'.join(isotopes_sub_content)

        # Constract the parameter dict.
        for ele in valid_batch_pars[self.run_type]:
            if ele in args.keys():
                content = args[ele]
            elif ele in batch_pars_default.keys():
                content = batch_pars_default[ele]
            else:
                raise ValueError('Please provide {} as a argument in prepare_file().'.format(ele))

            if type(content) not in [str, int, float]:
                if ':' in ele:
                    # For the special paras whose name is defined in pymoog but not MOOG.
                    MOOG_contant.append(para_format[batch_pars_format[ele]].format(*content))
                else:
                    
                    MOOG_contant.append(para_format[batch_pars_format[ele]].format(ele, *content))
            else:
                MOOG_contant.append(para_format[batch_pars_format[ele]].format(ele, content))
            

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

        if self.run_type == 'weedout':
            MOOG_run = subprocess.run([MOOG_path + '/MOOGSILENT'], stdout=subprocess.PIPE, input=bytes('{}'.format(self.kappa_ratio), 'utf-8'), cwd=self.rundir_path)
        elif self.run_type == 'cog':
            MOOG_run = subprocess.run([MOOG_path + '/MOOGSILENT'], stdout=private.subprocess.PIPE, input=bytes('n', 'utf-8'), cwd=self.rundir_path)
        elif self.run_type == 'binary':
            MOOG_run = subprocess.run([MOOG_path + '/MOOGSILENT'], stdout=private.subprocess.PIPE, input=bytes('q', 'utf-8'), cwd=self.rundir_path)
        else:
            MOOG_run = subprocess.run([MOOG_path + '/MOOGSILENT'], stdout=subprocess.PIPE, cwd=self.rundir_path)
        
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

    def remove_rundir(self):
        private.subprocess.run(['rm', '-r', self.rundir_path])