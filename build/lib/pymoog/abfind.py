#!/usr/bin/python
from . import moog_structure
from . import private

MOOG_path = '{}/.pymoog/moog_nosm/moog_nosm_NOV2019/'.format(private.os.environ['HOME'])
# self.rundir_path = '{}/.pymoog/rundir/'.format(private.os.environ['HOME'])
MOOG_file_path = '{}/.pymoog/files/'.format(private.os.environ['HOME'])

class abfind(moog_structure.moog_structure):
    def __init__(self, teff, logg, m_h, vmicro=2, mass=1, line_list='ges', prefix='', vmicro_mode='flexible'):
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
        super(abfind, self).__init__('abfind', prefix=prefix)
        self.teff = teff
        self.logg = logg
        self.m_h = m_h
        self.vmicro = vmicro
        self.mass = mass
        self.line_list_in = line_list
        self.vmicro_mode = vmicro_mode
    
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
            self.remove_rundir()
        
        self.abfind_res = abfind_dict
