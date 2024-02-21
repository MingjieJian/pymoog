#!/usr/bin/python
from . import moog_structure
from . import private
from . import line_data
from . import model
from . import rundir_num

MOOG_path = '{}/.pymoog/moog_nosm/moog_nosm_NOV2019/'.format(private.os.environ['HOME'])
# self.rundir_path = '{}/.pymoog/rundir/'.format(private.os.environ['HOME'])
MOOG_file_path = '{}/.pymoog/files/'.format(private.os.environ['HOME'])

class blends(moog_structure.moog_structure):
    def __init__(self, teff, logg, m_h, start_wav, end_wav, EW, ele, vmicro=2, mass=1, line_list='ges', prefix='', vmicro_mode='flexible', edge_width=0.5, step=0.005):
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
        super(blends, self).__init__('blends', prefix=prefix)
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
        self.vmicro_mode = vmicro_mode
        self.edge_width = edge_width
        self.step = step

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
            self.remove_rundir()