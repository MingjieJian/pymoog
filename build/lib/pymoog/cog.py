#!/usr/bin/python
from . import moog_structure
from . import private
from . import line_data
from . import model
from . import rundir_num

MOOG_path = '{}/.pymoog/moog_nosm/moog_nosm_NOV2019/'.format(private.os.environ['HOME'])
MOOG_file_path = '{}/.pymoog/files/'.format(private.os.environ['HOME'])

class cog(moog_structure.moog_structure):
    def __init__(self, teff, logg, m_h, line_list, vmicro=2., mass=1, cog_low=-7.5, cog_up=-3.5, cog_step=0.05, lp_step=0, prefix='', vmicro_mode='flexible'):
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
        mass : float, default 1
            The stellar mass of the input model. Only used when the model type is MARCS spherical.
        cog_low : float, default -7.5
            The log(W/lambda) lower limit for the cog. If larger than -6.5 then it will be set as -6.5 in MOOG.
        cog_up : float, default -3.5
            The log(W/lambda) upper limit for the cog. When larger than -3 step jump of cog will present.
        cog_step : float, default 0.05
            The log(W/lambda) step size for the cog.
        lp_step : float, default 0
            An explicit line profile wavelength step size, if desired.
        prefix : str, default ''.
            The prefix to be added to the name of rundir. Convenient when you want to find a specified rundir if there are many.
        edge_width : float, defalut 0.5
            The width to be included in the fitting around the central wavelength.
        step : float
            The wavelength step size of the syntheses.
        '''
        super(cog, self).__init__('cog', prefix=prefix)
        self.teff = teff
        self.logg = logg
        self.m_h = m_h
        self.line_list_in = line_list
        self.cog_low = cog_low
        self.cog_up = cog_up
        self.cog_step = cog_step
        self.lp_step = lp_step
        self.vmicro = vmicro
        self.mass = mass
        self.vmicro_mode = vmicro_mode

    def read_output(self, remove=True):
        '''
        Read the output of cog.

        Parameters
        ----------
        remove : bool, default True
            Whether remove the working folder after this function.

        Returns
        ---------
        self.loggf : a numpy array
            An array of loggf for cog.
        self.logrw : a numpy array
            An array of log(W/lambda) for cog.
        '''
        with open(self.rundir_path + 'MOOG.out2', 'r') as file:
            cog_content = file.readlines()

        res_df = private.pd.read_csv(self.rundir_path + 'MOOG.out2', skiprows=6, names=['loggf', 'logrw'])
        # cog_single_content = [ele.replace(',', '').split() for ele in cog_content[6:]]
        # cog_single_content = [item for sublist in cog_single_content for item in sublist]
        # cog_array = private.np.array(cog_single_content, dtype=float).reshape(-1,2)
            
        self.loggf = res_df['loggf'].values
        self.logrw = res_df['logrw'].values
        
        if remove:
            self.remove_rundir()
            
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
            input_loggf = self.line_list['loggf'].values[0]
        
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
            
            private.plt.axvline(sat_loggf_sep, ls='--', c='C3', label='linear-saturated separation')
            private.plt.axvline(dam_loggf_sep, ls='--', c='C4', label='saturated-damped separation')
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
                ax.axvline(sat_loggf_sep, ls='--', c='C3', label='linear-saturated separation')
                ax.axvline(dam_loggf_sep, ls='--', c='C4', label='saturated-damped separation')
                ax.axvline(input_loggf, c='C5', label='Input $\log{(gf)}$')

            ax1.set_title('Line status: {}'.format(line_status))
            ax1.legend(fontsize=7, loc=4)
            
        self.line_status = [line_status, sat_loggf_sep, dam_loggf_sep]