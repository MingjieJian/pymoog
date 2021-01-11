# Keywords

## freeform 
indicates whether or not the linelist will be read in under format control (7e10.3) or will be free-form.  If freeform = 0, the default value, then the old-style formatted input will be used; If freeform = 1, unformatted read will be used, BUT the user must then give values for all quantities (that is, explicit zeros will need to be put instead of blank spaces.

## RUN
signals that there are either multiple syntheses being done or multiple comparisons with observed spectra.

## standard_out
controls the name of the verbose standard output.

## summary_out
controls the name of either the EW summary or the raw synthesis output.

## hardpost_out
controls the name of a postscript plot output.

## speccomp_out
controls the name of a text file containing the comparisons (wavelength shifts, sigmas, etc.) between observed and synthetic spectra

## bin_raw_out
controls the name of a file containing the raw synthesis of a spectroscopic binary, with an appropriate velocity  difference and luminosity ratio dialed in

## bin_smo_out
controls the name of a file containing the smoothed synthesis of a spectroscopic binary

## summary_in
controls the name of the raw synthesis file, created previously, that will be read in for plotting purposes

## smoothed_out
controls the name of the smoothed synthesis output

## keeplines_out
controls the name of the list of kept lines for future synthetic spectrum runs

## tosslines_out
controls the name of the list of discarded lines that are too weak to keep in future synthetic spectrum runs

## iraf_out
controls the name of the optional IRAF output

## model_in
controls the name of input model atmosphere file

## lines_in
controls the name of the input line list

## stronglines_in
controls the name of the input strong line list

## observed_in
controls the name of the input observed spectrum

## table_in
controls the name of the extra input instruction file

## table_out
controls the name of the extra input instruction file

## popsyn_out
controls the name of the extra input instruction file

## rawbin_out 
controls the name of the input observed spectrum

## smoobin_out
controls the name of the input observed spectrum

## atmosphere
controls the output of atmosphere quantities:
- 0 = do not print out the atmosphere
- 1 = print out the standard things about an atmsophere
- 2 = print standard things and additional stuff like continuous opacities, etc.

## molecules 
controls the molecular equilibrium calculations:
- 0 = do not do molecular equilibrium
- 1 = do molecular equilibrium but do not print results
- 2 = do molecular equilibrium and print results

## molset
controls the choice of which set of molecules will be used in molecular equilibrium calculations.
- 1 = the small set involving H, C, N, O, Mg, Ti (DEFAULT)
- 2 = the large set more useful for very cool stars

## deviations
controls whether, for synthetic spectrum computations, an 'obs-comp
plot will be made in addition to the normal spectrum plot
- 0 = do not plot the obs-comp plot
- 1 = plot the obs-comp plot

## lines     
controls the output of line data
- 0 = print out nothing about the input lines
- 1 = print out standard information about the input line list
- 2 = gory line data print (usually for diagnostic purposes)

## gfstyle   
controls the output of line data
- 0 = base-10 logarithms of the gf values (DEFAULT)
- 1 = straight gf values

## contnorm  
allows multiplicative adjustment of the continuum; useful probably only for batch syntheses the numbers employed should be around 1.0; default is 1.000000.

## plotpars  
allows you to set all of the plotting parameters if you know them in advance:
- 0 = none set (default); user can change in plotting routine
- 1 = given in following lines as follows 
    - `xlow         xhi         ylo       yhi`
    - `vshift       lamshift    obsadd    obsmult`
    - `smooth-type  FWHM-Gauss  vsini     L.D.C.    FWHM-Macro     FWHM-Loren`

## trudamp     
should moog use the detailed line damping for those transitions that have information stored in subroutine trudamp? (Default is *no*)

## veladjust   
shoud moog try to do a cross-correlation between observed and synthetic spectra and use that to align the spectra better in wavelength (Default is *no*)

## units      
controls the units in which moog outputs the final spectrum:
- 0 = angs
- 1 = microns
- 2 = 1/cm

## iraf       
allows the user to output a raw spectrum in a form suitable for IRAF's rtext input command
- 0 = don't do this, make output the normal way.
- 1 = make an IRAF-compatible output

## scat
allows the user to employ a source function which has both scattering and absorption components     
- 0 = NO scattering
- 1 = scattering

## flux/int  
choses integrated flux or central intensity
- 0 = integrated flux calculations
- 1 = central intensity calculations

## damping
here are the calculations to set up the damping; dampingopt is the damping parameter, and dampnum is the damping value read from the line list. 
For atomic lines there are several options:
- dampingopt = 0 and dampnum < 0 ---> gammav = 10^{dampnum(i)}*(T/10000K)^0.3*n_HI
- dampingopt = 0 and dampnum = 0 ---> c6 = Unsold formula
- dampingopt = 0 and dampnum > 10^(-10) ---> c6 =  (Unsold formula)*dampnum(i)
- dampingopt = 0 and dampnum(i) < 10^(-10) ---> c6 = dampnum(i)
- dampingopt = 1 ---> gammav = gamma_Barklem if possible, otherwise use dampingopt=0 options
- dampingopt = 2 ---> c6 = c6_Blackwell-group
- dampingopt = 3 and dampnum <= 10^(-10) ---> c6 = c6_NEXTGEN for H I, He I, H2
- dampingopt = 3 and dampnum > 10^(-10) ---> c6 = (c6_NEXTGEN for H I, He I, H2)*dampnum
For molecular lines (lacking a better idea) ---> c6 done as in dampingopt = 0

## obspectrum
controls the file type of the observed spectrum
- 0 = no observed spectrum is to be input
- 1 = read a true FITS file with internal read statements
- -1 = as if obspectrum = 1, but on a byte-swapping machine
- 2 = (not implemented yet)
- 3 = read a true Fits file with the FITSIO package
- 4 = (not implemented yet)
- 5 = read a special MONGO style (wavelength, flux pair) file

## histogram
makes histogram plots of observed spectra if histoyes = 1

## terminal  
gives the sm plotting window type smterm = a character string of the sm window type (see the appropriate sm manual for a list)

## plot      
decides whether or not to make a plot of results
c           0 = do not make a plot
c           For syntheses: 1 = plot only synthetic spectra
c                          2 = plot synthetic and observed spectra
c                          3 = smooth the syntheses but don't plot
c           For line analyses: # = the minimum number of lines of a 
c                                  species necessary to trigger a plot
c           For curves-of-growth: 1 = make plots
c           For flux curves: 1 = make plots

## abundances
gives the changes to be applied to the abundances
c           # = the number of different syntheses to run
c               (the next line gives the different abundance factors
c               to use)
c  minimum error check:  numatomsyn must equal numisosyn or code will stop

## isotopes   
gives the isotopes used in the line list and their abundance relative to the parent spiecies. minimum error check:  numatomsyn must equal numisosyn or code will stop

## lumratio
gives the ratio of the luminosity of two stars at a specific wavelength in a binary star system (used only with driver "binary")

## deltaradvel
gives the velocity difference between the stars binary star system (used only with driver "binary")

## synlimits 
gives the wavelength parameters for syntheses; start and sstop are beginning and ending wavelengths, step is the step size in the syntheses, and delta is the wavelength range to either side of a synthesis point to consider for line opacity calculations

## fluxlimits
gives the wavelength parameters for flux curves; start and sstop are beginning and ending wavelengths, and step is the step size in the flux curve

## blenlimits
gives the parameters for blended line abundance matches.  delwave is the wavelength offset  to the blue of first and to the red of the  last line in the blend to extend the syntheses;  step is the wavelength step size in the  computations; cogatom is the name of the element whose abundance should be varied to achieve an EW match with observations.

## coglimits 
gives the log(W/lambda) limits for curves-of-growth rwlow and rwhigh are the beginning and ending points of the log(red.width) values, rwstep is the step in log(red.width), cogatom is the declaration of which element will have its abundance varied (necessary only for spectrum synthesis curves-of-growth, and wavestep is a forced (if desired) step size in wavelength along the line (this applies to single line computations only

## limits    
old limits format...tell the user to change the keyword and quit.

## strong 
for lines which are to be considered for all of the synthesis

## opacit 
which takes the continuus opacity and scales it with the form of kaplam(i)= kaplam(i)*((factor*10000)/t(i)) in Opacit.f after it calulates the normal kaplam if value is <= 0 then it does not do it.