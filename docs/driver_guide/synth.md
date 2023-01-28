# synth

Now let's see what we can control in `synth`, the core part of `pymoog`.

The following code lists all parameters which can be changed 

```py
s = pymoog.synth.synth(5060, 4.7, 0, 10540, 10620, 50000, vmicro=2, mass=1, del_wav=0.02, 
                       line_list='vald_3000_24000', weedout=False, prefix='')
s.prepare_file(model_file=None, model_format='moog', loggf_cut=None, abun_change=None, molecules=None, 
               atmosphere=1, lines=1, smooth_para=None, model_type='marcs', model_chem='st', model_geo='auto')
s.run_moog(output=False)
s.read_spectra()
```

The first 5 parameters in `pymoog.synth.synth` refers to effective temperature $T_\mathrm{eff}$, logrithm of surface gravity $\log{g}$, metallicity $\mathrm{[M/H]}$, the start and end wavelength, and resolution of the synthetic spectra.
The following parameters are:
- `vmicro`: mircotrubulance velocity. Different treatment will be implemented depending on the `model_type`; see [the model page]() for more detail. 
- `mass`: the mass of the star. Only valid when `model_type` is `marcs`.
- `del_wav`: the wavelength step of synthetic spectra.
- `line_list`: the line list to be used. It can be either the internal line lists (), a file or a line list object.
- `weedout`: if set to a value from 0 to 1, then perform the weedout task to remove weak lines (see [weedout]()); if False then do not weedout.
- `prefix`: I don't think this parameter is used now.

Parameters in `s.prepare_file`:
- `model_file`: If specified then use the model file.
- `model_format`: 
- `loggf_cut`: The cut in loggf. If specified as a float number, then the 
- `abun_change`: Apply change in abundance ratios using dict, `{}`
- `molecules`: Molecules switch in MOOG. Controls molecular equlibrium calculations. 0: do not do; 1: do calculations. 
- `atmosphere`: 
- `lines`:
- `smooth_para`:
- `model_type`: 
- `model_chem`:
- `model_geo`:
 