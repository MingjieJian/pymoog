# synth: create synthetic spectra

Now let's see what we can control in `synth`, the core part of `pymoog`.

The example shown in quickstart only involve a few stellar parameters, i.e., $T_\mathrm{eff}$, $\log{g}$ and metallicity. 

## Changing elemental abundance ratios

The elemental abundance ratios are altered thruough `abun_change` keyword in `s.prepare_file`:

```py
s = pymoog.synth.synth(5000, 4.0,    0,       6000,     6200,          30000)
s.prepare_file(abun_change={14:0.1, 28:-0.4})
s.run_moog()
s.read_spectra()
```

`abun_change` is a dictionary with the keys as the atomic number, and value as the [X/Fe] value of the element.

## Providing your own model file

```py
s = pymoog.synth.synth(5000, 4.0,    0,       6000,     6200,          30000)
s.prepare_file(model_file='Yourmodel.mod', model_format='moog')
s.run_moog()
s.read_spectra()
```

You can provide your own model file into pymoog by specifying the `model_file` as your model file name. 
Note that the model file must be in the format of MOOG, ATLAS9, ATLAS12 or MARCS, and the `model_format` must be specified accordingly as "moog", "kurucz-atlas9", "kurucz-atlas12" or "marcs".

## Specifying the smoothing parameters 

```py
s = pymoog.synth.synth(5000, 4.0,    0,       6000,     6200,          30000)
s.prepare_file(smooth_para=['r', 0.1, 10, 0,  5, 0])
s.run_moog()
s.read_spectra()
```

The `smooth_para` kayword in `s.prepare_file` is used for specifying the smoothing parameters.
The arrangement of `smooth_para` follows the third line of `plotpars` in `batch.par` (refer to page 9 of MOOG Manual).
The first value is a one-character smoothing type for the synthetic spectra.
Possible types are: g (Gaussian), l (Lorentzian), v (rotational), m (macroturbulent), c=g+v, d=g+m, r=g+m+v
The following values are: 
- the full-width-at-half-maximum of a Gaussian smoothing function
- vsini of a rotational broadening function
- limb darkening coefficient of a rotational broadening function
- macroturbulence veloxity
- the full-width-at-half-maximum of a Lorentzian smoothing function

Note that if a type is specified, then only the value of this (or these) type(s) are valid.

`pymoog` use the full-width-at-half-maximum of a Gaussian smoothing function as the line-spread-function, and calculate the width from the input resolution.
If the second value in `smooth_para` is set as 0, then it will be over-written by the width, othewise not over-written.