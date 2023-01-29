# Find dominant line

It would be useful to have a function for finding the dominant line or lines within a wavelength range.

## Find single dominant line

```py
pymoog.line_data.find_single_dominant_line(10818.3, 5000, 4.0, 0, 28000, line_list='vald_3000_24000', weedout_switch=True)
```

This function find the most dominent line in a specified wavelength. 
`NaN`s will be returned if no line is dominant. 

## Find isolated lines of an element


```py
solar_spec = pd.read_csv('solar_spec.txt', sep=' ', names=['wav', 'obs_spec'])
solar_spec['obs_spec'] /= np.median(solar_spec['obs_spec'])

line_list_use = pymoog.mpfit.find_isolated_lines(5777, 4.0, 0, [[7400, 7500]], 26, 28000, solar_spec)
```

This function is located in `mpfit` module since it requires the observed spectra.