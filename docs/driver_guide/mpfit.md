# MPFIT

Under construction....

```py
solar_spec = pd.read_csv('solar_spec.txt', sep=' ', names=['wav', 'obs_spec'])
solar_spec['obs_spec'] /= np.median(solar_spec['obs_spec'])

line_list_use = find_isolated_lines(5777, 4.0, 0, [[7400, 7500]], 26, 28000, solar_spec)

line_wav = line_list_use.loc[3, 'wavelength']
del_wav = 20 / 3e5 * line_wav
wav_part = solar_spec.loc[np.abs(solar_spec['wav'] - line_wav) <= del_wav, 'wav'].values
flux_part = solar_spec.loc[np.abs(solar_spec['wav'] - line_wav) <= del_wav, 'obs_spec'].values

pymoog.mpfit.mpfit_main(wav_part, flux_part, 2, 0.4, 0, 0, {}, ['vbroad', 'm_h', 'rv'], 5777, 4.0,
                        fitting_boundary={'m_h':[-0.65, 0.65]})
```