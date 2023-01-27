# Quickstart

Here presents a minimum example for using pymoog to generate a synthetic spectra.
For more detailed usage and other functions of the code, please refer to [drivers' guide]().

Let's say that we want to generate a synthetic with $T_\mathrm{eff}=5000\,\mathrm{K}$, $\log{g}=4.0$, and metallicity $\mathrm{[M/H]}=0$ (these are the three stellar parameters you always need to provide). 
The spectra spans from $6000$ to $6200\,\mathrm{\AA}$ and with a resolution of 30000.

```py
s = pymoog.synth.synth(5000, 4.0,    0,       6000,     6200,          30000)
#                      Teff, logg, [Fe/H], wav_start(A), wav_end(A), resolution 
s.prepare_file()
s.run_moog()
s.read_spectra()
```
Then you are done! 
The synthetic spectra is stored in the object `s`:

```py
# Plot the synthesized spectra
plt.plot(s.wav, s.flux)
```

**There should be a figure here.**