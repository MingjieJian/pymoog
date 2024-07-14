# blends: force-fitting abundances to match blended-line equivalent widths

Blends driver fit the EW of a feature as a whole, only varing the abunance specified. 
The following example shows a fitting of Fe abundance of the feature from 8002.076 to 8003.076A, with an observed EW of 38.76mA.

```py
b = pymoog.blends.blends(5777, 4.0, 0, 8002.576-0.5, 8002.576+0.5, 38.76, 26)
b.prepare_file()
b.run_moog()
b.read_output()
```

The result is sotred in `b.blends_s_df`:

```py
   wavelength    ID    EP  logGF   EWin  logRWin  abund  delavg
0    8002.576  26.0  4.58 -1.618  38.76   -5.315  7.706     0.0
```