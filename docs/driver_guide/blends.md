# blends: force-fitting abundances to match blended-line equivalent widths

Blends driver fit the EW of a feature as a whole, only varing the abunance specified. 
The following example shows a fitting of Fe abundance of the feature from 10817.7 to 10818.9A, with an observed EW of 80mA.

```py
b = pymoog.blends.blends(5000, 4.0, 0, 10817.7, 10818.9, 80, 26, line_list='vald_3000_24000')
b.prepare_file()
b.run_moog()
b.read_output()
b.blends_s_df
```

```py
wavelength	ID	EP	logGF	EWin	logRWin	abund	delavg
0	10818.275	26.0	3.96	-1.948	80.0	-5.131	7.325	0.0
```