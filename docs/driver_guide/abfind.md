# abfind

force-fitting abundances to match single-line equivalent widths.

- Flow chart
![](img/pymoog/abfind_fc.PNG)
- Example

Suppose we have a list of strong lines (stored as `linelist` DataFrame) with their measured EW as follow:

```
	wavelength	id	EP	loggf	C6	D0	EW
0	11032.103	12.0	5.9460	-1.801	0.000	NaN	23.05
1	11033.662	12.0	5.9460	-2.073	0.000	NaN	13.33
2	11034.481	12.0	5.9460	-2.426	0.000	NaN	6.24
3	11005.127	14.0	6.2228	-2.098	-7.200	NaN	5.54
4	11013.703	14.0	6.2061	-0.983	-7.320	NaN	48.24
5	11017.966	14.0	6.2061	0.760	-7.320	NaN	281.02
6	11040.405	14.0	6.2061	-1.449	-7.320	NaN	22.53
7	11019.848	19.0	2.6700	-0.010	3073.329	NaN	11.13
8	11022.653	19.0	2.6703	-0.161	3073.329	NaN	7.98
9	11015.530	24.0	3.4493	-0.429	-7.530	NaN	43.10
10	11044.610	24.0	3.0112	-1.930	-7.780	NaN	5.74
11	11013.235	26.0	4.7955	-1.383	-7.550	NaN	35.22
12	11026.788	26.0	3.9434	-2.805	-7.810	NaN	11.92
13	11045.599	26.0	5.5870	-0.624	-7.500	NaN	35.73
14	11057.772	26.0	4.8349	-1.967	-7.550	NaN	11.55
15	11069.374	26.0	6.2223	-1.055	-7.140	NaN	5.06
16	11081.595	26.0	5.6425	-1.419	-7.470	NaN	7.37
17	11088.584	28.0	4.1647	-1.512	-7.560	NaN	16.36
18	11054.253	30.0	5.7958	-0.300	0.000	NaN	18.43
```

```{note}
The DataFrame has to be sorted in id and wavelength.
```

Then we can run abfind:
```py
a = pymoog.abfind.abfind(5000, 4.0, 0.0, line_list='line.list')
a.prepare_file()
pymoog.line_data.save_linelist(linelist, MOOG_run_path + 'line.list')
a.run_moog()
a.read_output()
```
The output of abfind will be stored as a dict with it keys as elements and values as DataFrame:
```
{12.0:    wavelength    ID     EP  logGF   EWin  logRWin  abund  delavg
 0   11032.103  12.0  5.946 -1.801  23.05   -5.680  7.376  -0.002
 1   11033.662  12.0  5.946 -2.073  13.33   -5.918  7.378   0.000
 2   11034.481  12.0  5.946 -2.426   6.24   -6.248  7.379   0.001,
 14.0:    wavelength    ID     EP  logGF    EWin  logRWin  abund  delavg
 0   11005.127  14.0  6.223 -2.098    5.54   -6.298  7.606  -0.030
 1   11013.703  14.0  6.206 -0.983   48.24   -5.359  7.639   0.003
 2   11017.966  14.0  6.206  0.760  281.02   -4.593  7.667   0.031
 3   11040.405  14.0  6.206 -1.449   22.53   -5.690  7.632  -0.004,
 19.0:    wavelength    ID    EP  logGF   EWin  logRWin  abund  delavg
 0   11019.848  19.0  2.67 -0.010  11.13   -5.996  4.866  -0.002
 1   11022.653  19.0  2.67 -0.161   7.98   -6.140  4.871   0.002,
 24.0:    wavelength    ID     EP  logGF   EWin  logRWin  abund  delavg
 0    11015.53  24.0  3.449 -0.429  43.10   -5.408  5.044   0.002
 1    11044.61  24.0  3.011 -1.930   5.74   -6.284  5.040  -0.002,
 26.0:    wavelength    ID     EP  logGF   EWin  logRWin  abund  delavg
 0   11013.235  26.0  4.795 -1.383  35.22   -5.495  7.129  -0.076
 1   11026.788  26.0  3.943 -2.805  11.92   -5.966  7.046  -0.159
 2   11045.599  26.0  5.587 -0.624  35.73   -5.490  7.242   0.037
 3   11057.772  26.0  4.835 -1.967  11.55   -5.981  7.171  -0.033
 4   11069.374  26.0  6.222 -1.055   5.06   -6.340  7.359   0.154
 5   11081.595  26.0  5.643 -1.419   7.37   -6.177  7.282   0.077,
 28.0:    wavelength    ID     EP  logGF   EWin  logRWin  abund  delavg
 0   11088.584  28.0  4.165 -1.512  16.36   -5.831  6.009     0.0}
```

![](img/pymoog/abfind.png)
