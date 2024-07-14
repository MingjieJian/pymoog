# abfind: force-fitting abundances to match single-line equivalent widths.

Suppose we have a list of lines (stored as as a DataFrame called `linelist`) with their measured EW as follow (appended to the line list as the last column):

```
       wavelength    id     EP  loggf    C6  D0      EW
85992    8002.576  26.0  4.580 -1.618 -7.79 NaN   22.80
86381    8026.940  14.0  6.261 -1.004 -6.98 NaN   34.24
86400    8027.941  26.0  3.252 -2.785 -7.73 NaN   28.46
86401    8028.313  26.0  4.473 -0.686 -7.54 NaN   76.66
86520    8035.618  14.0  5.984 -1.372 -7.16 NaN   24.95
86696    8046.047  26.0  4.415 -0.100 -7.55 NaN  125.55
86735    8047.617  26.0  0.859 -4.742 -7.85 NaN   50.73
86750    8049.366  14.0  6.269 -1.287 -6.98 NaN   20.15
87069    8071.283  14.0  6.099 -1.192 -7.06 NaN   31.16
87075    8072.164  26.0  2.424 -3.519 -7.82 NaN   32.43
87105    8073.028  14.0  6.274 -1.381 -6.98 NaN   16.29
87158    8075.150  26.0  0.915 -5.088 -7.85 NaN   30.34
87245    8080.545  26.0  3.301 -2.710 -7.73 NaN   29.76
87316    8085.172  26.0  4.446 -0.121 -7.55 NaN  121.43
87399    8089.354  26.0  5.067 -1.141 -7.73 NaN   23.30
87433    8093.232  14.0  5.863 -1.075 -7.19 NaN   47.38
87518    8096.875  26.0  4.076 -1.767 -7.80 NaN   37.74
```

```{note}
The DataFrame has to be sorted in id and wavelength.
```

This line list can be loaded using:

```py
ges_linelist = pymoog.line_data.read_linelist('ges')
ges_linelist = ges_linelist[(ges_linelist['wavelength'] > 8000) & (ges_linelist['wavelength'] < 8100) & (ges_linelist['EW'] > 15)]
ges_linelist = ges_linelist[(ges_linelist['id'] == 26) | (ges_linelist['id'] == 14)]
ges_linelist
```

Then abfind can run as:
```py
a = pymoog.abfind.abfind(5777, 4.0, 0, line_list=ges_linelist)
a.prepare_file()
a.run_moog()
a.read_output()
```
The output of abfind, `a.abfind_res` will be stored as a dict with it keys as elements and values as DataFrame:
```
{14.0:    wavelength    ID     EP  logGF   EWin  logRWin  abund  delavg
 0    8071.283  14.0  6.099 -1.192  31.16   -5.413  7.497   0.016
 1    8026.940  14.0  6.261 -1.004  34.24   -5.370  7.509   0.027
 2    8035.618  14.0  5.984 -1.372  24.95   -5.508  7.448  -0.033
 3    8049.366  14.0  6.269 -1.287  20.15   -5.601  7.504   0.023
 4    8093.232  14.0  5.863 -1.075  47.38   -5.233  7.432  -0.049
 5    8073.028  14.0  6.274 -1.381  16.29   -5.695  7.497   0.016,
 26.0:     wavelength    ID     EP  logGF    EWin  logRWin  abund  delavg
 0     8089.354  26.0  5.067 -1.141   23.30   -5.541  7.387   0.053
 1     8085.172  26.0  4.446 -0.121  121.43   -4.823  7.266  -0.068
 2     8080.545  26.0  3.301 -2.710   29.76   -5.434  7.376   0.042
 3     8075.150  26.0  0.915 -5.088   30.34   -5.425  7.373   0.039
 4     8002.576  26.0  4.580 -1.618   22.80   -5.545  7.387   0.053
 5     8047.617  26.0  0.859 -4.742   50.73   -5.200  7.293  -0.041
 6     8046.047  26.0  4.415 -0.100  125.55   -4.807  7.271  -0.063
 7     8028.313  26.0  4.473 -0.686   76.66   -5.020  7.238  -0.096
 8     8027.941  26.0  3.252 -2.785   28.46   -5.450  7.380   0.046
 9     8072.164  26.0  2.424 -3.519   32.43   -5.396  7.361   0.027
 10    8096.875  26.0  4.076 -1.767   37.74   -5.332  7.340   0.006}
```

![](../img/driver_guide/abfind.png)
