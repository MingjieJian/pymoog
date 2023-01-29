# weedout: removing weak lines

It is sometimes useful to remove the weak lines which do not affect the spectra so much, keeping the line list short and clean.
The `weedout` driver is for this task.

```py
w = pymoog.weedout.weedout(5000, 4.0, 0, 10830-15, 10830+15, kappa_ratio=0.2, line_list='vald_3000_11000',)
w.prepare_file()
w.run_moog()
w.compare(50000)
plt.plot(w.wav_all, w.flux_all, label='before weedout')
plt.plot(w.wav_keep, w.flux_keep, label='after weedout')
```

![](../img/driver_guide/weedout.png)

The lines with kappa_ratio (ratio of line to continuum opacity) smaller than the specified number is removed.
A value of 0.01 can usually remove the weak lines but keep the synthetic spectra identical.
Note that which lines are removed is subject to the stellar parameters. 