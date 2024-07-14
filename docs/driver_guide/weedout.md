# weedout: removing weak lines

It is sometimes useful to remove the weak lines which do not affect the spectra so much, keeping the line list short and clean.
The `weedout` driver is for this task.

After `w.read_linelist`, the strong lines are stored in `w.keep_list`, while the weak lines are sotred in `w.toss_list`.

```py
w = pymoog.weedout.weedout(5000, 2, 0, 10800, 10830, kappa_ratio=0.2)
w.prepare_file()
w.run_moog()
w.read_linelist(remove=False) # remove=False is used for w.compare(); see below.
```

The difference between the line lists can be seen by:

```py
w.compare(50000)

plt.figure(figsize=(12, 4), dpi=200)
ln1 = plt.plot(w.wav_all, w.flux_all, label='Before weedout')
ln2 = plt.plot(w.wav_keep, w.flux_keep, label='After weedout')

plt.xlabel('Wavelength ($\mathrm{\AA}$)')
plt.ylabel('Normalized flux (solid)')

plt.twinx()
ln3 = plt.plot(w.wav_keep, w.flux_keep - w.flux_all, ls='--', c='C3', label='Flux difference')
plt.ylim(top=0.1)
plt.ylabel('Flux difference (dashed)')

lns = ln1 + ln2 + ln3
labs = [l.get_label() for l in lns]
plt.legend(lns, labs, loc=4)
```

![](../img/driver_guide/weedout.png)

The lines with kappa_ratio (ratio of line to continuum opacity) smaller than the specified number is removed.
A value of 0.01 can usually remove the weak lines but keep the synthetic spectra identical.
Note that which lines are removed is subject to the stellar parameters. 