# binary

```py
b = pymoog.binary.binary([[5777, 4, 0], [5777, 2, 0]], 10540, 10600, 50000, line_list='vald_3000_24000')
b.prepare_file(deltaradvel=10.5, lumratio=1)
b.run_moog()
b.read_spectra()
```

Using `plt.plot(b.wav, b.flux)`:
