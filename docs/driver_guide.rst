Guide for each driver
=============

Here we presents the useage of each drivers imported to pymoog in detail.

MOOG provides 14 drivers used for various purpose, and 5 of them are imported in `pymoog`:

|driver|explanation|imported|
|---|---|---|
|synth | spectrum synthesis, varying atomic abundances|yes|
|weedout | segregation of very weak lines from stronger ones in a large line list|yes|
|abfind | force-fitting abundances to match single-line equivalent widths|yes|
|blends | force-fitting abundances to match blended-line equivalent widths |yes|
|cog | curve-of-growth creation for individual lines|yes|
|binary | spectrum synthesis of a binary star individual lines |no|
|gridsyn | mass production of synthetic spectra|no|
|cogsyn | curve-of-growth creation for blended features |no|
|ewfind | calculation of equivalent widths of individual lines |no|
|calmod | converting a BEGN tauross model to an equivalent tau5000 scale |no|
|doflux | plot the overall flux curve of the model atmosphere |no|
|synpop | spectrum synthesis for a stellar population|no|
|abpop | equivalent width matching for a stellar population|no|
|plotit | re-plotting of spectra that were created in a prior run|no|

There are two more functions which are available in `pymoog`:
- Contribution function: a function indicating where the lines are formed in the atmosphere
- MPFIT: an algorithm for self-consistent multi-parameter fitting for spectra, by [Takeda (1995)](https://ui.adsabs.harvard.edu/abs/1995PASJ...47..287T).

.. toctree::
   :maxdepth: 1

   synth <driver_guide/grid_points>
   weedout <driver_guide/line_lists>
   abfind <driver_guide/contrbution_function>
   blends <driver_guide/keywords>
   cog <driver_guide/cog>
   contribution function <driver_guide/contri_func>
   mpfit <driver_guide/mpfit>
   Advanced usage <driver_guide/advanced_usage>