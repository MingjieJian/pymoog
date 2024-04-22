Guide for each driver
=====================

Here we presents the useage of each drivers imported to pymoog in detail.

MOOG provides 14 drivers used for various purpose, and 9 of them are imported in `pymoog`:

+-----------------------+-----------------------+-----------------------+
| driver                | explanation           | imported              |
+=======================+=======================+=======================+
| synth                 | spectrum synthesis,   | yes                   |
|                       | varying atomic        |                       |
|                       | abundances            |                       |
+-----------------------+-----------------------+-----------------------+
| weedout               | segregation of very   | yes                   |
|                       | weak lines from       |                       |
|                       | stronger ones in a    |                       |
|                       | large line list       |                       |
+-----------------------+-----------------------+-----------------------+
| abfind                | force-fitting         | yes                   |
|                       | abundances to match   |                       |
|                       | single-line           |                       |
|                       | equivalent widths     |                       |
+-----------------------+-----------------------+-----------------------+
| blends                | force-fitting         | yes                   |
|                       | abundances to match   |                       |
|                       | blended-line          |                       |
|                       | equivalent widths     |                       |
+-----------------------+-----------------------+-----------------------+
| cog                   | curve-of-growth       | yes                   |
|                       | creation for          |                       |
|                       | individual lines      |                       |
+-----------------------+-----------------------+-----------------------+
| binary                | spectrum synthesis of | yes                   |
|                       | a binary star         |                       |
|                       | individual lines      |                       |
+-----------------------+-----------------------+-----------------------+
| doflux                | plot the overall flux | yes                   |
|                       | curve of the model    |                       |
|                       | atmosphere            |                       |
+-----------------------+-----------------------+-----------------------+
| synpop                | spectrum synthesis    | yes                   |
|                       | for a stellar         |                       |
|                       | population            |                       |
+-----------------------+-----------------------+-----------------------+
| abpop                 | equivalent width      | yes                   |
|                       | matching for a        |                       |
|                       | stellar population    |                       |
+-----------------------+-----------------------+-----------------------+
| gridsyn               | mass production of    | no                    |
|                       | synthetic spectra     |                       |
+-----------------------+-----------------------+-----------------------+
| cogsyn                | curve-of-growth       | no                    |
|                       | creation for blended  |                       |
|                       | features              |                       |
+-----------------------+-----------------------+-----------------------+
| ewfind                | calculation of        | no                    |
|                       | equivalent widths of  |                       |
|                       | individual lines      |                       |
+-----------------------+-----------------------+-----------------------+
| calmod                | converting a BEGN     | no                    |
|                       | tauross model to an   |                       |
|                       | equivalent tau5000    |                       |
|                       | scale                 |                       |
+-----------------------+-----------------------+-----------------------+
| plotit                | re-plotting of        | no                    |
|                       | spectra that were     |                       |
|                       | created in a prior    |                       |
|                       | run                   |                       |
+-----------------------+-----------------------+-----------------------+

There are two more functions which are available in `pymoog`:

- Contribution function: a function indicating where the lines are formed in the atmosphere
- MPFIT: an algorithm for self-consistent multi-parameter fitting for spectra, by `Takeda (1995) <https://ui.adsabs.harvard.edu/abs/1995PASJ...47..287T>`_ .

The following links provide the usage of each driver and function listed above. 
Only the keywords worth attention will be descripbed, and please refer to the help docstring for detail description of each keyword.

.. toctree::
   :maxdepth: 1
   :caption: Contents:

   synth <driver_guide/synth>
   weedout <driver_guide/weedout>
   abfind <driver_guide/abfind>
   blends <driver_guide/blends>
   cog <driver_guide/cog>
   find dominant line<driver_guide/find_dominant_line>
   contribution function <driver_guide/contri_func>
   mpfit <driver_guide/mpfit>
