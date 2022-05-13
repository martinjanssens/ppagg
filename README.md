# ppagg
Aggregating cumulus clouds

Minimal set of files needed to run postprocessing of DALES output, for tracking and explaining the evolution of mesooscale moisture fluctuations, and concurring cloud clusters.

## Overview:

### Files needed to run the postprocessing:
- stats3d_eco.py, can be called with command-line arguments and thus wrapped in batch scripts (see runfiles)
- functions.py, ppagg_io.py contain supporting routines

### Files needed to analyse the outcomes (for BOMEX run at dx=200 and L=100 km)
- notebooks/analyse_bomex
- notebooks/analyse_bomex_2d
- plot_fielddump.py

## Fair warning
This version of the code is intended as a documentation supporting the manuscript submission Janssens et al. (2022); "Non-precipitating shallow cumulus convection is intrinsically unstable to length-scale growth" to JAS. It is not yet ready for easy-to-use postprocessing of any LES model, though we are working towards this goal. Therefore, do get in touch if you'd like to use any of this code and I (Martin) will be happy to help set it up.
