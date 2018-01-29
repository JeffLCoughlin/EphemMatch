# EphemMatch

The code reads in the period, epoch, positional information, etc. of all the Kepler DR25 TCEs, as well as the cumulative KOI list, and lists of EBs from the Kepler Eclipsing Binary Working Group (http://keplerebs.villanova.edu) as well as several catalogs of EBs known from ground-based surveys. The code will then perform matching to identify two different objects that have a statistically identical period and epoch (within some tolerance) and perform logic to identify which is the real source (the parent) and which is a false posivite due to contamination from the parent (a child).

For a very detailed description, see Section A.6 of Thompson et al. 2017/2018 (http://adsabs.harvard.edu/abs/2017arXiv171006758T) for the DR25 version, and read http://adsabs.harvard.edu/abs/2014AJ....147..119C which first used the technique on earlier KOI catalogs.


## Compiling and Running the Code

The EphemMatch code is provided, along with the required input files.

### Prerequisites

The code is written in C++ and only requires the standard C++ library (specifically, the required libraries are iomainip, iostream, fstream, cmath, cstdlib, ssstream, and vector). It has been tested to work with the g++ compiler, but should work with any standard C++ compile.

As part of the code, a postscript file is produced that plots all of the matches on the focal plane / CCD array, produced using Gnuplot (http://www.gnuplot.info/), which will need  be installed on your system if you want this. Else, the code should run fine, with an error that Gnuplot is not found and no plot will be produced.


### Compiling

To compile the code, use your available C++ compiler, and a recommended O2 level of optimization. For example:

```
g++ -Wno-unused-result -O2 -o match /home/jeff/KeplerJob/DR25/Code/DR25EphemMatch.cpp
```

### Running

To run EphemMatch, all filenames are hardcoded, so just type:

```
./match
```




## Citing Model-Shift

If using Model-Shift, please cite the following papers:

http://adsabs.harvard.edu/abs/2014AJ....147..119C

http://adsabs.harvard.edu/abs/2017arXiv171006758T


## Future Updates

Ideally this code would be re-written in Python and generalized so it can be used for TESS. I might one day.

