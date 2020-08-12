# GORT

A C implementaton of the Geometric Optic Radiative Transfer model as described in Ni et al. (1999). 
It was originally used as the forward operator fr data assimilatin experiments described by 
Quaife et al. (2008).
This version coupled to PROSPECT-D to (optionally) provide leaf optical properties.

## Authors

Dr Tristan Quaife & Dr Wenge Ni-Meister

With the exception of the contents of the PROSPECT-D directory all code was written by Tristan based on the equations described in Ni et al. (1999).


## Usage
 
usage: `gortt \[options\] < angles.dat > output.dat`

A list of \[options\] and brief description of the formet of the angles.dat file can be found by typing:

gortt -u


## Example

The following example calculates the reflectance for a forest with an LAI of 4.0 at wavelengths of 450nm, 600nm and 1000nm for a view zenith of 10 degress, solar zenith of 30 degrees and relative azimuth of 20 degrees.

`echo -en "1 4 450 600 800 1000\n10 0 30 20\n"|./gortt -LAI 4.0`


## Install

I have only ever tried installing on Linux platforms. As long as you have the standard development tools installed (specifically `gcc`, `gfortran` and `make`) it should just be a case of typing `make`.


## References

Ni, W., Li, X., Woodcock, C. E., Caetano, M. R., & Strahler, A. H. (1999). An analytical hybrid GORT model for bidirectional reflectance over discontinuous plant canopies. IEEE Transactions on Geoscience and Remote Sensing, 37(2), 987-999.

Quaife, T., Lewis, P., De Kauwe, M., Williams, M., Law, B. E., Disney, M., & Bowyer, P. (2008). Assimilating canopy reflectance data into an ecosystem model with an Ensemble Kalman Filter. Remote Sensing of Environment, 112(4), 1347-1364.
