# cumulative-error
This code calculates the cumulative error for each umbrella.

Required files:
filename.inp: Contains the path to CONSTRAINT_XX, followed by umbrella mean and umbrella kappa in the next line.

NOTE:
CONSTRAINT_XX required to read the CV value.

USAGE:
gfortran estimate_error.f90

OUTPUT FILES:
Gives out variance.dat and delta_G.dat (contains cumulative error for each umbrella)
