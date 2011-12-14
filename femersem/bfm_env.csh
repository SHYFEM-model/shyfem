#!/bin/csh
# This script sets the BFM environmental variables used in the Makefile.
# # Modify the following according to your system settings
# # The script must be executed before compilation in the current shell
# # (e.g. > source ./bfm_env.csh) or added to the .cshrc script in $HOME
#
echo "Setting BFM environment variables for host" $HOSTNAME
# BFM environmental variables
setenv BFMDIR $HOME/bfm-2.5e-system-2.2
echo "Setting the BFMDIR to "$BFMDIR

# GOTM environmental variables #
setenv GOTMDIR $HOME/gotm
echo "Setting the GOTM ROOTDIR to "$GOTMDIR
setenv GOTM_CASES $HOME/bfm-run/gotm
echo "Setting the GOTM CASES to "$GOTMDIR

# NEMO environmental variables #

# POM environmental variables #

# Compilation #
setenv FORTRAN_COMPILER IFORT
echo "Current Fortran compiler is "$FORTRAN_COMPILER
setenv COMPILATION_MODE production
echo "Current compilation mode is "$COMPILATION_MODE
setenv NETCDFINC /opt/netcdf/include
setenv NETCDFLIBDIR /opt/netcdf/lib
echo "Linking NetCDF library from:"
echo $NETCDFLIBDIR
echo $NETCDFINC
