# Siwdm-class
A Self Interacting Warm Dark Matter (SIWDM) extension module built on top of CLASS (https://github.com/lesgourg/class_public) version 2.7.2 

## To install the modification

Simply run `bash install.sh` on a linux terminal. This will create a folder called `class_public-2.7` which contains the modified code. Make sure that you have the `python` command bound to a version of python 2.7.x. In addition, you would want to have another version of python 3.x to run the relaxation time calculator

## Getting Started

The code itself can be run through command line, using `./class ini_file.ini pre_file.pre`. The .ini file contains all the inizialization cosmology parameters, while the .pre contains precision and other internal code variables. 

We reccomend taking a look at `explanatory_plus_siwdm.ini` for further information about possible options.

Output is directed to the `/output` directory, where you can find power spectra, Cl, etc.

## Python wrapper `classy`

The install script also installs a custom class `Class` for importing from the module `classy` onto python 2.7.x. You can find a straightforward appication of this in a supplied python script `Simple_SIWDM_PS.py`, where it is used to calculate a (low-resolution) power spectra comparison between CDM, WDM and SIWDM.

## Relaxation Time Calculator

You can find a relaxation time calculator in the file `reltime.py`, in order to create custom runs with arbitrary coupling constants (requires python 3.x). Refer to the documentation in the file itself for instructions on its usage.
