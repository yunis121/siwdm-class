# Siwdm-class
A Self Interacting Warm Dark Matter (SIWDM) extension module built on top of CLASS (https://github.com/lesgourg/class_public) version 2.7.2 

## To install the modification

Simply run ~bash install.sh~ on a linux terminal. This will create a folder called ~class_public-2.7~ which contains the modified code.

## Getting Started

The code itself can be run through command line, using ~./class ini_file.ini pre_file.pre~. The .ini file contains all the inizialization cosmology parameters, while the .pre contains precision and other internal code variables. 

We reccomend taking a look at ~explanatory_plus_siwdm.ini~ for further information about possible options.

Output is directed to the ~/output~ directory, where you can find power spectra, $C_l$s, etc.
