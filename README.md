# Repository for the growth-mortality analysis

In this repository is stored the data and scripts to reproduce the results presented in the growth-mortality manuscript.

## Structure of the repository

`
_targets.R
-- data
-- R
-- figures
-- tables
`
* _targets.R : the analysis uses the targets package, after downloading the repo running `targets::tar_make()` from within the Rproject should run all the analysis
* data/ : the data tables used in the analysis. NOTE: the variables `poids` and `w` have dummy values, the real values being currently not open. Please contact the author for information on access to these data.
* R/ : the R scripts used in the analysis.
* figures/ : the subfolder where the figures are saved
* tables/ : the subfolder where the tables are saved

