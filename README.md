# Repository for the growth-mortality analysis

In this repository is stored the data and scripts to reproduce the results presented in the growth-mortality manuscript.

Questions / issues? Send a mail to lionel.hertzog_@_ign.fr or open an issue on the repo.

## Reproducing the analysis

To run the analysis using RStudio, first install, if needed, the needed R packages from CRAN:

```
install.packages(c("targets", "tarchetypes", "tidyverse", "crew", "sf", "broom", "pROC", "MLmetrics",
               "ggcorrplot", "ggh4x"))
```

The analysis uses INLA, to install INLA check the instructions [on this page](https://www.r-inla.org/download/index.html), make sure that the R version match the INLA version.

Then, double-click on `Mortalite.Rproj` to open the R project in RStudio. From the console run:

```
library(targets)
tar_make()
```

This should run the whole pipeline taking approx. 45min to run.

You can then explore individual models or outputs using:

```
m <- tar_read(model_03_abg.gtot_0.ba_tot)
```

Check the _targets/objects folder for all the intermediate objects that can be read with `tar_read`.

## Structure of the repository

```
_targets.R
-- data
-- R
-- figures
-- tables
```

* _targets.R : running the analysis pipeline
* data/ : the data tables used in the analysis. NOTE: the variables `poids` and `w` have dummy values, the real values being currently not open. Please contact the author for information on access to these data.
* R/ : the R scripts used in the analysis.
* figures/ : the subfolder where the figures are saved
* tables/ : the subfolder where the tables are saved

