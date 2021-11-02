# Alavi-et-al.-A-quantitative-framework-for-identifying-patterns-of-route-use-in-animal-movement-data
This repository contains r code related to Alavi et al. A quantitative framework for identifying patterns of route-use in animal movement data. Frontiers

The r code in this repository will eventually be updated, streamlined and packaged for ease of use.

 Original movement data are hosted on [Movebank](https://www.movebank.org/) (Processed data: Movebank ID 1120749252; Unprocessed data: Movebank ID 468460067)
 
 ## Instructions
There are various packages with complex dependencies used accross the scripts used in this study. In order to install every package, Windows users must first ensure that the latest version of [rtools](https://cran.r-project.org/bin/windows/Rtools/) is properly installed. Mac users must ensure that the latest version of [Xcode](https://developer.apple.com/xcode/) is properly installed.  This step cannot be skipped .

The following script installs the most recent versions of the required packages and their dependencies:

```
install.packages(c("devtools","miniCRAN","pacman"), dependencies = TRUE, type="source") 
devtools::install_github("wrathematics/getPass")
devtools::install_github("ctmm-initiative/ctmmh")


#Check if required packages and their dependencies need installation or updates
list_of_required_packages <- c("plyr", "dplyr", "lubridate", "data.table", "ggplot2", "ggpubr", "sp", 
"rgdal", "raster", "move", "spatstat", "mixR","ClusterR", "RANN", "StanHeaders", "rstan", "brms")

check_if_needs_install=as.character(miniCRAN::pkgDep(list_of_required_packages, suggests = TRUE, enhances = TRUE))
check_if_needs_update=as.character(pacman::p_update(FALSE))

new_packages <- check_if_needs_install[!(check_if_needs_install %in% installed.packages()[,"Package"])]
packages_to_update=check_if_needs_install[check_if_needs_install %in% check_if_needs_update]

packages_to_install=c(new_packages,packages_to_update)


install.packages(packages_to_install, type="source")
```
Package versions used at the time of analysis:
```getPass_0.2-3```, ```dplyr_1.0.7```, ```plyr_1.8.6```, ```data.table_1.14.2```, ```tidyr_1.1.4```, ```lubridate_1.8.0```, ```ggplot2_3.3.5```, ```ggpubr_0.4.0```, ```sp_1.4-5```, ```rgdal_1.5-27```, ```raster_3.5-2```, ```spatstat_2.2-0```, ```move_4.0.6```, ```ctmm_0.6.2```, ```mixR_0.2.0```, ```RANN_2.6.1```, ```ClusterR_1.2.5```, ```StanHeaders_2.26.1```, ```rstan_2.26.1```, ```brms_2.16.1```                               
            
Be sure to change the working directory to an appropriate one on your system. Also be sure to modify where write.csv() saves all outputs. 

## Script metadata

```route analysis frontiers.R```: This script contains majority of the code for the method presented in the paper. Path reconstruction, grid cell binning, calculation of grid level metrics, unsupervised clustering, and the calculation of the routineness score are all implimented here

```Cell DET calculations.R```: This script calculates the visits, recursion, and repeat values for each cell, which are used in ```route analysis frontiers.R```

