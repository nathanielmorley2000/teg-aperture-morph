# Overview

This data and code repository is supplemental to the manuscript:

Morley NED,* Johnston-Gramlich W, Risley Z, Baxter AG, Alexander J, Leighton LR. 2026. *Tegula funebralis* populations possess distinct aperture shapes corresponding to an energy and predation gradient. J Exp Mar Biol Ecol. 604:152223. https://doi.org/10.1016/j.jembe.2026.152223


All data and code required to replicate this analysis are contained within this repository.


# Software Requirements

"Data/TegulaOutlines.TPS" was created using `tpsUtil` v. 1.83 and `tpsDig` v. 2.32. 
"sliding-semilandmarks.R" was created using `R` v. 4.5.2. The following R packages 
are required to execute the code:
- `broom`
- `dplyr`
- `geomorph`
- `ggplot2`
- `gridExtra`
- `htmlwidgets`
- `openxlsx`
- `plotly`
- `RRPP`
- `stringr`
- `tibble`

The full references for these packages can be found in Supplementary Table 1 of the 
corresponding manuscript.


# Repository Guide

## Repository Root

The root of the repository contains the file `sliding-semilandmarks.R`, which is
used to perform the analysis and generate any plots and tables. When executed, this 
script will use data in the `Data/` subdirectory to reproduce the analyses in the `Results/` 
subdirectory.

Four additional files are included in the root of the repository to help with repository organization
and documentations: `.gitattributes`, `.gitignore`, `ReadMe.md`, and `teg-aperture-morph.Rproj`. These 
files do not contribute to the results in any way.


### Data/ subdirectory

The `Data/` subdirectory contains all of the data required to reproduce the analysis. 
All *Tegula funebralis* photos are located in the `Data/AllPhotos/` folder. `Data/TegulaOutlines.TPS` 
contains landmark data for the 30 best photographs from each locality, generated using the programs `TPSUtil` and 
`TPSDig`. `Data/Sliders.csv` contains the bounding coordinates for the sliding semi-landmarks. 
`Data/TegulaCharacteristicsList.csv` contains non-landmark data used in our analyses. 


## Results/ subdirectory

The `Results/` subdirectory contains all of the outputs from the R script. `Results/ProcrustesCoords.csv` 
contains coordinates for the Procrustes-transformed landmarks. The `Results/MorphologyResults/` folder contains
PCA and TPS outputs, as well as associated statistical tests. Interactive three-dimensional PCA plots are 
located in the `Results/MorphologyResults/*/Interactive_PCA.html` files. The 
`Results/AllometryResults/` folder contains graphical and statistical results 
relating to growth and allometry.


# Contact

If there are any issues with this code, please contact Nathaniel Morley via email 
at nmorley@ualberta.ca.

