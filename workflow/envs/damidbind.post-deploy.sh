#!/bin/bash

Rscript -e 'install.packages("BioVenn", repos="https://cloud.r-project.org"); stopifnot(requireNamespace("BioVenn", quietly=TRUE))'
Rscript -e 'remotes::install_github("marshall-lab/damidBind", upgrade="never"); stopifnot(requireNamespace("damidBind", quietly=TRUE))'
