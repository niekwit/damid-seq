#!/bin/bash

Rscript -e 'BiocManager::install("BioVenn", update=FALSE, ask=FALSE); stopifnot(requireNamespace("BioVenn", quietly=TRUE))'
Rscript -e 'remotes::install_github("marshall-lab/damidBind", ref="v1.1.1", upgrade="never"); stopifnot(requireNamespace("damidBind", quietly=TRUE))'
