# scHelper

This is a personal collection of helper functions for single-cell analyses.

# Installation

It's important to define branch, otherwise the package won't download:

devtools::install_github("rrydbirk/scHelper", ref = "dev", auth_token = <your-GitHub-PAT>)

# General information

For working with Conos:

- buildConosGraph
- createEmbeddings
- embedUMAP
- getConditionPerCell
- getConosCluster
- getConosDepth
- mitoFraction
- quickConos
- plotEmbeddingOverview

For working with Pagoda2 web applications:

- addEmbeddingP2Web
- checkDims

For working with annotations (factors):

- collapseAnnotation
- grepl.replace
- renameAnnotation

Miscellanious:

- sortCPDB # Sort CellPhoneDB results
- dotSize # Change ggplot2 geom_point size *post* definition of plot
- sget # sapply('[[', x) where x is a number
- lget # Same as sget, but lapply and x can be a character as well
- prepareObjectsForVelocity # Create objects to transfer data from R to Python
