# FunJaccR
Functional analysis of gene lists using clustering in R

## Installation

Currently FunJacc is available as a GitHub project and can be installed as follows:

1. Clone the repository

``


```
# Testing FunJacc package

library(devtools)
load_all(".")

# Regenerate docs
library(roxygen2)
roxygenise()


# Read in example gene list
gene_list <- as.vector(read.csv('test.list', header=FALSE))$V1
# This breaks things
#gene_list <- c("CDK1", "MCCA", "SOX9", "SOX13", "IL4")

# Run Funjacc
#funjacc_res <- funjacc(gene_list, data_types='all')

data_types=c('GO:CC', 'CORUM')
funjacc_res <- funjacc(gene_list, data_types='all', species='hsapiens', inflation=2)

#funjacc_res$annotation$cluster_names

create_cytoscape_network(funjacc_res, title="Test network")
```
