# FunJaccR
Functional analysis of gene lists using clustering in R

## Installation

Currently FunJaccR is available as a GitHub project and can be installed as follows:

1. Clone the repository

`git clone https://github.com/adamjamesreid/FunJaccR.git`

2. Install dependendies in R

```
install.packages('gRbase')
install.packages('MCL')
install.packages('gprofiler2')
install.packages("stringr")

# Optional for drawing Cytoscape networks
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("RCy3")
```

You can (optionally) get Cytoscape from [here](https://cytoscape.org/download.html)

## Create an R script to load the package and run the test

```
# Set working directory to the FunJaccR package
setwd("FunJaccR/")

# Load funjacc package
library(devtools)
load_all(".")

# Read in example gene list
gene_list <- as.vector(read.csv('test.list', header=FALSE))$V1

# Run Funjacc
funjacc_res <- funjacc(gene_list, data_types='all', species='hsapiens', inflation=2)

# Create Cytoscape network (requires RCy3 and an open Cytoscape session on the same machine as R)
create_cytoscape_network(funjacc_res, title="Test network")
```

## Usage

The only essential input is a vector of gene names suitable for gProfiler

See `?funjacc` and `?create_cytoscape_network` for more details
