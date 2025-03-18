# FunJaccR
Functional analysis of gene lists using clustering in R

## Installation

FunJaccR is available as an R package and can be installed as follows:

1. Download the tarball FunJaccR_0.1.0.tar.gz

2. Install dependendies in R

```
install.packages('gRbase')
install.packages('MCL')
install.packages('gprofiler2')
install.packages("stringr")

# For drawing Cytoscape networks
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("RCy3")
```

3. Install FunJaccR on the command line
   
`R CMD install FunJaccR_0.1.0.tar.gz`

4. You can download Cytoscape from [here](https://cytoscape.org/download.html)

## Create an R script to load the package and run the test

```
# Load funjacc package
library(FunJaccR)

# Run Funjacc
funjacc_res <- funjacc(gene_list, data_types='all', species='hsapiens', inflation=2)

# Create Cytoscape network (requires an open Cytoscape session on the same machine which R is running on)
create_cytoscape_network(funjacc_res, title="Test network")
```

## Usage

The only essential input is a vector of gene names suitable for gProfiler

See `?funjacc` and `?create_cytoscape_network` for more details
