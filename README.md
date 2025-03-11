# FunJaccR
Functional analysis of gene lists using clustering in R

## Installation

Currently FunJacc is available as a GitHub project and can be installed as follows:

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

**funjacc**

The first argument should be a vector of gene names

`species='hsapiens'`
The species used by gProfiler to interpret the gene list you provide e.g. hsapiens, mmusculus, dmelanogaster, 

`data_types=c('GO:BP', 'GO:MF', 'GO:CC', 'KEGG', 'REAC', 'TF', 'MIRNA', 'CORUM', 'HP', 'HPA', 'WP')`

The functional resources to include in the gProfiler analysis and subsequent clustering

Use 'all' to use the full list of resources

`p_cut=0.01`
The p-value cut off for gProfiler results to be included in clustering

`jaccard_cut=0.5`
The Jaccard index cut off for indentifying links between terms

`inflation=2`
The MCL inflation parameter for determining clusters

`default_node_label_size = 5`
Default size of nodes for cytoscape network

`cluster_label_size = 20`
Size of cluster label nodes for cytoscape network
