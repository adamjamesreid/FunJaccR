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

## Understanding the output

The `funjacc()` function returns a list with two tables 'annotation' and 'network'.

The `funjacc_res$annotation` table has term ids as row names, along with the following named columns:

* term_name - Full name of the term
* term_type - Source of the term e.g. GO:BP (Gene Ontology, Biological Process), CORUM etc. see [gProfiler docs](https://biit.cs.ut.ee/gprofiler/page/docs) for more info
* p_value - gProfiler p-value for term enrichment 
* cluster_number - FunJaccR cluster number
* cluster_colour - Colour for this cluster in Cytoscape plots
* cluster_names - Name of the FunJacc cluster (term name with lowest p-value in this cluster)
* genes - Genes from your gene list which are associated with this term
* node_label_size - Used for Cytoscape plotting

The `funjacc_res$network` table represents the network of terms based on the overlap of genes associated with them, it has three columns:

* n1 - A functional term
* n2 - Another functional term
* j - Jaccard index of the lists of genes associated with the two terms
