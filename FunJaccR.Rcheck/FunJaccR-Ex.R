pkgname <- "FunJaccR"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('FunJaccR')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("create_cytoscape_network")
### * create_cytoscape_network

flush(stderr()); flush(stdout())

### Name: create_cytoscape_network
### Title: Create a FunJacc network in Cytoscape
### Aliases: create_cytoscape_network

### ** Examples

create_cytoscape_network(funjacc_res, title="Test network")



cleanEx()
nameEx("funjacc")
### * funjacc

flush(stderr()); flush(stdout())

### Name: funjacc
### Title: Run FunJacc
### Aliases: funjacc

### ** Examples

funjacc(gene_list, data_types='all', species='hsapiens', inflation=2)



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
