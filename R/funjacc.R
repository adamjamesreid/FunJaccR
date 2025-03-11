# AUTHOR: Adam James Reid
# Copyright (C) 2025 University of Cambridge
# This program is distributed under the terms of the GNU General Public License

library(gprofiler2)
library(stringr)
library(MCL)
library(gRbase)

# Colour for annotating clusters
colours = c('#FFC312','#C4E538','#12CBC4','#FDA7DF','#ED4C67',
            '#F79F1F','#A3CB38','#1289A7','#D980FA','#B53471',
            '#EE5A24','#009432','#0652DD','#9980FA','#833471',
            '#EA2027','#006266','#1B1464','#5758BB','#6F1E51',
            '#40407a','#706fd3','#f7f1e3','#34ace0','#33d9b2',
            '#2c2c54','#474787','#aaa69d','#227093','#218c74',
            '#ff5252','#ff793f','#d1ccc0','#ffb142','#ffda79',
            '#b33939','#cd6133','#84817a','#cc8e35','#ccae62')
default_colour = '#FDE725FF' # For if we run out of colours

#' Run FunJacc
#'
#' Determine enriched functional term clusters with FunJacc
#' @param gene_list List of genes in which to look for enriched terms
#' @param p_cut P-value cutoff for gProfiler results
#' @param jaccard_cut Jaccard index cutoff for identifying high scoring links between terms
#' @param data_types List of functional terms to include e.g.c('GO:BP', 'GO:MF', 'GO:CC', 'KEGG', 'REAC', 'TF', 'MIRNA', 'CORUM', 'HP', 'HPA', 'WP'). Supplying 'all' will include all these terms
#' @param inflation Inflation parameter for MCL clustering ~0.5-2
#' @param species Organism to use for gProfiler to interpret gene list e.g. 'hsapiens', 'mmusculus'
#' @param default_node_label_size Node size to use for Cytoscape plotting
#' @param cluster_label_size Node size for cluster labels in Cytoscape plotting
#' @return list of results elements 'annotation' = clusters and their annotation, 'network' = network of clusters
#' @examples
#' funjacc(gene_list, data_types=data_types, species='hsapiens', inflation=3)
#' @export
funjacc <- function(gene_list, p_cut = 0.01, jaccard_cut=0.5, data_types=c('GO:BP'), inflation=2, species='hsapiens',
                    default_node_label_size=5, cluster_label_size=20){

  # This allows the use of 'all' to specify all the available data types, otherwise
  # we go with whatever set is supplied
  if (length(data_types) == 1) {
    if(data_types == 'all'){
    data_types <- c('GO:BP', 'GO:MF', 'GO:CC', 'KEGG', 'REAC', 'TF', 'MIRNA', 'CORUM', 'HP', 'HPA', 'WP')
    }
  }

  print(paste0('Data types = ', data_types))
  print(paste0('P-value cutoff = ', p_cut))
  print(paste0('Jaccard index cutoff = ', jaccard_cut))
  print(paste0('Inflation parameter = ', inflation))
  print(paste0('Species = ', species))

  # Run gProfiler
  print("Running gProfiler...")
  gprof_res <- run_gprofiler(gene_list, species)

  # Create network
  print("Creating network...")
  network <- create_edges(gprof_res, data_types, p_cut, jaccard_cut)

  # Run MCL
  print("Clustering...")
  term_clusters <- find_clusters(network, inflation=inflation)

  # Annotate clusters
  print("Annotating clusters...")
  cluster_annotation <- annotate_clusters(term_clusters, gprof_res, default_node_label_size=default_node_label_size,
                                          cluster_label_size=cluster_label_size)

  return(list(annotation=cluster_annotation, network=network))
}

run_gprofiler <- function(gene_list, organism) {
  # Run gProfiler (evcodes =TRUE gives us the genes associated with each term)
  gprof_res <- gprofiler2::gost(query = gene_list,
                    organism = organism, evcodes = TRUE)

  # Filter out difficult "parents" column (it is really a list)
  gprof_res <- gprof_res$result[, -which(colnames(gprof_res$result) %in% c("parents"))]

  # Write results to a file
  #write.table(gprof_res, file=output_file, quote=FALSE, sep="\t")

  return(gprof_res)
}


create_edges <- function(gprofiler_results, data_types, p_cut, jaccard_cut){

  # Filter results based on pvalue
  # Filter results for term types of interest
  gprof_res_filt <- subset(gprofiler_results, gprofiler_results$source %in% data_types & gprofiler_results$p_value <= p_cut)

  # Create dataframe for storing overlapping terms (network edges)
  columns = c("n1","n2","j")
  network_df = data.frame(matrix(nrow = 0, ncol = length(columns)))


  # Loop over terms to find connected ones
  for (g1 in 1:length(rownames(gprof_res_filt))) {
    for (g2 in 1:length(rownames(gprof_res_filt))) {
      g1_genes <- stringr::str_split(gprof_res_filt[g1,]$intersection, ',')[[1]]
      g2_genes <- stringr::str_split(gprof_res_filt[g2,]$intersection, ',')[[1]]
      inter = intersect(g1_genes, g2_genes)
      union <- union(g1_genes, g2_genes)
      j <- length(inter) / length(union)
      if (j > 0 & g1 != g2 & j >= jaccard_cut) {
        #print(paste(g1, g2, j))
        network_df <- rbind(network_df, c(gprof_res_filt[g1,]$term_id, gprof_res_filt[g2,]$term_id, j))
      }
    }
  }

  # Add header to dataframe
  colnames(network_df) = columns

  # Add in unconnected terms
  seen_terms <- union(network_df$n1, network_df$n2)
  orphan_terms <- subset(gprof_res_filt, !(gprof_res_filt$term_id %in% seen_terms))$term_id
  for (i in 1:length(orphan_terms)) {
    network_df <- rbind(network_df, c(orphan_terms[i], orphan_terms[i], 1))
  }

  return(network_df)
}

find_clusters <- function(network_df, inflation=2, max_iter=1000){

  # Subset dataframe for first two columns, then split the dataframe by rows
  network_list <- split(network_df[,c(1,2)], f=1:nrow(network_df))
  # Get adjacency matrix (gRbase)
  adj_mat <- ug(network_list, result="matrix")

  # Run MCL
  mcl_clust <- MCL::mcl(adj_mat, addLoops = TRUE, inflation=inflation, max.iter=max_iter, allow1=TRUE)

  term_clusters <- mcl_clust$Cluster
  names(term_clusters) <- rownames(adj_mat)

  return(term_clusters)
}

annotate_clusters <- function(mcl_output, gprof_results, default_node_label_size, cluster_label_size){
  # Annotate terms with term name, cluster name, cluster colour
  # Make a dataframe
  cluster_df <- as.data.frame(mcl_output)
  # Get subset of gprofiler results we are interested in
  gprof_res_term_subset <- subset(gprof_results, gprof_results$term_id %in% rownames(cluster_df))
  # Get term names in the correct order
  cluster_df$term_name <- gprof_res_term_subset[match(rownames(cluster_df), gprof_res_term_subset$term_id),]$term_name
  # Get term types
  cluster_df$term_type <- gprof_res_term_subset[match(rownames(cluster_df), gprof_res_term_subset$term_id),]$source
  # Get p-values types
  cluster_df$p_value <- gprof_res_term_subset[match(rownames(cluster_df), gprof_res_term_subset$term_id),]$p_value


  # Renumber clusters
  old_clusters <- unique(cluster_df$mcl_output)
  new_clusters <- 1:length(old_clusters)
  names(new_clusters) <- old_clusters
  cluster_df$cluster_number <- unlist(lapply(cluster_df$mcl_output, function(x) return(new_clusters[as.character(x)] )))

  # Add cluster colours
  cluster_df$cluster_colour <- unlist(lapply(cluster_df$mcl_output, function(x) return(colours[x] )))

  # Name clusters using the term in that cluster with the lowest p value
  cluster_names <- c()
  for (i in 1:length(unique(cluster_df$cluster_number))) {
    # Subset rows to cluster i
    df_subset <- subset(cluster_df, cluster_df$cluster_number == i)
    # order based on pvalue and take top term name
    cluster_name <- df_subset[order(df_subset$p_value),]$term_name[1]
    cluster_names <- append(cluster_names, cluster_name)
  }
  # Add cluster names to data frame
  cluster_df$cluster_names <- unlist(lapply(cluster_df$cluster_number, function(x) return(cluster_names[x] )))

  # Add gene lists
  cluster_df$genes <- gprof_res_term_subset[match(rownames(cluster_df), gprof_res_term_subset$term_id),]$intersection

  # Add node label sizes for Cytoscape
  cluster_df$node_label_size = lapply(cluster_df$term_name == cluster_df$cluster_names, function(x) ifelse(x==TRUE, cluster_label_size, default_node_label_size))

  return(cluster_df)
}

#' Create a FunJacc network in Cytoscape
#'
#' Draw the FunJacc results as a network in an open instance of Cytoscape running on the same machine
#' @param funjacc_results Funjacc results
#' @param title Network title
#' @param collection Network collection
#' @return Sets up a network in your Cytoscape application
#' @examples
#' create_cytoscape_network(funjacc_res, title="Test network")
#' @export
create_cytoscape_network <- function(funjacc_results, title="my first network", collection="Funjacc Networks"){
  # load cytoscape package
  library(RCy3)
  # Check cytoscape connection
  cytoscapePing ()
  cytoscapeVersionInfo ()

  # Create nodes data
  fj_nodes <- funjacc_results$annotation
  fj_nodes$id = rownames(funjacc_results$annotation)
  fj_nodes$mcl_output <-c()

  # Create edges data
  fj_edges <- data.frame(source=funjacc_results$network$n1, target=funjacc_results$network$n2, weight=as.numeric(funjacc_results$network$j),
                         stringsAsFactors=FALSE)

  # Create network
  createNetworkFromDataFrames(fj_nodes,fj_edges, title=title, collection=collection)

  set_cytoscape_style()
}

set_cytoscape_style <- function(name="FunJaccStyle"){
  style.name = name
  defaults <- list(NODE_SHAPE="circle",
                   EDGE_TRANSPARENCY=120,
                   NODE_LABEL_POSITION="W,E,c,0.00,0.00")
  nodeLabels <- mapVisualProperty('node label','term_name','p')
  nodeLabelSize <- mapVisualProperty('node label font size','node_label_size','p')
  nodeFills <- mapVisualProperty('node fill color', 'cluster_colour', mapping.type="passthrough")
  arrowShapes <- mapVisualProperty('Edge Target Arrow Shape','interaction','d',c("activates","inhibits","interacts"),c("Arrow","T","None"))
  edgeWidth <- mapVisualProperty('edge width','weight','p')

  createVisualStyle(style.name, defaults, list(nodeLabels,nodeFills,arrowShapes,edgeWidth, nodeLabelSize))
  setVisualStyle(style.name)
}


