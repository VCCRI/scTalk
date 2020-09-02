

##########################################################################
### Calculate marker (differentially expressed) genes for each cluster ###
##########################################################################
#'
#' Calculate expression metrics across clusters for a set of genes
#'
#' For a provided set of genes, go through each cluster and calculate the
#' average log2-expression of each gene in the cluster, the percentage of cells
#' expressing the gene in the cluster and the log2 FC difference between the cluster
#' and the remaining clusters. Return a data-frame of the results
#'
#' @param geneList a list of genes
#' @param seurat.obect a Seurat object with cluster identities
#' @param exp.threshold expression threshold for considering a gene expressed in a cell
#' @param threshold Threshold of percentage of cells expressing the gene in a cluster
#' for it to be considered expressed
#' @param option set of clusters (default: all in the Seurat object)
#'
#' @return a data-frame of results
#'
calculate_cluster_specific_expression<-function(geneList, seurat.object, exp.threshold=0,
                                             threshold=0.5, cluster.set = NULL) {
  geneList = unique(geneList)
  results = data.frame()
  expression.matrix = Seurat::GetAssayData(seurat.object)[geneList, ]

  if (is.null(cluster.set)) {
    cluster.set = names(table(Idents(seurat.object)))
  }
  for (i in cluster.set) {

    # Get the set of cells for the target cluster
    cluster = i
    cell.set = names(Idents(seurat.object)[Idents(seurat.object) == cluster])
    remainder.set = names(Idents(seurat.object)[Idents(seurat.object) != cluster])

    # Get genes expressed in the foreground set
    gene.proportions.cluster = apply(expression.matrix[, cell.set], 1, function(x)
    {return(sum(x > exp.threshold)/length(x))})
    names(gene.proportions.cluster) = rownames(expression.matrix)
    genes.cluster = names(gene.proportions.cluster[which(gene.proportions.cluster > threshold)])
    exp.percentages = gene.proportions.cluster[genes.cluster]

    # Calculate average log2 fold-changes for the genes in the cluster
    log2.fc = getLog2FCList(seurat.object, genes.cluster, cluster)

    # Calculate averaged Log2 expression for genes in the cluster
    aveLog2Expression = getLog2AvereragedRNA(seurat.object = seurat.object, thisCluster = cluster, c.id.used = TRUE)
    aveLog2Expression = aveLog2Expression[genes.cluster]

    # build a data-frame for the results
    this.result = data.frame(Cluster = rep(i, length(genes.cluster)),
                             Gene=genes.cluster,
                             Ave_Log2_exp = aveLog2Expression,
                             Pct_expressed = exp.percentages,
                             Log2FC = log2.fc)
    rownames(this.result) = paste0(i, ".", genes.cluster)
    results = rbind(results, this.result)
  }

  return(results)
}


#########################################################################################
### For each gene pair, test if one of them is expressed above threshold in a cluster ###
#########################################################################################
#'
#' Identify expression of ligands/receptors among clusters
#'
#' For a set of ligand-receptor pairs, evaluates whether one of the gene pairs is expressed
#' for each of the clusters and builds a table mapping each ligand-expressing cluster to
#' each corresponding receptor-expressing cluster.
#'
#' @param ligand.receptor.pairs table of ligand-receptor pairs
#' @param mapping.table human:mouse gene ID mappings
#' @param cluster.names set of cluster IDs to use
#' @param gene.de.results table of expression values produced by calculate_cluster_specific_expression
#'
#' @return a data-frame of results
#'
#' @importFrom magrittr "%>%"
#'
get_gene_pair_expression_values <- function(ligand.receptor.pairs,
                                            mapping.table,
                                            cluster.names,
                                            gene.de.results) {

  results <- c()

  for (n in 1:nrow(ligand.receptor.pairs)){

    ligand = as.character(mapping.table[as.character(ligand.receptor.pairs[n, 1]), 2])
    receptor = as.character(mapping.table[as.character(ligand.receptor.pairs[n, 2]), 2])

    gene.de.results$Label = rownames(gene.de.results)
    de.res.subset <- subset(gene.de.results, Gene %in% c(ligand, receptor))
    de.res.subset %>% dplyr::mutate(Class = ifelse(Gene==ligand, "Ligand", "Receptor")) ->
      de.res.subset

    receptor.subset <- subset(de.res.subset, Class == "Receptor")
    ligand.subset <- subset(de.res.subset, Class == "Ligand")

    if (nrow(ligand.subset) > 0 & nrow(receptor.subset) > 0) {

      for (i in 1:nrow(ligand.subset)) {
        ligand.row <- ligand.subset[i, ]

        this.result <- data.frame(Cluster1 = rep(ligand.row$Cluster, nrow(receptor.subset)),
                                  Gene1 = rep(ligand.row$Gene, nrow(receptor.subset)),
                                  Gene1.Ave_Log2_exp = rep(ligand.row$Ave_Log2_exp, nrow(receptor.subset)),
                                  Gene1.Log2_fold_change = rep(ligand.row$Log2FC, nrow(receptor.subset)),
                                  Gene1.Pct_expressed = rep(ligand.row$Pct_expressed, nrow(receptor.subset)),
                                  Cluster2 = receptor.subset$Cluster,
                                  Gene2 = receptor.subset$Gene,
                                  Gene2.Ave_Log2_exp = receptor.subset$Ave_Log2_exp,
                                  Gene2.Log2_fold_change = receptor.subset$Log2FC,
                                  Gene2.Pct_expressed = receptor.subset$Pct_expressed,
                                  stringsAsFactors = FALSE
        )
        results <- rbind(results, this.result)
      }

    }
  }

  return(results)

}

#################################################################################
### Get STRINGdb scores: given a list of ligands and list of receptors, pulls ###
### out associations from STRINGdb and builds a ligand-receptor scoring table ###
#################################################################################

#' Ligand:receptor association scores
#'
#' Generates a table of ligand:receptor association scores from the STRING database
#'
#' @param ligands a vector of ligands
#' @param receptors a vector of receptors
#' @param species species to use - either mouse (default) or human
#' @param dir.path output directory for storing STRINGdb data. If not provided uses current working directory
#' @param string.ver STRING version to use. Default is unspecified (NULL).
#' @param verbose whether to print additional information about run (default: FALSE)
#' @param string.receptors a list of receptor names that STRING recognises (if cannot map defaults)
#' @param string.ligands a list of ligand names that STRING recegnises (if cannot map defaults)
#' @param string.input.map named vector of STRING-compatable gene symbols with names corresponding to genes in dataset
#' @return a data-frame containing ligands, receptors and STRING association scores between them.
#'
#' @examples
#'
#' make_STRING_table(ligands, receptors, species="mouse")
#'
#' ## Or if STRING doesn't recognize some genes
#' string.input.map <- c("Ackr3)
#' names(string.input.map) <- c("Cxcr7")
#' string.receptors <- c("Cxcr7)
#' make_STRING_table(ligands, receptors, species="mouse", string.receptors=string.receptors, string.input.map=string.input.map)
#'
make_STRING_table <- function(ligands,
                              receptors,
                              species,
                              dir.path = NULL,
                              string.ver = NULL,
                              verbose = FALSE,
                              string.receptors = NULL,
                              string.ligands = NULL,
                              string.input.map = NULL) {

  ## If user does not provide a directory use the working directory
  if (is.null(dir.path)) {
    dir.path = getwd()
  }

  ## Check species and set species code and STRING data directory accordingly
  if (species == "mouse") {
    species.id = 10090
    string.dir = file.path(dir.path, "STRINGdb_mouse")
  } else if (species == "human") {
    species.id = 9606
    string.dir = file.path(dir.path, "STRINGdb_human")
  } else{
    stop("invalid species input - either use 'mouse' or 'human'")
  }

  ## Check if the directory already exists and create if not
  if (!dir.exists(string.dir)) {
    dir.create(string.dir)
    }

  if (verbose) print("Connecting to STRING...")
  ## Connect to the STRINGdb
  if (is.null(string.ver)) {
    string_db <- STRINGdb::STRINGdb$new(species=species.id, score_threshold=0, input_directory=string.dir)

  } else {
    string_db <- STRINGdb::STRINGdb$new(version=string.ver, species=species.id, score_threshold=0, input_directory=string.dir)
  }

  ## For this analysis only one gene that needs to be mapped
  gene.table = data.frame(Gene = as.character(unique(ligands)), Class = "ligand")
  gene.table = rbind(gene.table, (data.frame(Gene = as.character(unique(receptors)), Class = "receptor")))
  if (!is.null(string.ligands)) {gene.table = rbind(gene.table, data.frame(Gene=string.ligands, Class="ligand"))}
  if (!is.null(string.receptors)) {gene.table = rbind(gene.table, data.frame(Gene=string.ligands, Class="receptor"))}

  if (verbose) print("Getting STRING identifier map")
  gene.table$GeneOrig <- gene.table$Gene
  gene.table.mapped <- string_db$map(gene.table, "Gene", removeUnmappedRows = TRUE )
  if (verbose) print("Identifier map retrieved")
  gene.table.mapped %>% dplyr::distinct(STRING_id, .keep_all = TRUE) -> gene.table.mapped
  rownames(gene.table.mapped) = gene.table.mapped$STRING_id

  ## Switch original gene names back
  gene.table.mapped$Gene <- gene.table.mapped$GeneOrig

  if (verbose) print(paste0(nrow(gene.table.mapped), " STRING results retrieved. Some examples:"))
  if (verbose) print(head(gene.table.mapped))

  ## Get the list of STRING identifiers
  hits <- gene.table.mapped$STRING_id

  ## Interaction table for identifier list
  gene.interactions.table = string_db$get_interactions(hits)
  print(paste0(nrow(gene.interactions.table), " interactions identified from STRINGdb"))

  ## gene interactions table are using STRING ID so need to convert back to gene name
  source.names = unlist(lapply(gene.interactions.table$from, function(x) gene.table.mapped[x, ]$Gene))
  target.names = unlist(lapply(gene.interactions.table$to, function(x) gene.table.mapped[x, ]$Gene))

  gene.interactions.table$from.gene = source.names
  gene.interactions.table$to.gene = target.names

  gene.interactions.table.subset = gene.interactions.table[, c("from.gene", "to.gene", "combined_score")]

  if (!is.null(string.input.map)) {
    ## Replace STRING gene names with original gene names in the from column
    replace.index = gene.interactions.table.subset$from.gene %in% names(string.chromium.map)
    if (sum(replace.index) > 0) {
      gene.interactions.table.subset[replace.index, 1] = as.character(string.chromium.map[gene.interactions.table.subset[replace.index, 1]])
    }

    ## Replace STRING gene names with original gene names in the to column
    replace.index = gene.interactions.table.subset$to.gene %in% names(string.chromium.map)
    if (sum(replace.index) > 0) {
      gene.interactions.table.subset[replace.index, 2] = as.character(string.chromium.map[gene.interactions.table.subset[replace.index, 1]])
    }
  }

  all.evidence.labels <- c("neighborhood", "neighborhood_transferred",
                           "fusion", "cooccurence", "homology", "coexpression",
                           "coexpression_transferred", "experiments",
                           "experiments_transferred", "database",
                           "database_transferred", "textmining", "textmining_transferred",
                           "combined_score")

  all.evidence.labels <- intersect(all.evidence.labels, colnames(gene.interactions.table.subset))

  ## Now identify all ligand-receptor interactions from the input table
  lr_score_table = data.frame()

  ## re-order the to-from relationships to indicate ligand-receptor direction
  for (i in c(1:nrow(gene.interactions.table.subset))) {
    gene1 = gene.interactions.table.subset[i, 1]
    gene2 = gene.interactions.table.subset[i, 2]
    scores = gene.interactions.table.subset[i, all.evidence.labels]
    if ((gene1 %in% ligands) & (gene2 %in% receptors)) {
      ## Keep order
      this.row = data.frame(Ligand = gene1, Receptor = gene2)
      this.row <- cbind(this.row, scores)
      rownames(this.row) = paste0(gene1, ".", gene2)
      lr_score_table = rbind(lr_score_table, this.row)
    } else if ((gene2 %in% ligands) & (gene1 %in% receptors)) {
      ## Reverse the order
      this.row = data.frame(Ligand = gene2, Receptor = gene1)
      this.row <- cbind(this.row, scores)
      rownames(this.row) = paste0(gene2, ".", gene1)
      lr_score_table = rbind(lr_score_table, this.row)
    } ## Ignore other combinations
  }
  nrow(lr_score_table)
  head(lr_score_table)

  return(lr_score_table)
}


##############################################################################
### Given a set of weighted edges, build a network and find weighted paths ###
### between query source nodes and all other nodes                         ###
##############################################################################
#' Determine weights paths between query source and target populations
#'
#' Given a table of weights for source:ligands, ligands:receptors and receptors:targets and a query
#' source and target population, build a weighted, directed graph (using the igraph package) and pull
#' out summed paths between the query source and target populations. Retain paths above a user-provided
#' minimum weight.
#'
#' @param edge.weight.table a table of edges and weights for building the source:ligand:receptor:target graph
#' @param source.population source (ligand-expressing) population
#' @param target.population target (receptor-expressing) population
#' @param min.weight minimum weight for retaining a source:target path (deafult is 1.5)
#' @param print.num.connections whether to print the number of paths passing the min weight threshold (False by default)
#' @return a vector of source:target path weights
#' @examples
#' getWeightedPaths(edge.weight.table, source.population, target.population)
get_weighted_paths <- function(edge.weight.table, source.population, target.population,
                             min.weight = 1.5, print.num.connections = FALSE) {
  ## put edges in format for creating a graph
  all.weights = edge.weight.table$weight
  all.edges = edge.weight.table[, c("source", "target")]

  edges.list = unlist(lapply(t(all.edges), function(x) c(x[1])))

  ## Make directed graph and add edge weights
  lr.plot <- igraph::make_graph(edges.list, directed = TRUE)
  igraph::E(lr.plot)$weight <- all.weights

  ## Calculate all shortest paths (ignoring weight) between source and target node
  all.paths = igraph::all_shortest_paths(lr.plot,
                                         from = source.population,
                                         to = target.population,
                                         mode = "out",
                                         weights = NA)

  ## Calculate sum of weights for all shortest paths
  path.weight.table = matrix(0, ncol = 5, nrow = length(all.paths$res))
  colnames(path.weight.table) = c("Source", "Ligand", "Receptor", "Target", "Weight")
  path.weight.table = as.data.frame(path.weight.table)
  for (i in c(1:length(all.paths$res))) {
    node.names = names(unlist(all.paths$res[i]))
    this.epath = igraph::E(lr.plot, path=unlist(all.paths$res[i]))
    this.weight = sum(igraph::E(lr.plot)$weight[this.epath])
    this.elem = c(node.names[1], node.names[2], node.names[3], node.names[4], this.weight)
    path.weight.table[i, ] = this.elem
  }
  path.weight.table$Weight = as.numeric(path.weight.table$Weight)

  ## Filter paths by a minimum weight
  path.weight.table = path.weight.table[path.weight.table$Weight >= min.weight, ]
  if (print.num.connections) {
    print(paste0(nrow(path.weight.table), " paths for ", source.population, " to ", target.population))
  }

  return(path.weight.table)
}


####################################################################################
### Given a set of real ligand-receptor connections and weights, randomly        ###
### select fold-changes from the cluster-ligand and receptor-cluster connections ###
####################################################################################
#' Randomise source:ligand and receptor:target weights
#'
#' Generates a randomisation of the network by shuffling the source:ligand and receptor:target weights,
#' which are dervied from fold-change expression. The ligand:receptor weights derived from STRINGdb
#' remain unshuffled. The weights are then summed to represent randomised source:target path weights.
#'
#' @param ppi.weights a vector of ligands
#' @param cluster.ligand.table a vector of receptors
#' @param receptor.cluster.table species to use - either mouse (default) or human
#' @return a vector of source:target path weights
#' @examples
#' RandomiseFCWeights(ppi.weights, cluster.ligand.table, receptor.cluster.table)
randomise_FC_weights <- function(ppi.weights, cluster.ligand.table, receptor.cluster.table) {

  sample.size = length(ppi.weights)

  ## Get the set of randomly selected ligands
  rand.weights1 = sample(cluster.ligand.table$weight, size = sample.size, replace = FALSE)

  ## Get the set of randomly selected receptors
  rand.weights2 = sample(receptor.cluster.table$weight, size = sample.size, replace = FALSE)

  ## Add the weights together and return
  combined.weights = rand.weights1 + ppi.weights + rand.weights2

  return(combined.weights)
}

##################################
### Calculate Log2 fold-change ###
##################################
getLog2FC = function(x, y) {
  ## Assumes data has previous been processed as log(x+1)
  return(log2(mean(exp(x))/mean(exp(y))))
}

#########################################################################################
### Get a list of Log2 fold-change values for a gene list between a specified cluster ###
### and either all remaining cells (default) or an alternative cluster                ###
#########################################################################################
getLog2FCList <- function(seurat.object, geneList, cluster1, cluster2=NULL) {

  ### Need to check whether the input is a cluster label or vector of cell names
  if (length(cluster1) == 1){ # cluster identity used as input
    foreground.set = names(Idents(seurat.object)[Idents(seurat.object)==cluster1])
  } else { # cell identity used as input
    foreground.set = cluster1
  }
  if (is.null(cluster2)) {
    remainder.set = setdiff(colnames(seurat.object), foreground.set)
  } else {
    if (length(cluster2) == 1) { # cluster identity used as input
      remainder.set = names(Idents(seurat.object)[Idents(seurat.object)==cluster2])
    } else { # cell identity used as input
      remainder.set = cluster2
    }
  }
  log2.fc.values = apply(Seurat::GetAssayData(seurat.object)[geneList, c(foreground.set, remainder.set)], 1, function(x)
    getLog2FC(x[foreground.set], x[remainder.set]))
  return(log2.fc.values)
}

Log2ExpMean <- function (x) {
  return(log2(x = mean(x = exp(x = x) - 1) + 1))
}

####################################################################
### return an average of the log2-transformed expression data ###
####################################################################
getLog2AvereragedRNA <- function(seurat.object, thisCluster, c.id.used = FALSE){

  if (length(thisCluster) == 0) {
    detection.rate = rep(0, nrow(Seurat::GetAssayData(seurat.object)))
    names(detection.rate) = rownames(Seurat::GetAssayData(seurat.object))
    return(detection.rate)
  }

  if (c.id.used == TRUE){ # cluster identity used as input
    cell.set = names(Idents(seurat.object)[Idents(seurat.object)==thisCluster])
  } else { # cell identity used as input
    cell.set = thisCluster
  }

  if (length(cell.set) == 1) {
    return(Seurat::GetAssayData(seurat.object)[, cell.set])
  }

  average.rna = apply(Seurat::GetAssayData(seurat.object)[, cell.set], 1, function(x) Log2ExpMean(x))
  return(average.rna)
}



