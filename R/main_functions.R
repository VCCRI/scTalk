
#######################################################################
### Generates files of cluster:ligand:receptor:cluster edge weights ###
#######################################################################
#'
#' Generates files of cluster:ligand:receptor:cluster edge weights
#'
#' For a Seurat object
#'
#' @param seurat.object a list of genes
#' @param out.label a Seurat object with cluster identities
#' @param species expression threshold for considering a gene expressed in a cell
#' @param populations.use Threshold of percentage of cells expressing the gene in a cluster
#' for it to be considered expressed
#'
#' @return NULL - results written to file
#'
GenerateEdgeWeights <- function(seurat.object,
                                out.label,
                                species,
                                populations.use = NULL) {

  if (is.null(populations.use)) {
    populations.use <- names(table(Idents(seurat.object)))
  }

  ## Read in mouse to human orthologue mappings
  extdata.path <- system.file("extdata", package = "SCTalk")
  if (species == "human") {
    mapping.file <- paste0(extdata.path, "/human_human_ensembl_gene_names.txt")
  } else if (species == "mouse") {
    mapping.file <- paste0(extdata.path, "/human_mouse_ensembl_gene_mappings.txt")
  } else{
    stop(paste("species", species, "not currently a valid option: please specify mouse or human"))
  }

  mappings = read.csv(mapping.file, header=TRUE, stringsAsFactors = FALSE)

  ## Read in the ligand-receptor pairs
  ## These were taken from Ramilowski et al. (2015) Nature Communications
  ligand.receptor.pairs = read.csv("All.Pairs-Table 1.csv", header=TRUE, row.names=1, stringsAsFactors = FALSE)

  ## Keep pairs that are either literature supported or putative but filter out those annotated as being incorrect
  ligand.receptor.pairs = ligand.receptor.pairs[ligand.receptor.pairs[ ,"Pair.Evidence"] %in% c("literature supported", "putative"), ]
  receptors = as.character(ligand.receptor.pairs[, 3])
  ligands = as.character(ligand.receptor.pairs[, 1])
  all.genes = c(receptors, ligands)

  ## Map the human gene names to mouse
  ligand.receptor.mappings = mappings[mappings[, 1] %in% all.genes, ]
  mouse.names = rownames(seurat.object)
  ligand.receptor.mappings = ligand.receptor.mappings[as.character(ligand.receptor.mappings[, 2]) %in% mouse.names, ]

  ### Remove non-unique human -> mouse mappings
  counts = table(ligand.receptor.mappings[, 1])
  counts = counts[counts==1]
  ligand.receptor.mappings = ligand.receptor.mappings[as.character(ligand.receptor.mappings[, 1]) %in% names(counts), ]
  rownames(ligand.receptor.mappings) = ligand.receptor.mappings[ , 1]

  ### Overlap pairs with the genes expressed in the data
  ligand.receptor.pairs = ligand.receptor.pairs[, c(1, 3)]
  ligand.receptor.pairs.expressed = ligand.receptor.pairs[as.character(ligand.receptor.pairs[, 1]) %in% rownames(ligand.receptor.mappings)
                                                          & as.character(ligand.receptor.pairs[, 2]) %in% rownames(ligand.receptor.mappings), ]

  ### Produce a table of weighted edges for the following:
  ### Cluster -> Ligand
  ### Ligand -> Receptor
  ### Receptor -> Cluster

  unique.genes = unique(c(as.character(ligand.receptor.pairs.expressed[,1]), as.character(ligand.receptor.pairs.expressed[, 2])))
  unique.mouse.genes = unlist(lapply(unique.genes, function(x) as.character(ligand.receptor.mappings[x, 2])))
  print("Calculating cluster-specific ligand/expression characteristics")
  de.results = calculate_cluster_specific_expression(geneList = unique.mouse.genes,
                                                     seurat.object = seurat.object,
                                                     threshold=0.1,
                                                     cluster.set = populations.use)

  ## Build a table for each ligand-receptor pair as expressed in each cluster
  print("Identifying all potential cluster:ligand:receptor:cluster paths")
  all.de.results = get_gene_pair_expression_values(ligand.receptor.pairs = ligand.receptor.pairs.expressed,
                                                   mapping.table = ligand.receptor.mappings,
                                                   cluster.names = populations.use,
                                                   gene.de.results = de.results)
  dim(all.de.results)


  ### Pull out connections for clusters of interest and add weights
  indicies = all.de.results$Cluster1 %in% populations.use & all.de.results$Cluster2 %in% populations.use

  ligand.receptor.edges = unique(all.de.results[indicies , c("Gene1", "Gene2")])
  colnames(ligand.receptor.edges) = c("source", "target")
  pair.identifiers = as.character(apply(ligand.receptor.edges, 1, function(x) paste0(x[1], ".", x[2])))
  rownames(ligand.receptor.edges) = pair.identifiers

  ## Get weights using the STRING data-base
  ligands = as.character(ligand.receptor.edges[, 1])
  receptors = as.character(ligand.receptor.edges[, 2])

  ### Here use the STRING data-base to give mouse-specific scores to ligand-receptor relationships
  lr_score_table = make_STRING_table(ligands, receptors)
  head(lr_score_table) ## print out some interactions

  overlapping.genes = intersect(pair.identifiers, rownames(lr_score_table))
  print(paste0(length(overlapping.genes), " overlaps between ligand-receptor map and STRINGdb"))

  weights = lr_score_table[overlapping.genes, ]$Combined_score
  weights = weights/1000

  ligand.receptor.edges.overlap = ligand.receptor.edges[overlapping.genes, ]
  ligand.receptor.edges.overlap$weight = weights
  ligand.receptor.edges.overlap$relationship = "ligand.receptor"
  ligands = as.character(ligand.receptor.edges[, 1])
  receptors = as.character(ligand.receptor.edges[, 2])

  ## Mapping 'source' population to corresponding ligand
  cluster.ligand.edges = unique(all.de.results[indicies , c("Cluster1", "Gene1", "Gene1.Log2_fold_change")])
  colnames(cluster.ligand.edges) = c("source", "target", "weight")
  cluster.ligand.edges$relationship = "cluster.ligand"
  cluster.ligand.edges$source = paste0("S:", cluster.ligand.edges$source)

  ## Mapping 'target' population to corresponding receptor
  receptor.cluster.edges = unique(all.de.results[indicies , c("Gene2", "Cluster2", "Gene2.Log2_fold_change")])
  colnames(receptor.cluster.edges) = c("source", "target", "weight")
  receptor.cluster.edges$relationship = "receptor.cluster"
  receptor.cluster.edges$target = paste0("T:", receptor.cluster.edges$target)

  ## Build a table of edges,
  col.order =  c("source", "target", "relationship", "weight")
  ligand.receptor.edges.overlap = ligand.receptor.edges.overlap[, col.order]
  receptor.cluster.edges = receptor.cluster.edges[, col.order]
  cluster.ligand.edges = cluster.ligand.edges[, col.order]

  all.edges = rbind(as.matrix(ligand.receptor.edges.overlap), trimws(as.matrix(receptor.cluster.edges)), trimws(as.matrix(cluster.ligand.edges)))

  ### Write the network edges to file
  write.csv(all.edges, file = paste0(out.lab, "_all_ligand_receptor_network_edges.csv"),
            row.names = FALSE, quote = FALSE)

  ## Also write out just the weights based on expression values
  expression.edges = rbind(trimws(as.matrix(receptor.cluster.edges)), trimws(as.matrix(cluster.ligand.edges)))
  write.csv(expression.edges, file = paste0(out.lab, "_ligand_receptor_weights.csv"),
            row.names = FALSE, quote = FALSE)

}

#######################################################################
### Generates files of cluster:ligand:receptor:cluster edge weights ###
#######################################################################
#'
#' Generates files of cluster:ligand:receptor:cluster edge weights
#'
#' For a Seurat object
#'
#' @param seurat.object a list of genes
#' @param out.label a Seurat object with cluster identities
#' @param species expression threshold for considering a gene expressed in a cell
#' @param populations.use Threshold of percentage of cells expressing the gene in a cluster
#' for it to be considered expressed
#'
#' @return NULL - results written to file
#'
GenerateNetworkEdges <- function(file.label = out.lab,
                                 min.weight = 1.5) {

  edge.score.file = paste0(file.label, "_all_ligand_receptor_network_edges.csv")
  score.table = read.csv(edge.score.file)

  all.weights = score.table$weight
  all.edges = score.table[, c("source", "target")]

  populations.test = clusters.use

  ## Go through and calculate summed path weights from source to target populations
  ## Minimum weight of 1.5 to select for paths with some up-regulation (in ligand, receptor or both)
  complete.path.table = c()
  for (s.pop in populations.test) {
    for (t.pop in populations.test) {
      source.population = paste0("S:", s.pop)
      target.population = paste0("T:", t.pop)
      this.path.table = get_weighted_paths(score.table, source.population = source.population,
                                           target.population = target.population, min.weight = min.weight)
      complete.path.table = rbind(complete.path.table, this.path.table)
    }
  }
  dim(complete.path.table)

  write.csv(complete.path.table, file = paste0(file.label, "_network_paths_weight", min.weight, ".csv"))

  ## Generate a background set of paths for permutation testing
  ## Set mimimum weight to -100 to capture all paths
  background.path.table = c()
  for (s.pop in populations.test) {
    for (t.pop in populations.test) {
      source.population = paste0("S:", s.pop)
      target.population = paste0("T:", t.pop)
      this.path.table = get_weighted_paths(score.table, source.population = source.population,
                                           target.population = target.population, min.weight = -100)
      background.path.table = rbind(background.path.table, this.path.table)
    }
  }
  dim(background.path.table)

  write.csv(background.path.table, file = paste0(file.label, "_background_paths.csv"))

}



