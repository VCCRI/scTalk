
#######################################################################
### Generates files of cluster:ligand:receptor:cluster edge weights ###
#######################################################################
#'
#' Generates files of cluster:ligand:receptor:cluster edge weights
#'
#' For a Seurat object, calculates
#'
#' @param seurat.object a list of genes
#' @param file.label a Seurat object with cluster identities
#' @param species expression threshold for considering a gene expressed in a cell
#' @param populations.use threshold of percentage of cells expressing the gene in a cluster for it to be considered expressed
#' @param string.dir directory for storing STRING data (defaults to working directory)
#' @param string.ver version of STRING data-base to use. Not set by default.
#' @param verbose whether to print additional information about run (default: FALSE)
#'
#' @return NULL - results written to file
#'
#' @export
#'
GenerateEdgeWeights <- function(seurat.object,
                                file.label,
                                species,
                                populations.use = NULL,
                                string.dir = NULL,
                                string.ver = NULL,
                                verbose = FALSE) {

  if (is.null(populations.use)) {
    populations.use <- names(table(Idents(seurat.object)))
  }

  ## Read in mouse to human orthologue mappings
  extdata.path <- system.file("extdata", package = "scTalk")
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
  ligand.receptor.pairs = read.csv(paste0(extdata.path, "/All.Pairs-Table 1.csv"),
                                   header=TRUE, row.names=1, stringsAsFactors = FALSE)

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
  if (verbose) {
    print(paste(nrow(all.de.results), " potential paths in data. Some examples:"))
    print(head(all.de.results))
    }


  ### Pull out connections for clusters of interest and add weights
  indicies = all.de.results$Cluster1 %in% populations.use & all.de.results$Cluster2 %in% populations.use

  ligand.receptor.edges = unique(all.de.results[indicies , c("Gene1", "Gene2")])
  colnames(ligand.receptor.edges) = c("source", "target")
  pair.identifiers = as.character(apply(ligand.receptor.edges, 1, function(x) paste0(x[1], ".", x[2])))
  rownames(ligand.receptor.edges) = pair.identifiers

  ## Get weights using the STRING data-base
  ligands = as.character(ligand.receptor.edges[, 1])
  receptors = as.character(ligand.receptor.edges[, 2])

  if (verbose) print(paste0(length(unique(ligands)), " ligands to score in STRING"))
  if (verbose) print(paste0(length(unique(receptors)), " receptors score in STRING"))

  ### Here use the STRING data-base to give mouse-specific scores to ligand-receptor relationships
  lr_score_table = make_STRING_table(ligands = ligands,
                                     receptors = receptors,
                                     dir.path = string.dir,
                                     string.ver = string.ver,
                                     verbose = verbose)
  if (verbose) print(head(lr_score_table)) ## print out some interactions

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
  write.csv(all.edges, file = paste0(file.label, "_all_ligand_receptor_network_edges.csv"),
            row.names = FALSE, quote = FALSE)

  ## Also write out just the weights based on expression values
  expression.edges = rbind(trimws(as.matrix(receptor.cluster.edges)), trimws(as.matrix(cluster.ligand.edges)))
  write.csv(expression.edges, file = paste0(file.label, "_ligand_receptor_weights.csv"),
            row.names = FALSE, quote = FALSE)

}

######################################################################
### Generates files of weigted source:ligand:receptor:target paths ###
######################################################################
#'
#' Generates files of weighted source:ligand:receptor:target paths
#'
#' Uses files generated with the GenerateEdgeWeights function to generate
#' weighted paths connecting source(clusters):ligands:receptors:target(clusters).
#' Paths are weighted by summing the three edges - source:ligand, ligand:receptor and
#' receptor:target. Paths passing a minimum weight threshold are retained. For
#' following permuation testing, a 'background' set of paths without thresholding
#' are also written out.
#'
#'
#' @param file.label a list of genes
#' @param populations.use the cell populations to use (default is all)
#' @param min.weight minimum weight for retaining a path
#' @param ncores number of cores for parallelisation
#'
#' @return NULL - results written to file
#'
#' @export
#'
#' @importFrom foreach "%dopar%"
#'
GenerateNetworkPaths <- function(file.label,
                                 populations.use = NULL,
                                 min.weight = 1.5,
                                 ncores = 1)
  {

  edge.score.file <- paste0(file.label, "_all_ligand_receptor_network_edges.csv")
  score.table <- read.csv(edge.score.file, stringsAsFactors = FALSE)

  if (is.null(populations.use)) {
    cluster.source.table <- subset(score.table, relationship == "cluster.ligand")
    source.clusters <- unique(sub(".:(.*)", "\\1", cluster.source.table$source))

    cluster.target.table <- subset(score.table, relationship == "receptor.cluster")
    target.clusters <- unique(sub(".:(.*)", "\\1", cluster.target.table$target))

    populations.use <- union(source.clusters, target.clusters)
  }

  all.weights <- score.table$weight
  all.edges <- score.table[, c("source", "target")]

  populations.test <- populations.use

  # Set up multiple workers
  system.name <- Sys.info()['sysname']
  new_cl <- FALSE
  if (system.name == "Windows") {
    new_cl <- TRUE
    cluster <- parallel::makePSOCKcluster(rep("localhost", ncores))
    doParallel::registerDoParallel(cluster)
  } else {
    doParallel::registerDoParallel(cores=ncores)
  }

  SeuratObjectFilter <- function(x) class(get(x)) == "Seurat"
  seurat.objs <- ls()[unlist(lapply(ls(), SeuratObjectFilter))]

  ## Go through and calculate summed path weights from source to target populations
  ## Minimum weight of 1.5 to select for paths with some up-regulation (in ligand, receptor or both)
  complete.path.table <- foreach::foreach(s.pop = populations.test,
                                   .combine = 'rbind') %dopar%
  {
    target.path.table <- c()
    for (t.pop in populations.test) {
      source.population <- paste0("S:", s.pop)
      target.population <- paste0("T:", t.pop)
      this.path.table <- get_weighted_paths(edge.weight.table = score.table,
                                           source.population = source.population,
                                           target.population = target.population,
                                           min.weight = min.weight)
      target.path.table <- rbind(target.path.table, this.path.table)
    }
    return(target.path.table)
  }

  write.csv(complete.path.table, file = paste0(file.label, "_network_paths_weight", min.weight, ".csv"))

  ## Generate a background set of paths for permutation testing
  ## Set mimimum weight to -100 to capture all paths
  background.path.table <- foreach::foreach(s.pop = populations.test,
                                          .combine = 'rbind',
                                          .noexport = seurat.objs) %dopar%

  {
    target.path.table <- c()
    for (t.pop in populations.test) {
      source.population <- paste0("S:", s.pop)
      target.population <- paste0("T:", t.pop)
      this.path.table <- get_weighted_paths(score.table,
                                            source.population = source.population,
                                           target.population = target.population,
                                           min.weight = -100)
      target.path.table <- rbind(target.path.table, this.path.table)
    }
    return(target.path.table)
  }
  dim(background.path.table)

  if (new_cl) { ## Shut down cluster if on Windows
    ## stop cluster
    parallel::stopCluster(cluster)
  }

  write.csv(background.path.table, file = paste0(file.label, "_background_paths.csv"))

}

####################################################
### Perform permutation testing on network edges ###
####################################################
#'
#' Perform permutation testing on network edges
#'
#' Uses permutation testing to identify cell-cell connections that have
#' path weights greater than would be expected by change.
#'
#'
#' @param file.label label for input and output files
#' @param populations.test set of cell populations to test for (default: all in file)
#' @param num.permutations number of permutations to perform
#' @param return.results whether to return the results table (default: FALSE)
#' @param ncores number of cores to use in parallelisation
#'
#' @return NULL - results written to file
#'
#' @export
#'
#' @importFrom foreach "%dopar%"
#'
EvaluateConnections <- function(file.label,
                                populations.test = NULL,
                                num.permutations = 100000,
                                return.results = FALSE,
                                ncores = 1) {

  ## Read in the individual weights for the edges in the network
  weights.file = paste0(file.label, "_all_ligand_receptor_network_edges.csv")
  weights.table = read.csv(weights.file, stringsAsFactors = FALSE)

  if (is.null(populations.test)) {
    cluster.source.table <- subset(weights.table, relationship == "cluster.ligand")
    source.clusters <- unique(sub(".:(.*)", "\\1", cluster.source.table$source))

    cluster.target.table <- subset(weights.table, relationship == "receptor.cluster")
    target.clusters <- unique(sub(".:(.*)", "\\1", cluster.target.table$target))

    populations.test <- union(source.clusters, target.clusters)
  }

  ligand.receptor.table = weights.table[weights.table$relationship == "ligand.receptor", ]
  rownames(ligand.receptor.table) = paste0(ligand.receptor.table$source, "_", ligand.receptor.table$target)

  receptor.cluster.table = weights.table[weights.table$relationship == "receptor.cluster", ]
  rownames(receptor.cluster.table) = paste0(receptor.cluster.table$source, "_", receptor.cluster.table$target)

  cluster.ligand.table = weights.table[weights.table$relationship == "cluster.ligand", ]
  rownames(cluster.ligand.table) = paste0(cluster.ligand.table$source, "_", cluster.ligand.table$target)

  ## Read in the background network file
  completeFile = paste0(file.label, "_background_paths.csv")
  background.table = read.csv(completeFile, row.names = 1)
  dim(background.table)

  ## Read in the filtered network file
  thisFile = paste0(file.label, "_network_paths_weight1.5.csv")
  complete.path.table = read.csv(thisFile, row.names = 1)

  # Set up multiple workers
  system.name <- Sys.info()['sysname']
  new_cl <- FALSE
  if (system.name == "Windows") {
    new_cl <- TRUE
    cluster <- parallel::makePSOCKcluster(rep("localhost", ncores))
    doParallel::registerDoParallel(cluster)
  } else {
    doParallel::registerDoParallel(cores=ncores)
  }

  SeuratObjectFilter <- function(x) class(get(x)) == "Seurat"
  seurat.objs <- ls()[unlist(lapply(ls(), SeuratObjectFilter))]

  ## Permutation testing for determing signifcant cell-cell connections
  ## Iterate through each combination of populations and do random selections
  ## of fold changes for ligands and receptors. Calculate empirical P-value.
  pvalue.table <- foreach::foreach(s.pop = populations.test,
                                   .combine = 'rbind',
                                   .export = c("complete.path.table", "background.table"),
                                   .noexport = seurat.objs) %dopar%
    {
      this.source = paste0("S:", s.pop)
      table.subset <- c()
      for (t.pop in populations.test) {
        this.target = paste0("T:", t.pop)

        ## Add weights together
        s.indicies = which(complete.path.table$Source == this.source)
        t.indicies = which(complete.path.table$Target == this.target)
        sub.table = complete.path.table[intersect(s.indicies, t.indicies), ]

        num.paths = nrow(sub.table)
        path.sum = sum(sub.table$Weight)

        ## Get number of total paths from background table
        s.indicies = which(background.table$Source == this.source)
        t.indicies = which(background.table$Target == this.target)
        sub.table = background.table[intersect(s.indicies, t.indicies), ]

        num_total = nrow(sub.table)

        ## Get weights from the background table
        s.indicies.background = which(background.table$Source == this.source)
        t.indicies.background = which(background.table$Target == this.target)
        combined.indicies = intersect(s.indicies.background, t.indicies.background)
        sub.background.table = background.table[combined.indicies, ]
        ligand.receptor.set = paste0(sub.background.table$Ligand, "_", sub.background.table$Receptor)

        ## Now get the real weights of ligand-receptor connections
        ligand.receptor.sub.table = ligand.receptor.table[ligand.receptor.set, ]
        ppi.weights = ligand.receptor.sub.table$weight

        random.paths = rep(NA, num_total)
        random.weight.sums = rep(NA, num_total)
        for (i in 1:num.permutations) {
          random.weights = randomise_FC_weights(ppi.weights, cluster.ligand.table, receptor.cluster.table)
          random.weights = random.weights[random.weights >= 1.5]

          this.random.weight.sum = sum(random.weights)
          this.random.path.count = length(random.weights)

          random.paths[i] = this.random.path.count
          random.weight.sums[i] = this.random.weight.sum

        }
        p.sum = sum(random.weight.sums >= path.sum)/length(random.weight.sums)

        thisLine = data.frame(Source_population = s.pop,
                              Target_population = t.pop,
                              Num_paths = num.paths,
                              Sum_path = path.sum,
                              Sum_path_pvalue = p.sum)

        table.subset <- rbind(table.subset, thisLine)

      }
      return(table.subset)
  }

  if (new_cl) { ## Shut down cluster if on Windows
    ## stop cluster
    parallel::stopCluster(cluster)
  }

  ## Do P-value adjustment for multiple testing
  pval.adj = p.adjust(pvalue.table$Sum_path_pvalue, method = "BH")
  pvalue.table$Sum_path_padj = pval.adj
  pvalue.table <- pvalue.table[order(pvalue.table$Sum_path_padj, decreasing = FALSE), ]

  ## Write test results to file
  out.file = paste0("Permutation_tests_", file.label, "_network.csv")
  write.csv(pvalue.table, file = out.file, row.names = FALSE)

  if (return.results) return(pvalue.table)

}


