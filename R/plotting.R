
#################################################################
### Generate a figure showing inbound vs outbound connections ###
#################################################################

#' Plot inbound vs outbound weights
#'
#' Generates a dot plot showing inbound vs outbound weights for all populations.
#'
#' @param input.file a table of source:ligand:receptor:target paths.
#' @param col.set optional paramater to set colours of the populations.
#' @param p.labels cell-type labels (defaults to labels in input.file)
#' @param p.adj.thresh minimum adjusted P-value for including a cell-cell connection
#' @param lab.weight.thresh minimum weight for labelling populations in the plot.
#' @param return.plot whether to return the ggplot2 object (default: TRUE)
#' @return a ggplot2 object.
#' @examples
#' InboundOutboundPlot(path.table = example.table)
#'
#' @import ggplot2
#'
#' @export
#'
InboundOutboundPlot <- function(input.file,
                                col.set = NULL,
                                p.labels = NULL,
                                p.adj.thresh = 0.05,
                                lab.weight.thresh = 1000,
                                return.plot = TRUE) {

  path.sig.file = input.file

  sig.paths = read.csv(path.sig.file, sep = ",", stringsAsFactors = FALSE)

  ## Subset according to adjusted P-value
  sig.paths <- subset(sig.paths, Sum_path_padj < p.adj.thresh)

  if (is.null(p.labels)) {
    p.labels <- union(sig.paths$Source_population, sig.paths$Target_population)
    p.labels <- p.labels[order(p.labels)]
  }

  if (is.null(col.set)) {
    col.set = scales::hue_pal()(length(p.labels))
  }

  col.set <- col.set[p.labels]

  edge.sources = factor(sig.paths$Source_population, levels = p.labels)
  edge.targets = factor(sig.paths$Target_population, levels = p.labels)

  source.counts = table(edge.sources)
  target.counts = table(edge.targets)

  outgoing.weights = c()
  incoming.weights = c()
  for (this.pop in p.labels) {
    outgoing = sum(sig.paths[which(sig.paths$Source_population == this.pop), ]$Sum_path)
    outgoing.weights = append(outgoing.weights, outgoing)

    incoming = sum(sig.paths[which(sig.paths$Target_population == this.pop), ]$Sum_path)
    incoming.weights = append(incoming.weights, incoming)
  }

  ggData = data.frame(Population = factor(p.labels, levels = p.labels), Outbound = as.numeric(source.counts),
                      WeightsOut = outgoing.weights, Inbound = as.numeric(target.counts), WeightsIn = incoming.weights)
  ggData$Population = factor(ggData$Population, levels = p.labels)
  pl <- ggplot(ggData, aes(x = WeightsIn, y = WeightsOut, colour = Population)) + geom_point(size = 5) +
    scale_colour_manual(values = col.set) + theme_bw(base_size = 15) + ylab("Total outgoing weights") +
    xlab("Total incoming weights") +
    geom_text(data=subset(ggData, (WeightsOut > lab.weight.thresh | WeightsIn > lab.weight.thresh)),
              aes(x=WeightsIn, y = WeightsOut, label=Population), hjust=0.4, vjust=-0.6, colour = "black") +
    guides(color = guide_legend(override.aes = list(size=4))) +
    theme(legend.position = "right")

  if (return.plot) return(pl)
}

#########################################################################################
### For a population generate a bar plot showing ligands contributing greatest weight ###
#########################################################################################

#' Plot a bar graph of the top ligands
#'
#' @param input.file a table of source:ligand:receptor:target paths.
#' @param cell.identity the name of a cell population/cluster to use
#' @param col.use optional paramater to set colour of the bar graph.
#' @param min.weight minimum weight for including paths in the top ligand calculations.
#' @param max.plot.num the maximum number of ligands to plot
#' @return a ggplot2 object.
#' @examples
#' PlotTopLigands(path.table = example.table, cell.identity = '1')
#'
#' @import ggplot2
#'
#' @export
#'
PlotTopLigands <- function(input.file,
                           cell.identity,
                           col.use = "#757575",
                           min.weight = 2,
                           max.plot.num = 15,
                           title.use = NULL) {

  path.table = read.csv(input.file, sep = ",", row.names = 1, stringsAsFactors = FALSE)

  source.population = paste0("S:", cell.identity)
  paths.subset = subset(path.table, Source == source.population)

  # SET MIN WEIGHT
  paths.subset.high = subset(paths.subset, Weight > min.weight)
  ligand.set = unique(paths.subset.high$Ligand)
  ligand.table = table(paths.subset.high$Ligand)

  ligand.weights.table = c()
  for (gene in ligand.set) {
    table.subset = subset(paths.subset.high, Ligand == gene)
    summed.weight = sum(table.subset$Weight)
    this.row = data.frame(Ligand = gene, Weight = summed.weight)
    ligand.weights.table = rbind(ligand.weights.table, this.row)
  }

  ligand.weights.table = ligand.weights.table[order(ligand.weights.table$Weight, decreasing = TRUE), ]

  # Check number of ligands - subset table if greater than max plot number
  if (nrow(ligand.weights.table) > max.plot.num) {
    ligand.weights.table = ligand.weights.table[1:max.plot.num, ]
  }

  pl <- ggplot(ligand.weights.table, aes(x=reorder(Ligand, -Weight), y=Weight)) +
    geom_bar(stat="identity", fill = col.use) + theme_classic(base_size = 16) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    xlab("Ligand") + ylab("Summed weight")
  if (is.null(title.use)) {
    pl <- pl + ggtitle(paste0("Summed weights for top outbound ligands"))
  } else{
    pl <- pl + ggtitle(title.use)
  }

  return(pl)
}


###########################################################################################
### Rescales the X-coordinates of nodes in a tree graph to space them more aestheically ###
### Input: 'graph' is required to be a ggraph object                                    ###
###########################################################################################
#' Rescale nodes in a ggraph tree for even distribution
#'
#' Takes as input a ggraph object containing a tree graph with 1 or more source nodes at the root
#' connecting to targets in the leaf nodes with connecting ligands and receptors as intermediates.
#' Rescales the X-coordinates of the nodes to ensure an even spacing. Rescales according to a reference
#' row - ligands by default, but probably works best when specified as the row with largets number
#' of nodes. Can select rows to rescale, though all are done by default.
#'
#' @param graph a ggraph object containing a tree graph.
#' @param scale.by a reference row to scale to. Ligand by default.
#' @param do.update select graph rows to updates. All by default.
#' @return a ggraph object with node X-coordinates rescaled.
#' @examples
#' rescale_target_node_coordinates(tree.graph)
rescale_node_coordinates <- function(graph, scale.by = "Ligand", do.update = c("Source", "Ligand", "Receptor", "Target")) {
  graph.data = graph$data
  scale.indicies = which(graph.data$Level == scale.by)
  these.coords = graph.data[which(graph.data$Level == scale.by), ]$x

  ## First get the set of ligand X-coordinates
  if ("Ligand" %in% do.update) {
    ligand.indicies = which(graph.data$Level == "Ligand")
    ligand.coords = graph.data[which(graph.data$Level == "Ligand"), ]$x

    if (length(ligand.coords) == 1) {
      new.ligand.coordinates = median(these.coords)
    } else{
      new.ligand.coordinates <- seq(from = min(these.coords), to = max(these.coords),
                                    length.out = length(ligand.coords))
    }
    graph$data[ligand.indicies, ]$x = new.ligand.coordinates
  }

  ## Second get the set of receptor X-coordinates
  receptor.indices = which(graph.data$Level == "Receptor")
  receptor.coords = graph.data[receptor.indices, ]$x

  ## Re-scale the Receptor coordinates to be in-line with ligands
  #new.receptor.coordinates <- rescale(receptor.coords, to = c(min(ligand.coords), max(ligand.coords)))
  new.receptor.coordinates <- seq(from = min(ligand.coords), to = max(ligand.coords),
                                  length.out = length(receptor.coords))

  ## Update the x-coordinates for the receptor nodes
  graph$data[receptor.indices, ]$x = new.receptor.coordinates

  ## Now pull out the current target node X-coordinates
  target.indices = which(graph.data$Level == "Target")
  target.coords = graph.data[target.indices, ]$x

  ## Re-scale the Target coordinates to be in-line with receptors
  #new.target.coordinates <- rescale(target.coords, to = c(min(new.receptor.coordinates), max(new.receptor.coordinates)))
  new.target.coordinates <- seq(from = min(new.receptor.coordinates), to = max(new.receptor.coordinates),
                                length.out = length(target.coords))

  ## Update the x-coordinates for the target nodes
  graph$data[target.indices, ]$x = new.target.coordinates

  ## Finally set the X-coordinate for the 'root' node/s
  if ("Source" %in% do.update) {
    scale.indicies = which(graph$data$Level == scale.by)
    these.coords = graph$data[which(graph$data$Level == scale.by), ]$x

    source.indicies = which(graph.data$Level == "Source")
    source.coords = graph.data[which(graph.data$Level == "Source"), ]$x
    if (length(source.coords) == 1) {
      new.source.coordinates = median(these.coords)
    } else{
      new.source.coordinates <- seq(from = min(these.coords), to = max(these.coords),
                                    length.out = length(source.coords))
    }
    graph$data[which(graph.data$Level == "Source"), ]$x = new.source.coordinates
  }

  return(graph)
}


###################################################################
### Generate a circle plot of significant cell-cell connections ###
###################################################################
#'
#' Generate a circle plot of significant cell-cell connections
#'
#' Given a file generated by the EvaluateConnections function,
#' subsets according to adjusted P-value and plots the significant
#' cell-cell connections. Uses the igraph package for network plotting.
#'
#' @param input.file input file
#' @param adj_pval_thresh adjusted p-value threshold for defining significance
#' @param col.set named vector of colour codes to fill node colours
#' @param lab.weight.thresh minimum weight for labelling populations in the plot.
#' @return a ggplot2 object.
#' @examples
#' PlotTopLigands(path.table = example.table, this.population = '1')
#'
#' @export
#'
#' @importFrom igraph layout_in_circle
#'
CellCirclePlot <- function(input.file,
                           adj_pval_thresh = 0.05,
                           col.set = NULL,
                           arrow.size = 0.6,
                           arrow.width = 2.0,
                           edge.multi = 0.02)
  {

  path.sig.file = input.file

  sig.paths = read.csv(path.sig.file, sep = ",", stringsAsFactors = FALSE)
  #sig.paths$Source_population = sub(".:(.*)", "\\1", sig.paths$Source_population)
  #sig.paths$Target_population = sub(".:(.*)", "\\1", sig.paths$Target_population)

  ## Define population colours for plotting if not provided
  if (is.null(col.set)) {
    all.nodes <- union(sig.paths$Source_population, sig.paths$Target_population)
    col.set <- scales::hue_pal()(length(all.nodes))
    names(col.set) <- all.nodes
  }

  ### Set P-value threshold to use
  sig.paths <- subset(sig.paths, Sum_path_padj < adj_pval_thresh)
  rownames(sig.paths) <- paste0(sig.paths$Source_population, "_", sig.paths$Target_population)

  ## sub-set the collapsed weight table by the edges found to be significant
  collapsed.path.table <- sig.paths[, c("Source_population", "Target_population", "Sum_path")]
  colnames(collapsed.path.table) <- c("Source", "Target", "Weight")

  ## put edges in format for creating a graph
  all.weights <- collapsed.path.table$Weight
  all.edges <- collapsed.path.table[, c("Source", "Target")]
  edge.sources <- all.edges$Source
  edges.list <- unlist(lapply(t(all.edges), function(x) c(x[1])))

  ## Make directed graph and add edge weights
  lr.plot <- igraph::make_graph(edges.list, directed = TRUE)
  igraph::E(lr.plot)$weight <- all.weights

  #population.cluster.map = names(new.ident)
  #names(population.cluster.map) = as.character(new.ident)

  populations.use <- igraph::as_ids(igraph::V(lr.plot))
  col.set <- col.set[populations.use]
  cluster.colors <- data.frame(cluster=populations.use, color=col.set)
  rownames(cluster.colors) <- cluster.colors$cluster

  ## Generate a color pallete for the edges
  colfunc <- grDevices::colorRampPalette(c("grey", "black"))

  vertex.colors <- c()
  for (this.vertex in populations.use) {
    if (this.vertex %in% cluster.colors$cluster) {
      this.col <- as.character(cluster.colors[cluster.colors$cluster == this.vertex, ]$color)
      vertex.colors <- append(vertex.colors, this.col)
    } else{
      vertex.colors <- append(vertex.colors, "#e6e6e6")
    }
  }

  edge.colours <- c()
  for (this.source in edge.sources) {
    if (this.source %in% cluster.colors$cluster) {
      this.col <- as.character(cluster.colors[cluster.colors$cluster==this.source, ]$color)
      edge.colours <- append(edge.colours, this.col)
    } else{
      edge.colours <- append(edge.colours, "#e6e6e6")
    }
  }

  arrow.size <- arrow.size
  arrow.width <- arrow.width
  edge.multi <- edge.multi
  igraph::E(lr.plot)$width <- as.numeric(all.weights)*edge.multi

  layout = "layout_in_circle"
  l <- do.call(layout, list(lr.plot))
  par(mar=c(0,0,0,0)+.1)
  plot(lr.plot, edge.curved = 0.2, vertex.color=vertex.colors,
       layout=l, vertex.label.cex=1.1, edge.arrow.size=arrow.size, edge.arrow.width=arrow.width, vertex.size=18,
       vertex.label.font=2, vertex.label.color="black", edge.color = edge.colours)


}

#####################################################################################
### Generate a network plot linking source:ligand:receptor:targets in tree format ###
#####################################################################################
#'
#' Generate a network plot linking source:ligand:receptor:targets in tree format
#'
#' Given a file generated by the EvaluateConnections function,
#' subsets according to adjusted P-value and plots the significant
#' cell-cell connections. Uses the igraph package for network plotting.
#'
#' @param path.file input file
#' @param edge.score.file adjusted p-value threshold for defining significance
#' @param source.population named vector of colour codes to fill node colours
#' @param target.populations minimum weight for labelling populations in the plot.
#' @param ligand.set set of ligands to plot
#' @param receptor.set set of receptors to plot
#' @param population.cols option colour set for cell populations. provide in order of ligand, receptors
#' @param ligand.col colour for the row of ligands
#' @param receptor.col colour for the row of receptors
#' @return a ggplot2 object.
#' @examples
#' NetworkTreePlot(path.file = example.file)
#'
#'
#' @export
#'
#' @importFrom igraph layout_in_circle
#'
NetworkTreePlot <- function(path.file,
                            edge.score.file,
                            source.population,
                            target.populations,
                            source.marker.genes = NULL,
                            target.marker.genes = NULL,
                            ligand.set = NULL,
                            receptor.set = NULL,
                            population.cols = NULL,
                            ligand.col = "#c7d8e8",
                            receptor.col = "#d2bcff") {

  path.table <- read.csv(path.file, row.names = 1, stringsAsFactors = FALSE)

  ## if a set of receptors has been provided, subset the path table accordingly
  if (!is.null(receptor.set)) {
    path.table <- path.table[path.table$Receptor %in% receptor.set, ]
  } else if (!is.null(target.marker.genes)) {
    path.table <- path.table[path.table$Receptor %in% target.marker.genes, ]
  }

  ## if a set of ligands has been provided, subset the path table accordingly
  if (!is.null(ligand.set)) {
    path.table <- path.table[path.table$Ligand %in% ligand.set, ]
  } else if (!is.null(source.marker.genes)) {
    path.table <- path.table[path.table$Ligand %in% source.marker.genes, ]
  }

  population.labels <- union(source.population, target.populations)

  ## subset according to the source population
  source.population <- paste0("S:", source.population)
  path.table <- path.table[path.table$Source == source.population, ]

  ## If target populations have been provided, subset the table
  target.populations <- paste0("T:", target.populations)
  this.path.table <- path.table[path.table$Target %in% target.populations, ]




  #############################################################
  ### Given a path table, plot a network in a 'tree' format ###
  #############################################################

  ## Define a paths table
  path.sub.table <- this.path.table

  ## Read in score file
  score.table <- read.csv(edge.score.file)
  score.labels <- apply(score.table[, c("source", "target")], 1, function(x) paste0(x[1], ".", x[2]))
  rownames(score.table) <- score.labels

  ## Pull out the receptor examples and scale values
  receptor.foldchange.indicies <- which(score.table$relationship == "receptor.cluster")
  receptor.foldchange.table <- score.table[receptor.foldchange.indicies, ]
  receptor.foldchange.table$weight <- scales::rescale(receptor.foldchange.table$weight)

  ## Now Pull out ligand examples and scale values
  ligand.foldchange.indicies <- which(score.table$relationship == "cluster.ligand")
  ligand.foldchange.table <- score.table[ligand.foldchange.indicies, ]
  ligand.foldchange.table$weight <- scales::rescale(ligand.foldchange.table$weight)

  ligands <- path.sub.table$Ligand
  receptors <- path.sub.table$Receptor

  ## First layer of edges: origin cluster -> ligand
  layer1.edges <- c()
  for (i in c(1:nrow(path.sub.table))) {
    this.edge <- c(as.character(path.sub.table[i, 1]), as.character(path.sub.table[i, 2]))
    layer1.edges <- rbind(layer1.edges, this.edge)
  }
  layer1.edges <- unique(layer1.edges)
  layer.labels <- as.character(apply(layer1.edges, 1, function(x) paste0(x[1], ".", x[2])))
  rownames(layer1.edges) <- layer.labels

  layer1.weights <- ligand.foldchange.table[layer.labels, ]$weight

  ### Second layer of edges: ligand -> receptor
  layer2.edges <- c()
  for (i in c(1:nrow(path.sub.table))) {
    this.edge <- c(as.character(path.sub.table[i, 2]), as.character(path.sub.table[i, 3]))
    layer2.edges <- rbind(layer2.edges, this.edge)
  }
  layer2.edges <- unique(layer2.edges)
  layer.labels <- as.character(apply(layer2.edges, 1, function(x) paste0(x[1], ".", x[2])))
  rownames(layer2.edges) <- layer.labels

  layer2.weights <- score.table[layer.labels, ]$weight

  ### Third layer of edges: receptor -> target cluster
  layer3.edges <- c()
  for (i in c(1:nrow(path.sub.table))) {
    this.edge <- c(as.character(path.sub.table[i, 3]), as.character(path.sub.table[i, 4]))
    layer3.edges <- rbind(layer3.edges, this.edge)
  }
  layer3.edges <- unique(layer3.edges)
  layer.labels <- as.character(apply(layer3.edges, 1, function(x) paste0(x[1], ".", x[2])))
  rownames(layer3.edges) <- layer.labels

  layer3.weights <- receptor.foldchange.table[layer.labels, ]$weight

  ## Combined all edges together for graph input
  all.edges <- c()
  for (i in c(1:nrow(layer1.edges))){
    all.edges <- append(all.edges, layer1.edges[i,])
  }
  for (i in c(1:nrow(layer2.edges))){
    all.edges <- append(all.edges, layer2.edges[i,])
  }
  for (i in c(1:nrow(layer3.edges))){
    all.edges <- append(all.edges, layer3.edges[i,])
  }

  ## Combine all weights
  all.weights <- c(layer1.weights, layer2.weights, layer3.weights)

  #net = network(all.edges, directed = TRUE)
  lr.plot <- igraph::graph(all.edges)
  graph <- tidygraph::as_tbl_graph(lr.plot)

  named.p.labels <- population.labels
  names(named.p.labels) <- population.labels

  # Specific vertex class - either a cluster or a gene
  vertex.class <- c()
  tree.level <- c()
  node.names <- c()
  for (this.vertex in names(igraph::V(lr.plot))) {
    if (startsWith(this.vertex, "S:")) {
      this.cluster <- sub(".:(.*)", "\\1", this.vertex)
      this.cluster <- named.p.labels[this.cluster]
      vertex.class <- append(vertex.class, this.cluster)
      node.names <- append(node.names, this.cluster)
      tree.level <- append(tree.level, "Source")
    } else if (startsWith(this.vertex, "T:")) {
      this.cluster <- sub(".:(.*)", "\\1", this.vertex)
      this.cluster <- named.p.labels[this.cluster]
      vertex.class <- append(vertex.class, this.cluster)
      node.names <- append(node.names, this.cluster)
      tree.level <- append(tree.level, "Target")
    } else{
      if (this.vertex %in% ligands) {
        vertex.class <- append(vertex.class, "Ligand")
        node.names <- append(node.names, this.vertex)
        tree.level <- append(tree.level, "Ligand")
      } else if (this.vertex %in% receptors) {
        vertex.class <- append(vertex.class, "Receptor")
        node.names <- append(node.names, this.vertex)
        tree.level <- append(tree.level, "Receptor")
      } else {
        warning(paste0(this.vertex, " not in ligand or receptor list"))
      }
    }
  }
  igraph::V(lr.plot)$Class <- vertex.class
  igraph::V(lr.plot)$Label <- node.names
  igraph::V(lr.plot)$Level <- tree.level

  igraph::E(lr.plot)$Weight <- all.weights

  v.cols = c("#e6e6e6")

  # Use ggraph to plot in a tree layout
  graph <- tidygraph::as_tbl_graph(lr.plot)
  net.pl <- ggraph::ggraph(graph, 'igraph', algorithm = "tree") +
    ggraph::geom_node_point() +
    ggplot2::theme_void() +
    ggraph::geom_edge_link(arrow = arrow(20, unit(.25, "cm"), type = "closed"), end_cap = ggraph::circle(4.8, "mm"), width = 0.75) +
    ggraph::geom_node_point(aes(colour = Class, size = 5)) +
    ggplot2::scale_size(range = c(4, 14)) +
    ggraph::geom_node_text(aes(label = Label), size = 3) +
    ggplot2::theme(legend.position = "bottom") + ggplot2::guides(size = FALSE)
  net.pl <- rescale_node_coordinates(net.pl, scale.by = "Ligand")

  ## Now adjust the node colours
  if (is.null(population.cols)) {
    col.set <- scales::hue_pal()(length(population.labels))
  } else{
    col.set <- population.cols
  }

  col.table = data.frame(Population = population.labels, Col = col.set)
  rownames(col.table) = population.labels

  node.names = levels(net.pl$data$name)
  labels.names = c()
  node.cols = c()
  for (thisName in levels(factor(net.pl$data$Class))) {
    if (thisName %in% rownames(col.table)) {
      node.cols = append(node.cols, as.character(col.table[thisName, 2]))
      labels.names = append(labels.names, as.character(col.table[thisName, 1]))
    } else {
      if (thisName == "Ligand") {
        node.cols = append(node.cols, ligand.col)
        labels.names = append(labels.names, thisName)
      } else {
        node.cols = append(node.cols, receptor.col)
        labels.names = append(labels.names, thisName)
      }
    }
  }

  net.pl <- net.pl + scale_colour_manual(values = node.cols) + theme(legend.position = "right") +
    guides(colour = guide_legend(override.aes = list(size=4)))

  return(net.pl)
}

