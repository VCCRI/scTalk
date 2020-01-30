
#################################################################
### Generate a figure showing inbound vs outbound connections ###
#################################################################

#' Plot inbound vs outbound weights
#'
#' Generates a dot plot showing inbound vs outbound weights for all populations.
#'
#' @param path.table a table of source:ligand:receptor:target paths.
#' @param this.population a cluster ID contained in the path table.
#' @param col.set optional paramater to set colours of the populations.
#' @param lab.weight.thresh minimum weight for labelling populations in the plot.
#' @return a ggplot2 object.
#' @examples
#' PlotTopLigands(path.table = example.table, this.population = '1')
#'
#' @export
#'
InboundOutboundPlot <- function(sig.paths, this.population,
                                col.set = NULL, p.labels = NULL, lab.weight.thresh = 1000) {
  col.set = hue_pal()(length(table(tip.aggr@ident)))
  p.labels = names(table(tip.aggr@ident))

  edge.sources = factor(all.edges$Source, levels = p.labels)
  edge.targets = factor(all.edges$Target, levels = p.labels)

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
    theme(legend.position = c(0.7, 0.75)) + guides(color = guide_legend(override.aes = list(size=4), ncol = 3))
  return(pl)
}

#########################################################################################
### For a population generate a bar plot showing ligands contributing greatest weight ###
#########################################################################################

#' Plot a bar graph of the top ligands
#'
#' @param path.table a table of source:ligand:receptor:target paths.
#' @param this.population a cluster ID contained in the path table.
#' @param col.use optional paramater to set colour of the bar graph.
#' @param min.weight minimum weight for including paths in the top ligand calculations.
#' @return a ggplot2 object.
#' @examples
#' PlotTopLigands(path.table = example.table, this.population = '1')
#'
#' @export
#'
PlotTopLigands <- function(path.table, this.population, col.use = NULL, min.weight = 2,
                           max.plot.num = 15) {
  p.labels = unique(as.character(tip.aggr@ident))
  col.set = scales::hue_pal()(length(p.labels))

  path.file = paste0(out.lab, "_network_paths_weight1.5.csv")
  path.table = read.csv(path.file, sep = ",", row.names = 1)

  source.population = paste0("S:", this.population)
  this.col = col.set[which(p.labels == this.population)]
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

  pl <- ggplot(ligand.weights.table, aes(x=reorder(Ligand, -Weight), y=Weight)) + geom_bar(stat="identity") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + xlab("Ligand") + ylab("Summed weight") +
    ggtitle(paste0("Summed weights for top outbound cluster ", this.population, " ligands"))

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

  sig.paths = read.csv(path.sig.file, sep = ",")
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

  populations.use <- unique(as.character(igraph::V(lr.plot)))
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




