#' Converts the pcalgo graph object into an igraph object and plots it.
#'
#' Constructs igraph object from graphnel type and adds weights (if passed in).
#' Default layout is in a circle so graphs are easily comparable.
#'
#' @param g A graph object of class graph
#' @param weights A list of weights for edges in the graph
#' @param layout A layout to organize the graph with
#' @param vertex.size A node size for the graph
#' @param vertex.color Specify the color of nodes in the graph
#' @return
igplot <- function(g, weights=TRUE, layout=igraph::layout_in_circle,
                   vertex.size=50, vertex.color="lightblue",...){

  g <- igraph::graph_from_graphnel(as(g, "graphNEL"))
  op <- par(mar=c(1,1,1,1))

  if (weights == TRUE){
    ew <- round(igraph::get.edge.attribute(g,"weight"),2)
    igraph::plot.igraph(g,layout=layout,edge.label=ew,edge.label.cex=0.9, vertex.size=vertex.size,vertex.color=vertex.color,...)
  }
  else {
    igraph::plot.igraph(g,layout=layout,vertex.size=vertex.size,vertex.color=vertex.color,...)
  }

  par(op) # Reset to previous settings
}

#' Calculate the graph of variable relationships at a particular quantile level.
#'
#' It assumes the input dataframe has at least two columns (columns are variables) and that
#' the quantile level is between 0 and 1 (non-inclusive). Additionally it assumes all
#' column names are unique.
#'
#' This function calls the pc algorithm with the quantile association hypothesis test.
#' Afterwards, it constructs the marginal associations of all variables and uses their
#' marginal associations as edge weights in the returned graph.
#'
#' An igraph plot is plotted with the marginal associations at the quantile level
#' prescribed.
#'
#' @param data A data frame of data with unique column/rownames
#' @param tau A particular quantile level (0 to 1, not inclusive)
#' @return
#' @examples
#' plot.quantile.dag(df, tau=0.5)
plot.quantile.dag <- function(data, tau, verbose=FALSE) {
  if(tau <= 0 | tau >= 1) {
    return("Tau must be within 0 and 1 (non-inclusive)")
  }

  saveRDS(tau, "tau.rds")
  data <- jitter.columns(data, factor=0.1)
  pc_graph <- pcalg::pc(data, indepTest = quantile.ztest, labels = colnames(data), alpha = 0.05, verbose = verbose)

  lookup_table <- marginal.test(data, tau=tau)

  colnames(lookup_table) <- colnames(data)
  rownames(lookup_table) <- colnames(data)

  data_var = pc_graph@graph@edgeData@data
  edgenames = names(data_var)

  for(i in 1:length(edgenames)) {
    curredge = edgenames[i]
    varnames = strsplit(curredge, split="|", fixed=TRUE)
    weightval = lookup_table[ varnames[[1]][1], varnames[[1]][2] ]
    data_var[[ curredge ]]$weight = weightval
  }

  pc_graph@graph@edgeData@data = data_var
  unlink("tau.rds")

  igplot(pc_graph@graph, weights=TRUE)
}
