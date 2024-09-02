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
#' @return Nothing
qgm.igraph.plot <- function(g, weights=TRUE, layout=igraph::layout_in_circle,
                   vertex.size=30, vertex.color="lightblue",...){

  g <- igraph::graph_from_graphnel(as(g, "graphNEL"))
  op <- par(mar=c(1,1,1,1))

  if (weights == TRUE){
    ew <- round(igraph::get.edge.attribute(g,"weight"),2)
    igraph::plot.igraph(g,layout=layout,edge.label=ew,vertex.size=vertex.size,vertex.color=vertex.color,...)
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
#' @param weights Weights will either be "marginal" (marginal association)
#' or "conditional" (association conditioned on everything) - default is
#' "marginal" (otherwise they will not be used)
#' @param quacc If you want to calculate the linear QuACC statistic versus the general no train/test statistic
#' @return Nothing
#' @examples plot.quantile.dag(df, tau=0.5)
plot.quantile.dag <- function(data, tau, m.max=Inf, weights="marginal", verbose=FALSE, quacc=TRUE, linear=TRUE, correl=FALSE, rho=FALSE, alpha=0.05) {
  if(tau <= 0 | tau >= 1) {
    return("Tau must be within 0 and 1 (non-inclusive)")
  }

  #data <- jitter.columns(data, factor=0.1)
  pc_graph <- calculate.skeleton(data, tau, m.max=m.max, quacc=quacc, linear=linear, correl=correl, verbose=verbose, adj_vector=FALSE, alpha=alpha)

  n <- length((colnames(data)))
  lookup_table <- matrix(0, nrow=n, ncol=n)
  use.weights <- tolower(weights) %in% c("conditional", "marginal")
  if(use.weights) {
    lookup_table <- pairwise.test(data, tau=tau, weights=tolower(weights), quacc=quacc, rho=rho)
  }

  colnames(lookup_table) <- colnames(data)
  rownames(lookup_table) <- colnames(data)

  data_var = pc_graph@graph@edgeData@data
  edgenames = names(data_var)

  if(is.null(edgenames)) {
    qgm.igraph.plot(pc_graph@graph, weights=FALSE)
  }
  else {
    for(i in 1:length(edgenames)) {
      curredge = edgenames[i]
      varnames = strsplit(curredge, split="|", fixed=TRUE)
      weightval = lookup_table[ varnames[[1]][1], varnames[[1]][2] ]
      data_var[[ curredge ]]$weight = weightval
    }

    pc_graph@graph@edgeData@data = data_var
    qgm.igraph.plot(pc_graph@graph, weights=use.weights)
  }
  return(pc_graph@graph)
}


#' Calculates the table of marginal relationships at a particular quantile level.
#'
#' It assumes the input dataframe has at least two columns (columns are variables) and that
#' the quantile level is between 0 and 1 (non-inclusive). Additionally it assumes all
#' column names are unique.
#'
#' This function uses calculations from the quantile concordance statistic
#' in pairwise.test and plots it in a nicer table format.
#'
#' If weights is equal to "cor" the pairwise Pearson correlation table is
#' returned.
#'
#' @param data A data frame of data with unique column/rownames
#' @param tau A particular quantile level (0 to 1, not inclusive)
#' @param weights Weights will either be "marginal" (marginal association)
#' or "conditional" (association conditioned on everything) - default is
#' "marginal" (otherwise they will not be used)
#' @param quacc If you want to calculate the linear QuACC statistic versus the standard no train/test split statistic
#' @return Nothing
#' @examples plot.marginal.table(df, tau=0.5, weights="marginal")
plot.marginal.table <- function(df, tau, weights="marginal", quacc=TRUE, rho=FALSE) {
  if (weights == "cor") {
    table.data <- cor(df)
  }
  else{
    table.data <- pairwise.test(df, tau=tau, weights=weights, quacc=quacc, rho=rho)
  }

  data <- expand.grid(X=colnames(df), Y=colnames(df))
  data$P <- as.vector(table.data)
  ggplot2::ggplot(data, ggplot2::aes(X, Y, fill=P)) +
    ggplot2::geom_tile(show.legend = FALSE) +
    ggplot2::geom_text(ggplot2::aes(label=round(P, digits=2)), size=4) +
    ggplot2::scale_fill_gradientn(limits=c(-1, 1), colors=c("blue", "white", "red")) +
    ggplot2::coord_equal() +
    ggplot2::theme(axis.title.x=element_blank(),
                   axis.title.y=element_blank(),
                   axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   axis.text.y = ggplot2::element_text(angle = 45, hjust = 1)
    )
}
