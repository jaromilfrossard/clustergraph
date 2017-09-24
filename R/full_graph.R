# Compute time dependant graph
#
# @description Extend a graph of a scalp to t time point
# @param graph A igraph object representing the adjacency of electrod in a scalp.
# @param t A numeric representing the time point
# @return An igraph object with t \times V(graph) vertices.
#' @importFrom igraph as_adjacency_matrix graph_from_adjacency_matrix
#' @importFrom Matrix bdiag Diagonal Matrix cBind
full_graph = function(graph, t){
  full_graph = as_adjacency_matrix(graph, type = c("upper"))
  nv = ncol(full_graph)
  names_df = expand.grid(rownames(full_graph),1:t)
  names = paste(names_df[,1],"_",names_df[,2],sep="")
  full_graph = lapply(1:t,function(i)full_graph)
  full_graph = bdiag(full_graph)

  full_graph = full_graph +rbind(cBind(Matrix(0,nrow = nv*(t-1),ncol = nv),Diagonal(nv*(t-1))),
                                 cBind(Matrix(0,ncol = nv*(t),nrow = nv)))
  rownames(full_graph) = names
  colnames(full_graph) = names
  full_graph = graph_from_adjacency_matrix(full_graph,mode="upper")
  full_graph = set.vertex.attribute(full_graph,name = "electrode", value = as.character(names_df[,1]))
  full_graph = set.vertex.attribute(full_graph,name = "time", value = as.numeric(as.character(names_df[,2])))

  full_graph
}
