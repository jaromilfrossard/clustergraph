# Compute the troendle multiple comparison procedre
#
#
# @description Compute the Maris Oostenveld pvalue correction for an array
# @param distribution An 3d array representing the null distribution of multiple signal. The first dimension is the permutations, the second the time points, the third is the elecrodes.
# @param graph A igraph object representing the adjacency of electrod in a scalp.
# @return graph A igraph object with vertices attributes statistic and pvalue. data,A data frame containing the vraible elecrod, time, statistic, pvalue. cluster the result of the connected componant search on the observed graph enhanced with the mass and pvalue of the cluster. Distribution the cluster mass null distribution.
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom grDevices gray
compute_troendle_array = function(distribution,graph,alpha = 0.05, ...){
  distribution_mat = matrix(distribution,nrow=dim(distribution)[1],
                            ncol=prod(dim(distribution)[-1]))

  distribution_rank = apply(-distribution_mat,2,function(coli){
    ceiling(rank(coli))
  })
  rank_uncorr <- rank(distribution_rank[1,],ties.method = "min")
  order_test = integer(0)
  p_corrected <- numeric(length = 0)

  urank = sort(unique(rank_uncorr))[1]
  ntest = length(sort(unique(rank_uncorr)))
  pb = txtProgressBar(min = 1, max = ntest, style = 3)
  testi = 0
  for(urank in sort(unique(rank_uncorr))){
    testi=testi+1
    setTxtProgressBar(pb, testi)
    which_test <- which(urank==rank_uncorr)
    order_test <- c(order_test,which_test)
    pvali <- distribution_rank[,which(urank<=rank_uncorr),drop=F]
    distr_min <- apply(pvali,1,min)
    p_corrected <- c(p_corrected,
                     permuco:::compute_pvalue(distribution = distr_min,
                                              stat = matrix(distribution_rank[,which_test],nrow=1),
                                              alternative = "less"))
  }
  p_corrected = cummax(p_corrected)[order(order_test)]
  p_corrected = matrix(p_corrected,nrow= dim(distribution)[2],ncol= dim(distribution)[3])
  ### in graph
  graph = full_graph(graph, t = dim(distribution)[2])

  graph = set_vertex_attr(graph, name = "statistic",
                          value = as.numeric(t(distribution[1, , ])))

  graph = set_vertex_attr(graph, name = "pvalue",
                          value = as.numeric(t(p_corrected)))

  g = delete_vertices(graph, V(graph)[get.vertex.attribute(graph,
                                                           "pvalue") >= alpha])
  cc = clusters(g, mode = "weak")
  cc$mass_statistic = rep(NA,length(cc$csize))
  cc$pvalue = rep(paste0("<",alpha),length(cc$csize))




  graph = set.vertex.attribute(graph, "custer_id",
                               value = NA)

  graph = set.vertex.attribute(graph, "mass_statistic",
                               value = NA)

  graph = set.vertex.attribute(graph, "custer_id", index = names(cc$membership),
                               value = cc$membership)

  df = data.frame(electrode = get.vertex.attribute(graph, name = c("electrode")),
                  time = as.numeric(get.vertex.attribute(graph, name = c("time"))),
                  statistic = as.numeric(get.vertex.attribute(graph, name = c("statistic"))),
                  pvalue = as.numeric(get.vertex.attribute(graph, name = c("pvalue"))),
                  cluster_id = as.numeric(get.vertex.attribute(graph, name = c("custer_id"))),
                  mass_statistic = as.numeric(get.vertex.attribute(graph, name = c("mass_statistic"))))

  df$cluster_id[is.na(df$cluster_id)] = 0
  df$mass_statistic[is.na(df$mass_statistic)] = 0
  return(list(graph = graph, data = df, cluster = cc, distribution = NA,
              threshold = NA))


}