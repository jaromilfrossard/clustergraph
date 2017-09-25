# Compute the Maris Oostenveld pvalue correction
#
# @description Compute the Maris Oostenveld pvalue correction for an array
# @param distribution An 3d array representing the null distribution of multiple signal. The first dimension is the permutations, the second the time points, the third is the elecrodes.
# @param threshold The threshold used to compute the clusters.
# @param aggr_FUN The function that aggregate the cluster into a scalar (cluster mass).
# @param graph A igraph object representing the adjacency of electrod in a scalp.
# @param t A numeric representing the time point
# @return graph A igraph object with vertices attributes statistic and pvalue. data,A data frame containing the vraible elecrod, time, statistic, pvalue. cluster the result of the connected componant search on the observed graph enhanced with the mass and pvalue of the cluster. Distribution the cluster mass null distribution.
#' @importFrom igraph set_vertex_attr delete_vertices clusters get.vertex.attribute V set.vertex.attribute
compute_maris_oostenveld_array = function(distribution, threshold,aggr_FUN, graph){
  graph = full_graph(graph, t = dim(distribution)[2])
  ###null
  mass_distribution = apply(distribution,c(1),function(stat){
    gi = set_vertex_attr(graph, name="statistic", value = as.numeric(t(stat)))
    gi = delete_vertices(gi,V(gi)[get.vertex.attribute(gi, "statistic")<=threshold])
    cc = clusters(gi,mode = "weak")
    if(length(cc$membership)==0){
      return(0)}else{
        return(max(sapply(1:max(1,max(cc$membership)),function(i){
          aggr_FUN(get.vertex.attribute(gi,name = "statistic", index = names(cc$membership)[cc$membership == i]))
        })))}

  })

  ##observed
  graph = set_vertex_attr(graph, name="statistic", value = as.numeric(t(distribution[1,,])))
  g = delete_vertices(graph,V(graph)[get.vertex.attribute(graph, "statistic")<=threshold])
  cc = clusters(g,mode ="weak")

  mass_statistic = sapply(1:max(cc$membership), function(i) {
    aggr_FUN(get.vertex.attribute(g,name = "statistic", index = names(cc$membership)[cc$membership == i]))
  })
  pvalue = sapply(mass_statistic, function(mi) permuco:::compute_pvalue(stat = mi,
                                                                        distribution = mass_distribution, laterality = "right"))

  cc$mass_statistic = mass_statistic
  cc$pvalue = pvalue


  graph = set.vertex.attribute(graph,"pvalue", index = names(cc$membership),value = cc$pvalue[cc$membership])
  graph = set.vertex.attribute(graph,"custer_id", index = names(cc$membership),value = cc$membership)
  graph = set.vertex.attribute(graph,"mass_statistic", index = names(cc$membership),value = cc$mass_statistic[cc$membership])

  df = data.frame(electrode = get.vertex.attribute(graph,name=c("electrode")),
                  time =as.numeric(get.vertex.attribute(graph,name=c("time"))),
                  statistic = as.numeric(get.vertex.attribute(graph,name=c("statistic"))),
                  pvalue = as.numeric(get.vertex.attribute(graph,name=c("pvalue"))),
                  cluster_id = as.numeric(get.vertex.attribute(graph,name=c("custer_id"))),
                  mass_statistic = as.numeric(get.vertex.attribute(graph,name=c("mass_statistic"))))

  df$pvalue[is.na(df$pvalue)]=1
  df$cluster_id[is.na(df$cluster_id)]=0
  df$mass_statistic[is.na(df$mass_statistic)]=0

  return(list(graph = graph, data = df,cluster = cc, distribution = mass_distribution))
}

