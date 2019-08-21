clustergraph_table = function(x){
  ct = lapply(1:length(x),function(i){
    effect = x[[i]]
    if(is.null(x[[i]])){append_to_name = ": effect not tested"}else{append_to_name = ""}
    tab= data.frame(size = effect[[2]]$cluster$csize,
                    mass = effect[[2]]$cluster$mass_statistic,
                    pvalue = effect[[2]]$cluster$pvalue)
    attr(tab,"threshold") = effect[[2]]$threshold
    attr(tab,"effect_name") = paste0(names(x)[i],append_to_name)
    class(tab) = append("cluster_table",class(tab))
    tab
  })
  class(ct) = append("listof_cluster_table",class(ct))
  names(ct) = names(x)
  ct
}