clustergraph_table = function(x){
  print(length(x))
  ct = lapply(1:length(x),function(i){
    effect = x[[i]]
    tab= data.frame(size = effect$maris_oostenveld$cluster$csize,
                    mass = effect$maris_oostenveld$cluster$mass_statistic,
                    pvalue = effect$maris_oostenveld$cluster$pvalue)
    attr(tab,"threshold") = effect$maris_oostenveld$threshold
    attr(tab,"effect_name") = names(x)[i]
    class(tab) = append("cluster_table",class(tab))
    tab
  })
  class(ct) = append("listof_cluster_table",class(ct))
  names(ct) = names(x)
  ct
}