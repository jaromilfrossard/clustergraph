#' @importFrom graphics axis box image
#' @importFrom tidyr spread
#' @importFrom igraph vertex_attr
image_maris_oostenveld = function(x, effect = 1, main = NULL,ylab = "", xlab = "",...){


  if(is.null(main)){
    main = names(x$multiple_comparison)[effect]
  }


  ## order from electrode position
  order_electrode = order(vertex_attr(x$graph,"x"),-vertex_attr(x$graph,"y"))
  enames = vertex_attr(x$graph,"name")[order_electrode]

  ## non corrected fvalue
  fvalue=x$multiple_comparison[[effect]]$uncorrected$statistic[order_electrode,]
  data=x$multiple_comparison[[effect]]$maris_oostenveld$data
  data = data[data$cluster_id!=0,]

  ## transform pvalue into matrix electrode-time
  na_mat = matrix(NA,ncol=ncol(fvalue),nrow=nrow(fvalue))
  pvalue = spread(data[,c(1,2,4)],"time","pvalue")
  pvalue = pvalue[order(match(as.character(pvalue$electrode),vertex_attr(x$graph,"name"))),]
  enames_pv = as.character(pvalue$electrode)
  pvalue = as.matrix(pvalue[,-1])
  enames_missing = enames[which(!(enames%in%enames_pv))]
  tnames_pv = as.numeric(colnames(pvalue))

  pvalue = rbind(pvalue,matrix(NA,nrow=length(enames_missing),ncol=dim(pvalue)[2]))
  pvalue = pvalue[match(enames,c(enames_pv,enames_missing)),]


  na_mat[,tnames_pv]<-pvalue
  pvalue = na_mat


  ## delete not in cluster data
  fvalue[is.na(pvalue)] = NA

  ## non significant time point in grey
  nsfvalue = fvalue
  nsfvalue[pvalue<0.05]=NA
  nsfvalue[pvalue>=0.05]=1

  ## plot only significant in color
  fvalue[pvalue>=0.05]=NA

  xlim= c(1,nrow(fvalue))
  ylim= c(1,ncol(fvalue))


  #rev matrix and plot
  fvalue = t(apply(fvalue, 2, rev))
  nsfvalue = t(apply(nsfvalue, 2, rev))

  image(1:nrow(fvalue),1:ncol(fvalue), nsfvalue,yaxt="n",col="grey", xlab=xlab, main=main,ylab=ylab,...)
  if(sum(!is.na(fvalue))!=0){image(1:nrow(fvalue),1:ncol(fvalue), fvalue,add=T)}
  box()
  axis(side=2,at = c(xlim[1]:xlim[2])[c(xlim[1]:xlim[2])%% 2 == 0],
       labels = rev(enames[c(xlim[1]:xlim[2])%% 2 == 0]))
  axis(side=4,at = c(xlim[1]:xlim[2])[c(xlim[1]:xlim[2])%% 2 != 0],
       labels = rev(enames[c(xlim[1]:xlim[2])%% 2 != 0]))
}