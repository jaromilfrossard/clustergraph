#'Computer clustermass test based on multiple signal.
#'
#'@description compute the maris oostenveld permutation test of on multiple signal. Used for full scalp EEG analysis.
#'
#'@param formula the formula of the model
#'@param data a dataframe containing the design
#'@param signal a 3 dimentional array. The the row are the observations, colomn the time and the third dimention are the nodes of the graph.
#'@param method a charcter string specifying the method
#'@param P a Pmat object from \code{permuco}
#'@param graph a igraph object of an undirected graph specifing the neighborgoods relationship between the nodes
#'@param multcomp the multiple comparison procedure only \code{"maris_oostenveld"} is available.
#'@param return_distribution If set to true return the null distribution by permutation.
#'@param coding_sum logical. If \code{TRUE}, set the coding of the design to sum, if \code{FALSE}, take the coding define in the dataframe.
#'@param threshold see \code{clusterlm}.
#'@param np a scalar indicating the number of permutations. Will be overwrite by \code{P} if specified.
#'@param aggr_FUN the function to aggregate individual statistics into cluster mass.
#'@param effect a number indicating the effect to test. Refer to the \code{assign} attribute of the \code{model.matrix} object. If \code{NULL} it will test all the effects.
#'@param ... further arguments
#'@importFrom stats update as.formula contr.sum model.frame contrasts<- model.matrix qf
#'@importFrom igraph permute.vertices
#'@import permuco
#'@export
clustergraph_rnd <- function(formula, data, signal, method, threshold, np, P, graph, effect = NULL, coding_sum = T,
                               aggr_FUN =sum, multcomp = "maris_oostenveld", return_distribution,...){

    if(is.null(method)){method = "Rd_kheradPajouh_renaud"}

    switch(
      paste(method,sep="_"),
      "Rd_kheradPajouh_renaud" = {
        funP = function(...) {
          graph_fisher_Rd_kheradPajouh_renaud_rnd(...)
        }},
      {funP=function(...){eval(parse(text=paste("graph_fisher_",method,"_rnd(...)",sep="",collpase="")))}})

  switch(multcomp,
         "maris_oostenveld" = {
           funMultComp = function(distribution,threshold,aggr_FUN,graph){
             compute_maris_oostenveld_array(distribution = distribution,threshold = threshold,
                                                           aggr_FUN = aggr_FUN,graph = graph)}},
         "troendle" = {funMultComp = function(distribution,threshold,aggr_FUN,graph){
           compute_troendle_array(distribution = distribution,graph = graph, alpha = 0.05, ...)
         }},{
           multcomp = "maris_oostenveld"
           funMultComp = function(distribution,threshold,aggr_FUN,graph){
             compute_maris_oostenveld_array(distribution = distribution,threshold = threshold,
                                                           aggr_FUN = aggr_FUN,graph = graph)}
         })

    if(class(signal) != "array"){
      stop("convert signal into a 3 dimentinal array")
    }
    dotargs = list(...)


    #Formula transforamtion
    terms<-terms(formula,special="Error",data=data)
    ind_error <- attr(terms, "specials")$Error
    error_term <- attr(terms, "variables")[[1 + ind_error]]
    formula_f <- update(formula, paste("~ .-",deparse(error_term, width.cutoff = 500L, backtick = TRUE)))
    e_term <- deparse(error_term[[2L]], width.cutoff = 500L,backtick = TRUE)

    #random/fix formula

    formula_allfixed <- as.formula(paste(c("~",formula_f[[2]],"+",e_term),collapse=""))
    formula_allfixed_design <- as.formula(paste(c("~",formula_f[[2]],"+",e_term),collapse=""))

    formula_within <- formula(paste("~", e_term, collapse = ""))
    formula_within<- formula(paste("~",deparse(error_term[[2]][[3]]),collapse=""))
    formula_id <- formula(paste("~",deparse(error_term[[2]][[2]]),collapse = ""))



    #Model frame
    mf <- model.frame(formula = formula_allfixed, data = data)
    mf_design <- model.frame(formula = formula_allfixed_design, data = data)
    if(coding_sum){mf_design <- permuco:::changeContrast(mf_design, contr = contr.sum)}

    mf_f <- model.frame(formula = formula_f, data = mf_design)
    mf_id <- model.frame(formula = formula_id, data = as.data.frame(lapply(mf_design,function(col){
      col = as.factor(col)
      contrasts(col) = contr.sum
      col})))

    ##response
    mf <- eval(mf, parent.frame(n=1))
    dim_y = dim(signal)
    dnames = dimnames(signal)
    dim(signal) = c(dim_y[1],dim_y[2]*dim_y[3])

    ##link fixed random
    link = link(formula_f=formula_f,formula_within=formula_within)

    ###model .matrix
    mm_f <- model.matrix(attr(mf_f, "terms"), data = mf_f)
    mm_id <- model.matrix(attr(mf_id, "terms"), data = mf_id)[,-1,drop=F]
    name <- colnames(mm_f)

    ##checkk data


    permuco:::checkBalancedData(fixed_formula = formula_f, data = cbind(mf))

    #compute permutation
    if (is.null(P)) {P = Pmat(np = np, n = dim(signal)[1])}
    np = permuco:::np.Pmat(P)

    ##distribution
    args <- list(mm = mm_f, mm_id = mm_id, link = link, P = P, y = signal)
    args = c(args,dotargs)
    if(is.null(effect)){effect = 1:max(attr(mm_f,"assign"))}


    multiple_comparison <- list()
    length(multiple_comparison) <- length(effect)
    names(multiple_comparison) <- attr(attr(mf_f, "terms"), "term.labels")[effect]


    ##adjust multiple threshold
    if(is.null(threshold)){
      df = permuco:::compute_degree_freedom_rnd(test = "fisher",mm = mm_f,assigni = attr(mm_f,"assign"),mm_id = mm_id,link = link)
      threshold = qf(p = 0.95, df1 = df[,1],df2 =df[,2])
    }else if(length(threshold)==1){threshold = rep(threshold,length(multiple_comparison))
    } else if(length(threshold)>1){
      threshold = as.numeric(matrix(threshold,nrow=length(multiple_comparison)))
    }

    if(is.null(dnames[[3]])){
      warning("Response and graph must match vertices and 3rd dimension.")
    }else{
      perm = match(get.vertex.attribute(graph,"name"),dnames[[3]])
      if( !isTRUE(all.equal(perm,1:length(perm)))){
        warning("reorder graph vertices to match signal")
        if(sum(is.na(perm))>=1){
          stop("Names of graph and electrodes must match.")
        }
          graph = permute.vertices(graph,perm )
      }

    }

    cat("Computing Effect:\n")


    for(i in effect){
      cat(i)
      cat("\n")
      args$i = i
      ###initialisze output
      distribution = funP(args = args)
      gc()
      pvalue <- apply(distribution,2,function(col){
        permuco:::compute_pvalue(distribution = col)})
      ##inarray
      pvalue = matrix(pvalue,nrow = dim_y[2],ncol = dim_y[3])
      dim(distribution) = c(np,dim_y[2],dim_y[3])
      dim(signal) = dim_y
      mci = which(effect==i)


      multiple_comparison[[mci]]=list()
      multiple_comparison[[mci]]$uncorrected = list(statistic = t(distribution[1,,]),pvalue = pvalue)
      if(return_distribution){multiple_comparison[[mci]]$uncorrected$distribution = distribution}
      multiple_comparison[[mci]][[2]] =
        compute_maris_oostenveld_array(distribution = distribution,
                                       threshold = threshold[i], aggr_FUN = aggr_FUN, graph = graph)
      names(multiple_comparison[[mci]])[2] = multcomp
    }


    dimnames(signal) = dnames
    attr(mf,"terms")=NULL


    table = clustergraph_table(multiple_comparison)


    out=list()
    out$y = signal
    out$model.matrix = mm_f
    out$model.matrix_id = mm_id = mm_id
    out$link = link
    out$P = P
    out$multiple_comparison = multiple_comparison
    out$table = table
    out$data=mf
    out$method = method
    out$multcomp = multcomp
    out$threshold = threshold
    out$effect = effect
    out$graph = graph
    class(out) <- "clustergraph"
    return(out)

  }
