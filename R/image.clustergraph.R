#' Plot a statistical map of one effect.
#'
#' @description plot the significant test for 1 effect, all timepoints and all electrodes. Electrode are in the y-axis and the time in the x-axis. Non-significant cluster are shown in grey and the significant one in color from yellow to red as a function of the individual statistics.
#' @param x a clustergraph model.
#' @param effect an integer indicating which effect to plot.
#' @param main a character string indicating the name of the graphics. The default, NULL, get the name of the effect.
#' @param ylab see par. Default is "".
#' @param xlab see par. Default is "".
#' @param ... other argument pass to image().,
#' @export
image.clustergraph = function(x, effect = 1, main = NULL,ylab = "", xlab = "",...){
  switch(names(x$multiple_comparison[[effect]])[2],
         "maris_oostenveld" = {image_maris_oostenveld(x = x, effect = effect, main = main,
                                                     ylab = ylab, xlab = xlab,...)},
         "troendle" =  image_troendle(x = x, effect = effect, main = main,
                                             ylab = ylab, xlab = xlab,...))
}






