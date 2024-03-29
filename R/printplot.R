#' Print method for proLikelihood object
#'
#' Print \sQuote{poi} (mean difference, total standard deviation ratio or within-subject standard deviation ratio) and \sQuote{maxLik} (corresponding profile likelihood) elements from a proLikelihood object.
#'
#' @param x proLikelihood object
#' @param \dots unused
#'
#' @return Output the mean difference, total standard deviation ratio or within-subject standard deviation ratio values (depending on the \sQuote{method}) with 
#' the calculated corresponding profile likelihood values.
#'
#' @examples
#' \donttest{
#' data(dat)
#' cols <- list(subject = 'subject', formula = 'formula', y = 'AUC')
#' l <- averageBE(dat, colSpec = cols, xlength = 300)
#' l
#' }
#'
#' @export

print.proLikelihood <- function(x, ...) { 
  print(cbind(x$poi, x$maxLik))
  invisible(x)
}

#' Plot method for proLikelihood object
#'
#' This function generates a plot of a standardized profile likelihood after running the proLikelihood() function.
#'
#' The function generates a plot of the standardized profile likelihood (the profile likelihood relative to the maximum) with
#' the maximum likelihood estimate and 1/8 and 1/32 likelihood intervals for the parameter of interest (mean difference, 
#' total standard deviation ratio or within-subject standard deviation ratio depending on the \sQuote{method}) printed inside the plot.
#'
#' @param x proLikelihood object
#' @param textx numeric value, position (x-axis) of label for the maximum likelihood estimate and the 1/8 and 1/32 likelihood intervals. 
#' @param texty numeric value, position (y-axis) of label the maximum likelihood estimate and the 1/8 and 1/32 likelihood intervals. 
#' @param textsize numeric value text size of the label.
#' @param \dots unused
#'
#' @return ggplot2 object, a plot of the standardized profile likelihood with the maximum likelihood estimate and 1/8 and 1/32
#' likelihood intervals printed inside the plot.
#'
#' @examples
#' \donttest{
#' data(dat)
#' cols <- list(subject = 'subject', formula = 'formula', y = 'AUC')
#' p4a <- averageBE(dat, colSpec = cols, xlength = 50)
#' p4t <- totalVarianceBE(dat, colSpec = cols, xlength = 50)
#' p4w <- withinVarianceBE(dat, colSpec = cols, xlength = 50)
#' plot(p4a)
#' plot(p4t)
#' plot(p4w)
#' # three period case
#' dd3 <- dat[dat$period < 4,]
#' p3a <- averageBE(dd3, colSpec = cols, xlength = 50)
#' plot(p3a)
#' # two period case
#' dd2 <- dat[dat$period < 3,]
#' p2a <- averageBE(dd2, colSpec = cols, xlength = 50)
#' plot(p2a)
#' }
#'
#' @export

plot.proLikelihood <- function(x, textx, texty = 0.9, textsize = 3, ...) {
  a6 <- round(x$MAX, 3)
  b6 <- round(x$LI['1/8 LI','lower'], 3)
  c6 <- round(x$LI['1/8 LI','upper'], 3)
  d6 <- round(x$LI['1/32 LI','lower'], 3)
  e6 <- round(x$LI['1/32 LI','upper'], 3)

  ## make a data frame##
  lik.norm <- x$maxLik / max(x$maxLik, na.rm = TRUE)
  profile <- data.frame(poi = x$poi, lik.norm)

  if(x$method == 'average') {
    vlinelow <- -0.223
    vlineup <- 0.223
    xlabel <- expression(mu[T]-mu[R])
    deftx <- 0
  } else if(x$method == 'total') {
    vlinelow <- 0.4
    vlineup <- 2.5
    xlabel <- bquote(sigma[TT]/sigma[TR])
    deftx <- 2
  } else if(x$method == 'within') {
    vlinelow <- 0.4
    vlineup <- 2.5
    xlabel <- bquote(sigma[WT]/sigma[WR])
    deftx <- 2
  }

  poi <- NA # avoid CHK warning
  if(missing(textx)) {
    textx <- deftx
  }
  l1 <- sprintf('Max at %s', a6)
  l2 <- sprintf('1/8 LI (%s,%s)', b6, c6)
  l3 <- sprintf('1/32 LI (%s,%s)', d6, e6)
  maxlabel <- paste(l1, l2, l3, sep = '\n')
  p6 <- ggplot(data = profile, aes(poi, lik.norm), colour="black")
  p6 <- p6 + geom_line(size = 0.2) +
    geom_text(x = textx, y = texty, label = maxlabel, size = textsize) +
    geom_segment(aes(x = b6, y = 1/8, xend = c6, yend = 1/8), size = 0.2) +
    geom_segment(aes(x = d6, y = 1/32, xend = e6, yend = 1/32), size = 0.2) +
    geom_vline(xintercept = c(vlinelow, vlineup), linetype = 2, size = 0.2) +
    ylab("Standardized profile likelihood") +
    xlab(xlabel)
  print(p6)
}
