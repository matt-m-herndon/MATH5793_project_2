#' Arrange several plots into a single view
#'
#' @param ... Set of plots to arrange.
#' @param plotlist List of plots
#' @param file Unused
#' @param cols Number of columns to arrange the plots into.
#' @param layout Layout matrix
#'
#'
#' @examples
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }

  if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#' Distribution of xbar and theta
#'
#' @param n Number of samples per iteration
#' @param iter Number of iterations to simulate
#' @param mu Expected value of the variables
#' @param sigma Covariance of the variables
#'
#'
#' @examples
xbarthetadist = function(n,iter,mu,sigma){
  library(mvtnorm)
  library(ggplot2)
  library(viridis)

  mat = matrix(NA, nr= iter, nc=3)
  colnames(mat)= c("xbar1","xbar2","theta")
  for(i in 1:iter){
    x = rmvnorm(n,mu,sigma)
    mat[i,c(1,2)] <- colMeans(x)
    s=cov(x)
    eig=eigen(s)
    theta =  acos(eig$vectors[,1][1])
    mat[i,3]<-theta
  }

  df=as.data.frame(mat)

  g = ggplot(df, aes(x=xbar1,y=xbar2))  + coord_equal()
  a = ggplot(df, aes(x=theta))

  gp = g + geom_point()
  #print(gp)
  gd = g + stat_density2d(aes(colour=..density..), geom='point', contour=F) + scale_color_viridis()
  #print(gd)

  ah = a + geom_histogram()
  ad = a + geom_density(fill="red")


  multiplot(gp, gd, ah, ad, cols=2)
  #print(ah)
  #print(ad)
  #head(mat)
}
# end of function
