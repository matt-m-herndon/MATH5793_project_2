library(future)
library(rprojroot)

root <- rprojroot::is_r_package
qq_critical_points <- future(
  read.csv(root$find_file("table_4_2_qq_critical.csv"))
)

#readLines(root$find_file("table_4_2_qq_critical.csv"), 3)

assess_normality = function(x, significance=5){
  # number of samples in variable
  n = length(x)
  
  # sample normal quantile function
  q = qnorm((seq(1,n)-1/2)/n)
  
  # mapping between q and sorted x
  mp = q ~ sort(x)
  
  # find linear model to dataset
  fit = lm(mp)
  
  return(fit)
}


#' Given a dataframe, assess normality of each numeric variable. 
#' When a variable exhibits non-normal behavior, attempt to
#' transform to normality. 
#' 
#' @param data Dataset to operate on
#'
#' @examples
normal = function(dataset)
{
  stopifnot(!missing(dataset))
  stopifnot(class(dataset)=="data.frame")
  stopifnot(length(dim(dataset)) < 3)
  
  # extract portion of dataset which is numeric
  dset_numeric = dataset[,sapply(dataset, function(x) is.numeric(x))]
  
  tn = assess_normality(dset_numeric[,1])
  
  obj = list(original=dataset, transformed_normal=tn)
  # obj is now a list of class `list` and `normal`
  class(obj) = append(class(obj),"normal") 
  return(obj)
}


summary.normal <- function (object, ...) {
  summary(object$dataset)
}


print.normal <- function (object, ...) {
  print(paste(
      sprintf('Dataset has %i variables: ',length(variable.names(object$dataset))),
      paste(variable.names(iris),collapse='; ')
  ))
}


plot.normal <- function (object, ...) {

}


# load iris dataset
data('iris')

# just use one species
setosa = subset(iris, Species == "setosa")[1:4]

iris_red  = normal(dataset=setosa)

summary(iris_red)
print(iris_red)
plot(iris_red)

