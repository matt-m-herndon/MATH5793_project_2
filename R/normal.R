library(future)
library(rprojroot)
library(pracma)

root <- rprojroot::is_r_package
# Table 4.2 from the book
qq_critical_points = read.csv(root$find_file("table_4_2_qq_critical.csv"))


boxcox = function(x,lambda){
  x = matrix(x)
  lambda = matrix(lambda,byrow=TRUE)
  print(dim(x))
  print(dim(lambda))
  # use box cox transformation
  x_l = lambda
  x_l[lambda==0] = log(x[lambda==0])
  x_l[lambda!=0] = (x[lambda!=0]^(lambda[lambda!=0])-1)/lambda[lambda!=0]
  return(x_l)
}

box_cox_maximization = function(x,lambda){
  # box-cox likelihood function, target for maximization
  n = length(x)
  m = length(lambda)
  xs = sum(log(x))
  x = repmat(matrix(x,nrow=n,nc=1),1,m)
  x_to_the_lambda = boxcox(x,repmat(matrix(lambda,nrow=1,nc=m),n,1))
  l = -n/2 * log(1/n * colSums((x_to_the_lambda-mean(x_to_the_lambda))^2)) + (lambda-1)*xs
  return(l)
}

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
  
  opt_interval = c(-5,5)
  lambda_opt = optimize(
    function(lambda) box_cox_maximization(iris[,1],lambda), 
    c(-5,5),
  )
  print(lambda_opt)
  
  obj = list(original=dataset, transformed_normal=tn, lambda=lambda)
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

