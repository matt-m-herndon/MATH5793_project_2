library(rprojroot)
library(pracma)
library(abind)
library(akima)
library(ggplot2)

root <- rprojroot::is_r_package
# Table 4.2 from the book
qq_critical_points = read.csv(root$find_file("table_4_2_qq_critical.csv"))


#' Spline generator along significance dimension---was going to use 2D splines,
#' but R wouldn't divulge its secrets. 
#' 
#' @param n Number of samples. In the table, the rows go from 
#' 
#' @examples 
#' # Number of samples
#' N = 10 
#' # Create generator function
#' q = qqspline(N)
#' # Create a set of critical points by sampling from the generator 
#' critical_points = q(seq(1,10,0.1))
#' 
qqspline = function(n){
  n_index = which(qq_critical_points[,1]>n)
  stopifnot(length(n_index)>0)
  return(splinefun(x=c(0.01,0.05,0.1), y=as.matrix(qq_critical_points[n_index[1],2:4]), method="fmm",  ties=mean))
}


#' Apply the boxcox transformation to a dataset for one or a series of lambdas. 
#' 
#' @param x One dimensional dataset. Must be coercible into a matrix.
#' @param lambda One dimensional container of lambdas. Must be coercible into a matrix.
#' @examples
#' 
#' # Load iris dataset 
#' data('iris')
#' 
#' # apply transformation to "Sepal.Length" for the Species "setosa"
#' transformed = boxcox(iris[1:25,1], 1.25)
#' 
#' # Now try it again for a series of lambdas
#' transformed = boxcox(iris[1:25,1], c(1.23, 1.25, 2.53))
#' 
boxcox = function(x,lambda){
  if (!is.matrix(x)) {
    x = matrix(x)
  }
  if (!is.matrix(lambda)) {
    lambda = matrix(lambda)
  }
  stopifnot(length(dim(x))<=2, length(dim(lambda))<=2)
  stopifnot(min(dim(x))==1, min(dim(lambda))==1)
  
  # generate mxn matrix
  m = prod(dim(x))
  n = prod(dim(lambda))
  x = repmat(matrix(x,nr=m,nc=1),1,n)
  lambda = repmat(matrix(lambda,nr=1,nc=n), m, 1)
  
  stopifnot(all(dim(x) == dim(lambda)))
  
  # use box cox transformation; mask to select which operation to perform
  x_l = x
  x_l[lambda==0] = log(x[lambda==0])
  x_l[lambda!=0] = (x[lambda!=0]^(lambda[lambda!=0])-1)/lambda[lambda!=0]
  return(x_l)
}


#' Boxcox likelihood function. Used internally to determine the best-performing
#' transformation from a search space of lambdas.
#' 
#' @param x One dimensional dataset to used for likelihood calculation
#' @param lambda One dimensional container of lambdas for testing.
#' 
#' @examples
#' 
#' # Load iris dataset 
#' data('iris')
#'
#' # Plot likelihood across a range of lambdas
#' lambdas = seq(-5,5,0.025)
#' plot(lambdas, box_cox_likelihood(iris[1:25,1], lambdas)
#' 
box_cox_likelihood = function(x,lambda){
  # box-cox likelihood function, target for maximization
  n = length(x)
  xs = sum(log(x))
  x_to_the_lambda = boxcox(x,lambda)
  l = -n/2 * log(1/n * colSums((x_to_the_lambda-mean(x_to_the_lambda))^2)) + (lambda-1)*xs
  return(l)
}


#' By using the QQ correlation test for normallity, attempt to determine 
#' whether or not a one dimensional dataset is normal within a specified
#' significance.  
#' 
#' @param x One dimensional dataset with which to apply the test.
#' @param significance Significance level percentage. Valid for any  
#' continuous value between 1 and 10 percent.
#' 
#' @examples
#' 
#' Load iris dataset 
#' data('iris')
#'
#' # Plot likelihood across a range of lambdas
#' lambdas = seq(-5,5,0.025)
#' assess_normality
#' 
_normal_assessment = function(x, significance=5){
  stopifnot(significance>=1,significance<=10)
  
  # number of samples in variable
  n = length(x)
  
  # We cannot safely proceed if sample size is too small
  stopifnot(n>=5)
  
  # If N exceeds the table's list of known critical points, just
  # use the maximum. The test will just become stronger.
  if (n >= (300)){
    # get spline for determining appropriate significance
    qq_crit = qqspline(299)
  }else{
    # get spline for determining appropriate significance
    qq_crit = qqspline(n)
  }
  
  # Estimate critical point
  critical_point = qq_crit(significance/100)
  
  # sort x values 
  sorted = sort(x)
  
  # sample normal quantile function
  q = qnorm((seq(1,n)-1/2)/n)
  
  # fit linear model
  model = lm(sorted~q)
  
  # model sampled using q points
  model_samples=model$coefficients[1] + q * model$coefficients[2]
  
  # correlate model with measurements
  correlation = cor(sorted, model_samples)
  
  obj = list(
    x=x,
    isnormal=correlation>=critical_point,
    model=model,
    sorted_samples=sorted,
    quartile_samples=q,
    model_samples=model_samples,
    correlation=correlation,
    critical_point=critical_point,
    significance=significance
  )
  # obj is now a list of class `list` and `normal`
  class(obj) = append(class(obj),"normal_assessment") 
  return(obj)
}

summary.normal_assessment <- function (object, ...) {
  nms = names(object$x)
  if (is.null(nms)){
    name = 'V1'
  }else{
    name = nms[[1]]
  }
  str = ''
  str = str + paste(sprintf('Variable: %s', name), collapse="\n")
  str = str + paste(summary(object$x), collapse="\n")
  # Initial test of normality
  
  
  cat(str)
}

# print.normal <- function (object, ...) {
#   print(paste(
#       sprintf('Dataset has %i variables: ',length(variable.names(object$dataset))),
#       paste(variable.names(iris),collapse='; ')
#   ))
# }
# 
# 
# plot.normal <- function (object, ...) {
#   assessments = object$assessments
#   first = assessments[[1]]
#   
#   if (all(is.na(first$transformed))){
#     assessment = first$assessment
#   }else{
#     lambda = first$transformed$lambda
#     print(sprintf('Transformed with lambda=%.2f',lambda))
#     assessment = first$transformed$assessment
#   }
#   
#   # build data frame
#   df <- data.frame(
#     q=assessment$quartile_samples,
#     sorted=assessment$sorted_sample,
#     predicted=predict(assessment$model),
#     residuals=residuals(assessment$model)
#   )
#   
#   model_a=assessment$model$coefficients[[1]]
#   model_b=assessment$model$coefficients[[2]]
#   
#   # QQ Plot showing residuals 
#   #
#   # Plot thanks to this fantastic tutorial:
#   # 
#   # https://drsimonj.svbtle.com/visualising-residuals
#   #
#   ggplot(df, aes(q, sorted)) +
#   geom_segment(aes(xend = q, yend = predicted), alpha = 0.4) +
#   geom_point(aes(color = residuals)) +
#   scale_color_gradient2(low = "blue", mid = "white", high = "red") +
#   guides(color = FALSE) + theme_bw() + 
#   geom_abline(intercept = model_a, slope = model_b, linetype="dashed", size=0.15)
# }







#' Given a dataframe, assess normality of each numeric variable. 
#' When a variable exhibits non-normal behavior, attempt to
#' transform to normality. 
#' 
#' @param dataset Validate normality of each variable within this dataset
#' @param significance Will attempt to verify each variable's normality to within 
#' this significance
#' @param boxcox_lambda_search When the significance test fails, will attempt to 
#' find the lambda from within this range which most effectively transforms 
#' the data to near-normality
#' 
#' @examples
normal = function(dataset, significance=5, boxcox_lambda_search=seq(-5,5,0.05))
{
  stopifnot(!missing(dataset))
  stopifnot(class(dataset)=="data.frame")
  stopifnot(length(dim(dataset)) < 3)
  stopifnot(length(dim(dataset)) > 0)
  
  # extract portion of dataset which is numeric
  dset_numeric = dataset[,sapply(dataset, function(x) is.numeric(x)),drop=FALSE]
  if (!all(dim(dset_numeric) != dim(dataset))){
    print('Ignoring non-numeric variables...')
  }
  
  n = dim(dset_numeric)[1]
  p = dim(dset_numeric)[2]
  print(sprintf('Checking normality of N=%i samples of the variables [%s]...', n, paste(names(dset_numeric), collapse = ', ')))
  
  assessments <- vector("list", p)
  for (i in 1:p){
    name = names(dset_numeric)[i]
    x = dset_numeric[,i]
    
    # Initial test of normality
    assessment=normal_assessment(x, significance=significance)
    
    if(assessment$isnormal){
      print(sprintf('- Variable %s is normal to within %.2f%% significance!', name, significance))
      transformed_x = NA
    }else{
      print(sprintf('- Variable %s failed the %.2f%% significance test!', name, significance))
      # Estimate lambda for boxcox transformation to attempt a transform to normality
      print('Attempting to find a boxcox transformation to normalize the dataset...')
      lambda = boxcox_lambda_search[which.max(
        box_cox_likelihood(x, boxcox_lambda_search)
      )]
      # could not get R's optimization tools to work, so just call this "close enough".
      # It's an estimate anyway.
      d = boxcox(x, lambda)
      transformed_x = list(
        data=d,
        lambda=lambda,
        assessment=normal_assessment(d)
      )
    }
    
    # store for later
    assessments[[i]] = list(
      assessment=assessment,
      transformed=transformed_x
    )
  }
  
  obj = list(
    original=dataset, 
    assessments=assessments
  )
  # obj is now a list of class `list` and `normal`
  class(obj) = append(class(obj),"normal") 
  return(obj)
}



summary.normal <- function (object, ...) {
  assessments = object$assessments
  dataset = object$original
  for (i in 1:length(assessments)){
    # Initial test of normality
    assessment = assessments[[i]]
    name = names(dataset)[[i]]
    significance = assessment$assessment$significance
    if(assessment$assessment$isnormal){
      print(sprintf('Variable %s is normal to within %.2f%% significance!', name, significance))
    }else{
      print(sprintf('Variable %s is not normal within %.2f%% significance test!', name, significance))
    }
  }
}

print.normal <- function (object, ...) {
  print(paste(
      sprintf('Dataset has %i variables: ',length(variable.names(object$dataset))),
      paste(variable.names(iris),collapse='; ')
  ))
}


plot.normal <- function (object, ...) {
  assessments = object$assessments
  first = assessments[[1]]
  
  if (all(is.na(first$transformed))){
    assessment = first$assessment
  }else{
    lambda = first$transformed$lambda
    print(sprintf('Transformed with lambda=%.2f',lambda))
    assessment = first$transformed$assessment
  }
  
  # build data frame
  df <- data.frame(
    q=assessment$quartile_samples,
    sorted=assessment$sorted_sample,
    predicted=predict(assessment$model),
    residuals=residuals(assessment$model)
  )
  
  model_a=assessment$model$coefficients[[1]]
  model_b=assessment$model$coefficients[[2]]
  
  # QQ Plot showing residuals 
  #
  # Plot thanks to this fantastic tutorial:
  # 
  # https://drsimonj.svbtle.com/visualising-residuals
  #
  ggplot(df, aes(q, sorted)) +
  geom_segment(aes(xend = q, yend = predicted), alpha = 0.4) +
  geom_point(aes(color = residuals)) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red") +
  guides(color = FALSE) + theme_bw() + 
  geom_abline(intercept = model_a, slope = model_b, linetype="dashed", size=0.15)
}

# load iris dataset
data('iris')

# just use one species
setosa = subset(iris, Species == "setosa")[1:5]

setosa_n  = normal(dataset=setosa)

plot(setosa_n)

# summary(iris_red)
# print(iris_red)
# plot(iris_red)