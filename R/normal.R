
critical_points = function(){
  fpath <- system.file("exdata", "table_4_2_qq_critical.csv", package="project2normal")
  # Table 4.2 from the book
  return(read.csv(fpath))
} 



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
  qq_critical_points = critical_points()
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
#' @export
boxcox = function(x,lambda){
  library(pracma)
  
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
  # x must always be greater than 0
  x = repmat(matrix(x - min(x) + 1,nr=m,nc=1),1,n)
  lambda = repmat(matrix(lambda,nr=1,nc=n), m, 1)
  
  stopifnot(all(dim(x) == dim(lambda)))
  # use box cox transformation; mask to select which operation to perform
  x_l = x
  zmask = all.equal(lambda,0,tolerance=0.1) == TRUE
  nzmask = !zmask
  x_l[zmask] = log(x[zmask])
  x_l[nzmask] = (x[nzmask]**lambda[nzmask]-1)/lambda[nzmask]
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
#' @export
box_cox_likelihood = function(x,lambda){
  # box-cox likelihood function, target for maximization
  n = length(x)
  # x must always be greater than 0
  x = x - min(x) + 1
  xs = sum(log(x))
  x_to_the_lambda = boxcox(x,lambda)
  l = -n/2 * log(1/n * colSums((x_to_the_lambda-colMeans(x_to_the_lambda))^2)) + (lambda-1)*xs
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
#' assessment = normal_assessment(iris[1:25,1], significance=5)
#' 
#' @export
normal_assessment = function(x, significance=5){
  stopifnot(significance>=1,significance<=10)
  # Only allow one dimensional datasets
  dims = dim(x)
  stopifnot(length(dims)<=2)
  if (length(dims) == 2){
    stopifnot(min(dims)==1)
  }
  name = colnames(x)[1]
  # reduce to 1d
  x = matrix(drop(as.matrix(x)),nc=1)
  colnames(x)[1] <- name
  
  # number of samples in variable
  n = length(drop(x))
  
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

#' Summarize normal_assessment object; called implicitely.  
#' 
#' @param object normal_assessment object
#' 
#' @examples
#' 
#' # Load iris dataset 
#' data('iris')
#'
#' summary(normal_assessment(iris[1:25,1], significance=5))
#' 
#' @export
summary.normal_assessment <- function (object, ...) {
  library(crayon)
  
  nms = colnames(object$x)
  if (is.null(nms)){
    name = 'V1'
  }else{
    name = nms[[1]]
  }
  dsum =  paste(capture.output(summary(drop(object$x))),collapse='\n')
  model = sprintf('y = %fx + %f',object$model$coefficients[[2]], object$model$coefficients[[1]])
  str = '\n'
  str = paste(str, sprintf('Variable Name: %s, Number of samples (n): %i\n\n', blue(name), length(object$x)), sep='')
  str = paste(str, 'Basic Data stats:\n', dsum, '\n', sep='')
  str = paste(str, 'Normality Assessment:\n', sep='')
  # Unicode characters apparently break the formatter, so I had to increase the alignment length by 1 hear
  str = paste(str, sprintf('%22s : %.2f%% \n','Significance (\U03B1)',object$significance), sep='')
  str = paste(str, sprintf('%21s : %.4f \n','QQ Critical Point',object$critical_point), sep='')
  str = paste(str, sprintf('%21s : [%s] \n','Normal Model',model), sep='')
  str = paste(str, sprintf('%21s : %.4f \n','Measured Correlation',object$correlation), sep='')
  if(object$isnormal){
    assessment = green("IS \U2714")
  }else{
    assessment = red("IS NOT \U2717")
  }
  str = paste(str, sprintf('\nDataset %s normal for significance level %.2f%%',assessment, object$significance), sep='')
  cat(str)
}

#' Summarize normal_assessment object; called implicitely.  
#' 
#' @param object normal_assessment object
#' 
#' @examples
#' 
#' # Load iris dataset 
#' data('iris')
#'
#' summary(normal_assessment(iris[1:25,1], significance=5))
#' 
#' @export
print.normal_assessment <- function (object, ...) {
  summary(object)
}

#' Create QQ plot from a normality_assessment object; called implicitely.  
#' 
#' @param object normal_assessment object
#' 
#' @examples
#' 
#' # Load iris dataset 
#' data('iris')
#'
#' plot(normal_assessment(iris[1:25,1], significance=5))
#' 
#' @export
plot.normal_assessment <- function (object, ...) {
  library(ggplot2)
  
  # build data frame
  df <- data.frame(
    q=object$quartile_samples,
    sorted=object$sorted_sample,
    predicted=predict(object$model),
    residuals=residuals(object$model)
  )
  
  name = colnames(object$x)[1]
  if (object$isnormal){
    msg = 'Normal'
  }else{
    msg = 'Not normal'
  }
  
  model_a=object$model$coefficients[[1]]
  model_b=object$model$coefficients[[2]]

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
  geom_abline(intercept = model_a, slope = model_b, linetype="dashed", size=0.15) + 
  ggtitle(sprintf('%s, %s; R=%.2f',name,msg,object$correlation)) + 
    xlab("Quantiles") + ylab("Observations")
}


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
#' data('iris')
#' irisnorm = normal(iris[1:25,])
#' 
#' @export
normal = function(dataset, significance=5, boxcox_lambda_search=seq(-5,5,0.01))
{
  library(crayon)
  
  stopifnot(!missing(dataset))
  # if this fails you're in trouble
  dataset = data.frame(dataset)
  stopifnot(length(dim(dataset)) < 3)
  stopifnot(length(dim(dataset)) > 0)
  
  str = 'Creating normal object for dataset!\n'
  
  # extract portion of dataset which is numeric
  dset_numeric = dataset[,sapply(dataset, function(x) is.numeric(x)),drop=FALSE]
  excluded = setdiff(names(dataset), names(dset_numeric))
  if (length(excluded)>0){
    str = paste(str, sprintf(
      red('Dataset has non-numeric variables [%s]!\n'),
      paste(excluded, collapse = ', ')), 
      ' * Please be sure each column of the dataset includes samples from ',
      red('only'),
      ' a single population.\n\n',
      sep=''
    )
  }
  n = dim(dset_numeric)[1]
  p = dim(dset_numeric)[2]
  str = paste(str, sprintf(
    'Checking normality of N=%i samples of the variables [%s]...\n',
    n,
    paste(names(dset_numeric), collapse = ', ')),
    sep=''
  )
  assessments <- vector("list", p)
  for (i in 1:p){
    name = colnames(dset_numeric)[i]
    x = as.matrix(dset_numeric[,i,drop=FALSE])
    # Initial test of normality
    assessment = normal_assessment(x, significance=significance)

    if(assessment$isnormal){
      str = paste(str, sprintf(
        '- Variable %s is normal to within %.2f%% significance!\n',
        name,
        significance),
        sep=''
      )
      transformed = NA
    }else{
      str = paste(str, sprintf(
        '- Variable %s failed the %.2f%% significance test!\n',
        name,
        significance),
        sep=''
      )
      str = paste(str, ' * Attempting to find a boxcox transformation to normalize the dataset...\n', sep='')
      # Estimate lambda for boxcox transformation to attempt a transform to normality
      lambda = boxcox_lambda_search[which.max(
        box_cox_likelihood(x, boxcox_lambda_search)
      )]
      
      str = paste(str, sprintf(' * \U03BB=%.2f maximizes the boxcox likelihood function!\n',lambda), sep='')
      # could not get R's optimization tools to work, so just call this "close enough".
      # It's an estimate anyway.
      d = boxcox(x, lambda)
      colnames(d) = colnames(x)
      #stop()
      transformed = list(
        lambda=lambda,
        assessment=normal_assessment(d, significance=significance)
      )
      if (transformed$assessment$isnormal){
        result = green('IS')
        using = 'Keeping transformation.'
      }else{
        result = red('IS NOT')
        transformed = NA
        using = 'Using original data'
      }
      str = paste(str, sprintf(' * Transformed data %s sufficiently normal! %s\n', result, using), sep='')
    }

    # store for later
    assessments[[i]] = list(
      assessment=assessment,
      transformed=transformed
    )
  }
  obj = list(
    dataset=dataset,
    assessments=assessments
  )
  cat(str)
  # obj is now a list of class `list` and `normal`
  class(obj) = append(class(obj),"normal")
  return(obj)
}

#' Run normal on a dataset and return a "normified" dataset, transformed when necessary 
#' 
#' @param data 2D Dataframe to normify.
#' 
#' @examples
#' 
#' # Load iris dataset 
#' data('iris')
#'
#' normifed_data = normify(iris[1:25,], significance=5)
#' 
#' @export
normify = function(data){
  UseMethod('normify', data)
}
#' @export
normify.default = function(data){
  return(normify(normal(data)))
}

#' Return "normified" dataset from a preallocated normal object 
#' 
#' @param data Dataset to normify
#' 
#' @export
normify.normal <- function(object){
  assessments = object$assessments
  datalist = list()
  for (i in 1:length(assessments)){
    assessment = assessments[[i]]$assessment
    transformed = assessments[[i]]$transformed
    lambda = ''
    if (!all(is.na(transformed))){
      # Data was transformed
      # Use transformed assessment
      assessment = transformed$assessment
      lambda = sprintf('.L.%.2f',transformed$lambda) 
    }
    datalist[[i]] = assessment$x
    colnames(datalist[[i]])[1] = paste(colnames(datalist[[i]])[1], lambda, sep='')
  }
  return(as.data.frame(datalist))
}

#' Summarize content of a normal object; called implicitely.  
#' 
#' @param object normal object
#' 
#' @examples
#' 
#' # Load iris dataset 
#' data('iris')
#'
#' summary(normal(iris[1:25,], significance=5))
#' 
#' @export
summary.normal <- function (object, ...) {
  library(crayon)
  
  assessments = object$assessments
  for (i in 1:length(assessments)){
    cat(yellow(sprintf('-----------------  Column %i  ------------------',i)))
    # Initial test of normality
    assessment = assessments[[i]]$assessment
    transformed = assessments[[i]]$transformed
    if (!all(is.na(transformed))){
      # Data was transformed
      cat(red(sprintf('\n\n !! Variable was transformed with \U03BB=%.2f !!\n',transformed$lambda)))
      # Use transformed assessment
      assessment = transformed$assessment
    }
    summary(assessment)
    cat('\n\n')
  }
}

#' Print some info about a normal object; called implicitely.  
#' 
#' @param object normal object
#' 
#' @examples
#' 
#' # Load iris dataset 
#' data('iris')
#'
#' print(normal(iris[1:25,], significance=5))
#' 
#' @export
print.normal <- function (object, ...) {
  library(crayon)
  
  assessments = object$assessments
  dataset = object$dataset
  cat(sprintf('\n%-20s| %10s | %-10s\n','Variable name', 'Normal?','Transformed'))
  cat('---------------------------------------------------\n')
  for (i in 1:length(assessments)){
    # Initial test of normality
    assessment = assessments[[i]]$assessment
    transformed = assessments[[i]]$transformed
    str = 'No'
    if (!all(is.na(transformed))){
      # Data was transformed
      str = red(sprintf('\U03BB=%.2f',transformed$lambda))
      # Use transformed assessment
      assessment = transformed$assessment
    }
    if (assessment$isnormal){
      isnormal = 'Normal'
    }else{
      isnormal = 'Not Normal'
    }
    cat(sprintf('%-20s| %10s | %-10s\n', colnames(dataset)[i] ,isnormal, str))
  }
}

#' Generate QQ plots for each variable in the dataset.  
#' 
#' @param object normal object
#' 
#' @examples
#' 
#' # Load iris dataset 
#' data('iris')
#'
#' plot(normal(iris[1:25,], significance=5))
#' 
#' @export
plot.normal <- function (object, ...) {
  library(gridExtra)
  
  assessments = object$assessments
  ps = list()
  for (i in 1:length(assessments)){
    assessment = assessments[[i]]
    if (all(is.na(assessment$transformed))){
      assessment = assessment$assessment
    }else{
      lambda = assessment$transformed$lambda
      assessment = assessment$transformed$assessment
    }
    ps[[i]] = plot(assessment)
  }
  grid.arrange(grobs=ps)
}