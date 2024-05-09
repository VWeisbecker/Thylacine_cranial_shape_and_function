# Function script

# Benjamini-Hochberg FDR procedure
# function(x, a) where x = p-values, a = False Discovery Rate (FDR)

benj.hoch = function(x, a){
  i = seq(along=x)
  m = length(x)
  k = i/m*a
  Q = max( which( sort(x) < k))
  cutoff <- sort(x)[Q]
  sig1 = which(x <= cutoff)
  sig2 = x[x <= cutoff]
  output = list(cutoff, sig1, sig2)
  names(output) = c("Significance cutoff", "H0 reject", "Significant p-values")
  return(output)
}


# From https://github.com/JeroenSmaers/evomap
#' GLS confidence intervals
#'
#' Computes confidence intervals for a GLS regression. When the variance covariance matrix is based on a phylogeny, this computes the confidence intervals for a phylogenetic regression (PGLS).
#' @param Y Dependent variable.
#' @param X Independent variable.
#' @param Sigma Matrix describing the covariance among data points. For the phylogenetic regression, this is based on a phylogeny (vcv(tree)).
#' @return GLS results.
#' @return Confidence intervals of the estimates ('CI').
#' @return Confidence intervals of unobserved values of Y ('CI.plot').
#' @references Smaers & Rohlf (2016) Testing species' deviations from allometric preductions using the phylogenetic regression. Evolution. 70 (5): 1145-1149.
#' @examples see https://smaerslab.com/software/ 
#' @export
gls.ci<-function(Y,X,Sigma){
  n<-length(X)
  #Standardize sigma to reduce rounding errors
  tr<-sum(diag(Sigma))
  Sigma<-n*Sigma/tr
  #Input
  invSigma<-solve(Sigma)
  #pgls mean and variance
  X1<-rep(1,n)
  q<-2          # adjust if multivariate!!
  C1<-solve(t(X1)%*%invSigma%*%X1)
  Y_PGLSmean<-c(C1%*%t(X1)%*%invSigma%*%Y)
  Y_PGLSdeviations = Y - Y_PGLSmean
  Y_PGLSvariance = (t(Y_PGLSdeviations)%*%invSigma%*%Y_PGLSdeviations)/(n-1)
  SE_Y_mean = sqrt(Y_PGLSvariance/n)
  #Regression
  XX<-cbind(rep(1,n),X)
  C<-solve(t(XX)%*%invSigma%*%XX)
  w<-C%*%t(XX)%*%invSigma
  B<-w%*%Y
  Yhat<-XX%*%B
  Yresid=Y-Yhat
  Y_MSEresid<-c((t(Yresid)%*%invSigma%*%Yresid)/(n-q))
  a<-B[1]
  b<-B[2]
  #SE of statistics
  SEa<-sqrt(diag(C)*Y_MSEresid)[1]
  SEb<-sqrt(diag(C)*Y_MSEresid)[2]
  #Save model
  intercept<-cbind(a,SEa)
  slope<-cbind(b,SEb)
  model<-rbind(intercept,slope)
  colnames(model)<-c("Estimate","Std.Error")
  rownames(model)<-c("intercept","slope")
  #confint
  SEYhat<-sqrt(diag(XX%*%C%*%t(XX))%*%((t(Yresid)%*%invSigma%*%Yresid)/(n-q)))
  CI<-cbind(X,Yhat,SEYhat)
  Lower2.5<-Yhat-qt(0.975,n)*SEYhat
  Lower5<-Yhat-qt(0.95,n)*SEYhat
  Upper5<-Yhat+qt(0.95,n)*SEYhat
  Upper2.5<-Yhat+qt(0.975,n)*SEYhat
  CI<-cbind(CI,Lower2.5)
  CI<-cbind(CI,Lower5)
  CI<-cbind(CI,Upper5)
  CI<-cbind(CI,Upper2.5)
  CI<-CI[order(CI[,1]),]
  colnames(CI)<-c("X","Yhat","SEYhat","Lower2.5","Lower5","Upper5","Upper2.5")
  CI<-as.data.frame(CI)
  #confint extrapolated
  Xi<-seq(c(min(X)-abs(max(X))),to=c(abs(max(X))*5),length.out=100)
  Z<-c(Xi)
  ZZ<-cbind(rep(1,length(Z)),Z)
  SEYhat.Xi<-sqrt(diag(ZZ%*%C%*%t(ZZ))%*%((t(Yresid)%*%invSigma%*%Yresid)/(n-q)))
  Yhat.Xi<-a+b*Xi
  Lower2.5.Yhat.Xi<-Yhat.Xi-qt(0.975,n)*SEYhat.Xi
  Lower5.Yhat.Xi<-Yhat.Xi-qt(0.95,n)*SEYhat.Xi
  Upper5.Yhat.Xi<-Yhat.Xi+qt(0.95,n)*SEYhat.Xi
  Upper2.5.Yhat.Xi<-Yhat.Xi+qt(0.975,n)*SEYhat.Xi
  CI.plot<-cbind(Xi,Yhat.Xi,SEYhat.Xi,Lower2.5.Yhat.Xi,Lower5.Yhat.Xi,Upper5.Yhat.Xi,Upper2.5.Yhat.Xi)
  colnames(CI.plot)<-c("X","Yhat","SEYhat","Lower2.5","Lower5","Upper5","Upper2.5")
  CI.plot<-as.data.frame(CI.plot)
  results<-list(model,CI,CI.plot)
  names(results)<-c("model","CI","CI.plot")
  return(results)
}


# From https://github.com/JeroenSmaers/evomap
#' GLS prediction intervals
#'
#' Computes prediction intervals for a GLS regression. When the variance covariance matrix is based on a phylogeny, this computes the prediction intervals for a phylogenetic regression (PGLS). 
#' @param Y Dependent variable. 
#' @param X Independent variable. 
#' @param Sigma Matrix describing the covariance among data points. For the phylogenetic regression, this is based on a phylogeny (vcv(tree)). 
#' @param k The sample size of the estimated future points.
#' @return GLS results.
#' @return Prediction intervals of the estimates ('PI').
#' @return Prediction intervals of unobserved values of Y ('PI.plot').
#' @references Smaers & Rohlf (2016) Testing species' deviations from allometric preductions using the phylogenetic regression. Evolution. 70 (5): 1145-1149.
#' @examples see https://smaerslab.com/software/
#' @export

gls.pi<-function(Y,X,Sigma,k){
  n<-length(X)
  #Standardize sigma to reduce rounding errors
  tr<-sum(diag(Sigma))
  Sigma<-n*Sigma/tr
  #Input
  invSigma<-solve(Sigma)
  #pgls mean and variance
  X1<-rep(1,n)
  q<-2          # correct if multivariate!!
  C1<-solve(t(X1)%*%invSigma%*%X1)
  Y_PGLSmean<-c(C1%*%t(X1)%*%invSigma%*%Y)
  Y_PGLSdeviations = Y - Y_PGLSmean
  Y_PGLSvariance = (t(Y_PGLSdeviations)%*%invSigma%*%Y_PGLSdeviations)/(n-1)
  SE_Y_mean = sqrt(Y_PGLSvariance/n)
  #Regression
  XX<-cbind(rep(1,n),X)
  C<-solve(t(XX)%*%invSigma%*%XX)
  w<-C%*%t(XX)%*%invSigma
  B<-w%*%Y
  Yhat<-XX%*%B
  Yresid=Y-Yhat
  Y_MSEresid<-c((t(Yresid)%*%invSigma%*%Yresid)/(n-q))
  a<-B[1]
  b<-B[2]
  #SE of statistics
  SEa<-sqrt(diag(C)*Y_MSEresid)[1]
  SEb<-sqrt(diag(C)*Y_MSEresid)[2]
  #Save model
  intercept<-cbind(a,SEa)
  slope<-cbind(b,SEb)
  model<-rbind(intercept,slope)
  colnames(model)<-c("Estimate","Std.Error")
  rownames(model)<-c("intercept","slope")
  #predint
  SEYhat<-sqrt((diag(XX%*%C%*%t(XX))+c(1/k))%*%((t(Yresid)%*%invSigma%*%Yresid)/(n-q)))
  PI<-cbind(X,Yhat,SEYhat)
  Lower2.5<-Yhat-qt(0.975,n)*SEYhat
  Lower5<-Yhat-qt(0.95,n)*SEYhat
  Upper5<-Yhat+qt(0.95,n)*SEYhat
  Upper2.5<-Yhat+qt(0.975,n)*SEYhat
  PI<-cbind(PI,Lower2.5)
  PI<-cbind(PI,Lower5)
  PI<-cbind(PI,Upper5)
  PI<-cbind(PI,Upper2.5)
  PI<-PI[order(PI[,1]),]
  colnames(PI)<-c("X","Yhat","SEYhat","Lower2.5","Lower5","Upper5","Upper2.5")
  PI<-as.data.frame(PI)
  #confint extrapolated
  Xi<-seq(c(min(X)-abs(max(X))),to=c(abs(max(X))*5),length.out=100)
  Z<-c(Xi)
  ZZ<-cbind(rep(1,length(Z)),Z)
  SEYhat.Xi<-sqrt((diag(ZZ%*%C%*%t(ZZ))+c(1/k))%*%((t(Yresid)%*%invSigma%*%Yresid)/(n-q)))
  Yhat.Xi<-a+b*Xi
  Lower2.5.Yhat.Xi<-Yhat.Xi-qt(0.975,n)*SEYhat.Xi
  Lower5.Yhat.Xi<-Yhat.Xi-qt(0.95,n)*SEYhat.Xi
  Upper5.Yhat.Xi<-Yhat.Xi+qt(0.95,n)*SEYhat.Xi
  Upper2.5.Yhat.Xi<-Yhat.Xi+qt(0.975,n)*SEYhat.Xi
  PI.plot<-cbind(Xi,Yhat.Xi,SEYhat.Xi,Lower2.5.Yhat.Xi,Lower5.Yhat.Xi,Upper5.Yhat.Xi,Upper2.5.Yhat.Xi)
  colnames(PI.plot)<-c("X","Yhat","SEYhat","Lower2.5","Lower5","Upper5","Upper2.5")
  PI.plot<-as.data.frame(PI.plot)
  results<-list(model,PI,PI.plot)
  names(results)<-c("model","PI","PI.plot")
  return(results)
}