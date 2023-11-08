#' @title Small Area Estimation using Hierarchical Bayesian under Binomial Distribution
#' @description This function is implemented to variable of interest \eqn{(y)} that assumed to be a Binomial Distribution using Logit normal model. The data is an accumulation from the Bernoulli process which there are exactly two mutually exclusive outcomes of a case.
#' @author Azka Ubaidillah [aut], Ika Yuni Wulansari [aut], Zaza Yuda Perwira [aut, cre], Jayanti Wulansari [aut, cre], Fauzan Rais Arfizain [aut,cre]
#' @param formula Formula that describe the fitted model
#' @param n.samp number of sample
#' @param iter.update Number of updates with default \code{3}
#' @param iter.mcmc Number of total iterations per chain with default \code{10000}
#' @param coef a vector contains prior initial value of Coefficient of Regression Model for fixed effect with default vector of \code{0} with the length of the number of regression coefficients
#' @param var.coef a vector contains prior initial value of variance of Coefficient of Regression Model with default vector of \code{1} with the length of the number of regression coefficients
#' @param thin Thinning rate, must be a positive integer with default \code{2}
#' @param burn.in Number of iterations to discard at the beginning with default \code{2000}
#' @param tau.u  Prior initial value of inverse of Variance of area random effect with default \code{1}
#' @param data The data frame
#'
#' @return This function returns a list of the following objects:
#'    \item{Est}{A vector with the values of Small Area mean Estimates using Hierarchical bayesian method, if there is non sampled area, then the result of estimation of mu in nonsampled areas are the probabilites}
#'    \item{refVar}{Estimated random effect variances}
#'    \item{coefficient}{A dataframe with the estimated model coefficient}
#'    \item{plot}{Trace, Dencity, Autocorrelation Function Plot of MCMC samples}
#'
#' @export Binomial
#'
#' @examples
#'
#' \donttest{
#' #Data Generation
#' set.seed(123)
#' m=30
#' x1=runif(m,0,1)
#' x2=runif(m,0,1)
#' b0=b1=b2=0.5
#' u=rnorm(m,0,1)
#' n.samp1=round(runif(m,10,30))
#' mu= exp(b0 + b1*x1+b2*x2+u)/(1+exp(b0 + b1*x1+b2*x2+u))
#' y=rbinom(m,n.samp1,mu)
#' vardir=n.samp1*mu*(1-mu)
#' dataBinomial=as.data.frame(cbind(y,x1,x2,n.samp=n.samp1,vardir))
#' dataBinomialNs = dataBinomial
#' dataBinomialNs$y[c(3,14,22,29,30)] <- NA
#' dataBinomialNs$vardir[c(3,14,22,29,30)] <- NA
#' dataBinomialNs$n.samp[c(3,14,22,29,30)] <- NA
#'
#'
#' ##Compute Fitted Model
#' ##y~x1+x2
#'
#'
#' ## For data without any nonsampled area
#' formula = y~x1+x2
#' n.s = "n.samp"
#' vc = c(1,1,1)
#' c = c(0,0,0)
#' dat = dataBinomial
#'
#'
#' ## Using parameter coef and var.coef
#' saeHBBinomial<-Binomial(formula,n.samp=n.s,iter.update=10,coef=c,var.coef=vc,data =dat)
#'
#' saeHBBinomial$Est                                 #Small Area mean Estimates
#' saeHBBinomial$refVar                              #Random effect variance
#' saeHBBinomial$coefficient                         #coefficient
#' #Load Library 'coda' to execute the plot
#' #autocorr.plot(saeHBBinomial$plot[[3]]) is used to generate ACF Plot
#' #plot(saeHBBinomial$plot[[3]]) is used to generate Density and trace plot
#'
#' ## Do not using parameter coef and var.coef
#' saeHBBinomial <- Binomial(formula,n.samp ="n.samp",data=dataBinomial)
#'
#'
#'
#' ## For data with nonsampled area use dataBinomialNs
#'
#' }
Binomial <- function(formula, n.samp, iter.update=3, iter.mcmc=10000,coef, var.coef,thin = 2, burn.in =2000, tau.u = 1, data){



  result <- list(Est = NA, refVar = NA, coefficient = NA,
                 plot=NA)


  formuladata <- model.frame(formula,data,na.action=NULL)
  if (any(is.na(formuladata[,-1])))
    stop("Auxiliary Variables contains NA values.")
  auxVar <- as.matrix(formuladata[,-1])
  nvar <- ncol(auxVar) + 1
  formuladata <- data.frame(formuladata, n.samp = data[,n.samp])

  if (!missing(var.coef)){

    if( length(var.coef) != nvar ){
      stop("length of vector var.coef does not match the number of regression coefficients, the length must be ",nvar)
    }

    tau.b.value = 1/var.coef
  } else {
    tau.b.value = 1/rep(1,nvar)
  }

  if (!missing(coef)){
    if( length(coef) != nvar ){
      stop("length of vector coef does not match the number of regression coefficients, the length must be ",nvar)
    }
    mu.b.value = coef
  } else {
    mu.b.value = rep(0,nvar)
  }

  if (iter.update < 3){
    stop("the number of iteration updates at least 3 times")
  }

  #Fungsi Tersampel
  if (!any(is.na(formuladata[,1]))){

    formuladata <- as.matrix(na.omit(formuladata))
    x <- model.matrix(formula,data = as.data.frame(formuladata))
    n <- nrow(formuladata)
    mu.b = mu.b.value
    tau.b = tau.b.value
    tau.ua = tau.ub = a.var= 1
    for (i in 1:iter.update){
      dat <- list("n"= n, "n.samp" = formuladata[,4], "nvar"= nvar, "y" = formuladata[,1], "x"=as.matrix(x[,-1]),
                  "mu.b"=mu.b, "tau.b"=tau.b, "tau.ua"=tau.ua, "tau.ub"=tau.ub)
      inits <- list(u = rep(0, n),b = mu.b, tau.u =tau.u)
      cat("model {
					for (i in 1:n) {

					    y[i]~dbin(p[i], n.samp[i])
              logit(p[i]) <-  b[1] + sum(b[2:nvar]*x[i,]) + u[i]
              u[i]~dnorm(0, tau.u)
              mu[i]<- n.samp[i]*p[i]


					}
          # prior
          for (k in 1:nvar){
					    b[k] ~ dnorm(mu.b[k],tau.b[k])
          }
					tau.u~dgamma(tau.ua, tau.ub)
          a.var <- 1/tau.u

		  	}", file="Binomial.txt")
      jags.m <- jags.model( file = "Binomial.txt", data=dat, inits=inits, n.chains=1, n.adapt=500 )
      file.remove("Binomial.txt")
      params <- c("mu","a.var","b", "tau.u")
      samps <- coda.samples( jags.m, params, n.iter=iter.mcmc, thin=thin)
      samps1 <- window(samps, start=burn.in+1, end=iter.mcmc)
      result_samps=summary(samps1)
      a.var=result_samps$statistics[1]
      beta=result_samps$statistics[2:(nvar+1),1:2]
      for (i in 1:nvar){
        mu.b[i]  = beta[i,1]
        tau.b[i] = 1/(beta[i,2]^2)
      }

      tau.ua = result_samps$statistics[(2+nvar+n),1]^2/result_samps$statistics[(2+nvar+n),2]^2
      tau.ub = result_samps$statistics[(2+nvar+n),1]/result_samps$statistics[(2+nvar+n),2]^2
    }
    result_samps=summary(samps1)
    b.varnames <- list()
    for (i in 1:(nvar)) {
      idx.b.varnames <- as.character(i-1)
      b.varnames[i] <-str_replace_all(paste("b[",idx.b.varnames,"]"),pattern=" ", replacement="")
    }

    result_mcmc <- samps1[,c(2:(nvar+1))]
    colnames(result_mcmc[[1]]) <- b.varnames

    a.var=result_samps$statistics[1]

    beta=result_samps$statistics[2:(nvar+1),1:2]
    rownames(beta) <- b.varnames

    mu=result_samps$statistics[(nvar+2):(1+nvar+n),1:2]

    Estimation=data.frame(mu)

    Quantiles <- as.data.frame(result_samps$quantiles[1:(2+nvar+n),])
    q_mu <- Quantiles[(nvar+2):(nvar+1+n),]
    q_beta <- (Quantiles[2:(nvar+1),])
    rownames(q_beta) <- b.varnames
    beta <- cbind(beta,q_beta)
    Estimation <- data.frame(Estimation,q_mu)
    colnames(Estimation) <- c("MEAN","SD","2.5%","25%","50%","75%","97.5%")
  } else {

    formuladata <- as.data.frame(formuladata)

    x <- as.matrix(formuladata[,2:nvar])
    n <- nrow(formuladata)

    mu.b =mu.b.value
    tau.b = tau.b.value
    tau.ua = tau.ub = a.var =1

    formuladata$idx <- rep(1:n)
    data_sampled <- na.omit(formuladata)
    data_nonsampled <- formuladata[-data_sampled$idx,]

    r=data_nonsampled$idx
    n1=nrow(data_sampled)
    n2=nrow(data_nonsampled)
    for (i in 1:iter.update){
      dat <- list("n1"= n1, "n2"=n2,"nvar"=nvar,
                  "n.samp_sampled" = data_sampled[,4], "n.samp_nonsampled" = data_nonsampled[,4],
                  "y_sampled" = data_sampled[,1],
                  "x_sampled"=as.matrix(data_sampled[,2:nvar]),
                  "x_nonsampled"=as.matrix(data_nonsampled[,2:nvar]),
                  "mu.b"=mu.b,"tau.b"=tau.b, "tau.ua"=tau.ua,"tau.ub"=tau.ub)
      inits <- list(u = rep(0,n1), v = rep(0,n2), b = mu.b, tau.u = tau.u)
      cat("model {
          for (i in 1:n1) {

              y_sampled[i]~dbin(p[i], n.samp_sampled[i])
              logit(p[i]) <- b[1] + sum(b[2:nvar]*x_sampled[i,])
              u[i]~dnorm(0, tau.u)
              mu[i]<- n.samp_sampled[i]*p[i]

          }

          for(j in 1:n2) {
              v[j]~dnorm(0, tau.u)
              y_nonsampled[j]~dbin(p.nonsampled[j], 1)
              logit(p.nonsampled[j]) <-  mu.b[1] + sum(mu.b[2:nvar]*x_nonsampled[j,])+ v[j]
              mu.nonsampled[j]<- p.nonsampled[j]
          }

					# prior
          for (k in 1:nvar){
					    b[k] ~ dnorm(mu.b[k],tau.b[k])
          }
          tau.u~dgamma(tau.ua, tau.ub)
          a.var <- 1/tau.u

			}", file="Binomial.txt")

      jags.m <- jags.model( file = "Binomial.txt", data=dat, inits=inits, n.chains=1, n.adapt=500 )
      file.remove("Binomial.txt")
      params <- c("mu","mu.nonsampled","a.var","b", "tau.u")
      samps <- coda.samples( jags.m, params, n.iter=iter.mcmc, thin=thin)
      samps1 <- window(samps, start=burn.in+1, end=iter.mcmc)
      result_samps=summary(samps1)
      a.var=result_samps$statistics[1]
      beta=result_samps$statistics[2:(nvar+1),1:2]
      for (i in 1:nvar){
        mu.b[i]  = beta[i,1]
        tau.b[i] = 1/(beta[i,2]^2)
      }

      tau.ua = result_samps$statistics[2+nvar+n,1]^2/result_samps$statistics[2+nvar+n,2]^2
      tau.ub = result_samps$statistics[2+nvar+n,1]/result_samps$statistics[2+nvar+n,2]^2
    }
    result_samps=summary(samps1)
    b.varnames <- list()
    for (i in 1:(nvar)) {
      idx.b.varnames <- as.character(i-1)
      b.varnames[i] <-str_replace_all(paste("b[",idx.b.varnames,"]"),pattern=" ", replacement="")
    }
    result_mcmc <- samps1[,c(2:(nvar+1))]
    colnames(result_mcmc[[1]]) <- b.varnames

    a.var=result_samps$statistics[1]

    beta=result_samps$statistics[2:(nvar+1),1:2]
    rownames(beta) <- b.varnames

    mu=round(result_samps$statistics[(nvar+2):(1+nvar+n1),1:2])
    mu.nonsampled =result_samps$statistics[(2+nvar+n1):(1+nvar+n),1:2]

    Estimation=matrix(rep(0,n),n,2)
    Estimation[r,]=mu.nonsampled
    Estimation[-r,]=mu
    Estimation = as.data.frame(Estimation)

    Quantiles <- as.data.frame(result_samps$quantiles[1:(2+nvar+n),])
    q_beta <- (Quantiles[2:(nvar+1),])
    q_mu <- (Quantiles[(nvar+2):(nvar+1+n1),])
    q_mu.nonsampled <- (Quantiles[(2+nvar+n1):(1+nvar+n),])
    q_Estimation <- matrix(0,n,5)
    for (i in 1:5){
      q_Estimation[r,i] <- q_mu.nonsampled[,i]
      q_Estimation[-r,i] <- q_mu[,i]
    }

    rownames(q_beta) <- b.varnames
    beta <- cbind(beta,q_beta)
    Estimation <- data.frame(Estimation,q_Estimation)
    colnames(Estimation) <- c("MEAN","SD","2.5%","25%","50%","75%","97.5%")
  }

  result$Est         = Estimation
  result$refVar      = a.var
  result$coefficient = beta
  result$plot        = list(graphics.off() ,par(mar=c(2,2,2,2)),autocorr.plot(result_mcmc,col="brown2",lwd=2),plot(result_mcmc,col="brown2",lwd=2))
  return(result)

}
