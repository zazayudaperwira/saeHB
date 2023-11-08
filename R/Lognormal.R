#' @title Small Area Estimation using Hierarchical Bayesian under Lognormal Distribution
#' @description This function is implemented to variable of interest \eqn{(y)} that assumed to be a Lognormal Distribution. The range of data is \eqn{(y > 0}
#' @param formula Formula that describe the fitted model
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
#'    \item{Est}{A vector with the values of Small Area mean Estimates using Hierarchical bayesian method }
#'    \item{refVar}{Estimated random effect variances}
#'    \item{coefficient}{A dataframe with the estimated model coefficient}
#'    \item{plot}{Trace, Dencity, Autocorrelation Function Plot of MCMC samples}
#'
#' @export Lognormal
#'
#' @examples
#' \donttest{
#' ##Data Generation
#' set.seed(123)
#' m=30
#' x1=runif(m,0,1)
#' x2=runif(m,1,2)
#' x3=runif(m,2,3)
#' b0=b1=b2=b3=0.5
#' u=rnorm(m,0,1)
#' mu=b0 + b1*x1+b2*x2+b3*x3+u
#' sig=1
#' y=rlnorm(m,mu,sig)
#' E=exp(mu+1/2*sig^2)
#' vardir=exp(2*mu+sig^2)*(exp(sig^2)-1)
#' dataLognormal=as.data.frame(cbind(y,x1,x2,x3,vardir))
#' dataLognormalNs=dataLognormal
#' dataLognormalNs$y[c(3,14,22,29,30)] <- NA
#' dataLognormalNs$vardir[c(3,14,22,29,30)] <- NA
#'
#'
#' ##Compute Fitted Model
#' ##y ~ x1 +x2 +x3
#'
#'
#' ## For data without any nonsampled area
#'
#' formula = y ~ x1 +x2 +x3
#' v = c(1,1,1,1)
#' c= c(0,0,0,0)
#'
#'
#' ## Using parameter coef and var.coef
#' saeHBLognormal <- Lognormal(formula,coef=c,var.coef=v,iter.update=10,data=dataLognormal)
#'
#' saeHBLognormal$Est                                 #Small Area mean Estimates
#' saeHBLognormal$refVar                              #Random effect variance
#' saeHBLognormal$coefficient                         #coefficient
#' #Load Library 'coda' to execute the plot
#' #autocorr.plot(saeHBLognormal$plot[[3]]) is used to generate ACF Plot
#' #plot(saeHBLognormal$plot[[3]]) is used to generate Density and trace plot
#'
#' ## Do not using parameter coef and var.coef
#' saeHBLognormal <- Lognormal(formula,data=dataLognormal)
#'
#'
#' ## For data with nonsampled area use dataLognormalNs
#'
#' }
#'
Lognormal <- function(formula,iter.update=3, iter.mcmc=10000, coef, var.coef, thin = 2, burn.in =2000, tau.u = 1, data){



  result <- list(Est = NA, refVar = NA, coefficient = NA,
                 plot=NA)


  formuladata <- model.frame(formula,data,na.action=NULL)
  if (any(is.na(formuladata[,-1])))
    stop("Auxiliary Variables contains NA values.")
  auxVar <- as.matrix(formuladata[,-1])
  nvar <- ncol(auxVar) + 1
  #formuladata <- data.frame(formuladata, n.samp = data[,n.samp])

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

    if (any(formuladata[,1]<=0)){
      stop("response variable must be " ,formula[2], " > 0")
    }

    x <- model.matrix(formula,data = as.data.frame(formuladata))
    n <- nrow(formuladata)
    mu.b = mu.b.value
    tau.b = tau.b.value
    tau.aa=tau.ab=tau.ba=tau.bb=1
    tau.ua=tau.ub=1
    a.var=1

    for (i in 1:iter.update){
      dat <- list("n"= n,  "nvar"= nvar, "y" = formuladata[,1], "x"=as.matrix(x[,-1]),
                  "mu.b"=mu.b, "tau.b"=tau.b,"tau.aa"=tau.aa,"tau.ab"=tau.ab,"tau.ba"=tau.ba,"tau.bb"=tau.bb,"tau.ua"=tau.ua,"tau.ub"=tau.ub)
      inits <- list(u = rep(0,n),b = mu.b, tau.u =tau.u)
      cat("model {
          for (i in 1:n) {
          		y[i] ~ dlnorm(mu.ln[i],tau[i])
							mu.ln[i] <- b[1] + sum(b[2:nvar]*x[i,]) + u[i]
							u[i] ~ dnorm(0,tau.u)
							mu[i] <- exp(mu.ln[i])
							tau[i] ~ dgamma(tau.a, tau.tb)
          }

          for (k in 1:nvar){
					    b[k] ~ dnorm(mu.b[k],tau.b[k])
          }
				  tau.a ~ dgamma(tau.aa, tau.ab)
				  tau.tb ~ dgamma(tau.ba, tau.bb)
					tau.u ~ dgamma(tau.ua, tau.ub)
					a.var <- 1 / tau.u
			}", file="lognormal.txt")

      jags.m <- jags.model(file = "lognormal.txt", data=dat, inits=inits, n.chains=1, n.adapt=500  )
      file.remove("lognormal.txt")
      params <- c("mu","a.var","b", "tau.u", "tau.a", "tau.tb")
      samps <- coda.samples( jags.m, params, n.iter=iter.mcmc, thin=thin)
      samps1 <- window(samps, start=burn.in+1, end=iter.mcmc)
      result_samps=summary(samps1)
      a.var=result_samps$statistics[1]
      beta=result_samps$statistics[2:(nvar+1),1:2]
      for (i in 1:nvar){
        mu.b[i]  = beta[i,1]
        tau.b[i] = 1/(beta[i,2]^2)
      }

      tau.aa =  result_samps$statistics[(n+nvar+2),1]^2/result_samps$statistics[(n+nvar+2),2]^2
      tau.ab =  result_samps$statistics[(n+nvar+2),1]/result_samps$statistics[(n+nvar+2),2]^2

      tau.ba =  result_samps$statistics[(n+nvar+3),1]^2/result_samps$statistics[(n+nvar+3),2]^2
      tau.bb =  result_samps$statistics[(n+nvar+3),1]/result_samps$statistics[(n+nvar+3),2]^2

      tau.ua = result_samps$statistics[(4+nvar+n),1]^2/result_samps$statistics[(4+nvar+n),2]^2
      tau.ub = result_samps$statistics[(4+nvar+n),1]/result_samps$statistics[(4+nvar+n),2]^2
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

    Quantiles <- as.data.frame(result_samps$quantiles[1:(3+nvar+n),])
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
    tau.aa=tau.ab=tau.ba=tau.bb=1
    tau.ua=tau.ub=1
    a.var=1

    formuladata$idx <- rep(1:n)
    data_sampled <- na.omit(formuladata)

    if (any(data_sampled[,1]<=0)){
      stop("response variable must be  " ,formula[2], " > 0")}


    data_nonsampled <- formuladata[-data_sampled$idx,]

    r=data_nonsampled$idx
    n1=nrow(data_sampled)
    n2=nrow(data_nonsampled)
    for (i in 1:iter.update){
      dat <- list("n1"= n1, "n2"=n2,"nvar"=nvar, "y_sampled" = data_sampled[,1],
                  "x_sampled"=as.matrix(data_sampled[,2:nvar]),
                  "x_nonsampled"=as.matrix(data_nonsampled[,2:nvar]),
                  "mu.b"=mu.b,"tau.b"=tau.b,
                  "tau.aa"=tau.aa,"tau.ab"=tau.ab,"tau.ba"=tau.ba,"tau.bb"=tau.bb,
                  "tau.ua"=tau.ua,"tau.ub"=tau.ub)
      inits <- list(u = rep(0,n1), v = rep(0,n2), b = mu.b, tau.u = tau.u)
      cat("model {
          for (i in 1:n1) {
          		y_sampled[i] ~ dlnorm(mu.ln[i],tau[i])
							mu.ln[i] <- b[1] + sum(b[2:nvar]*x_sampled[i,]) + u[i]
							u[i] ~dnorm(0,tau.u)
							mu[i] <- exp(mu.ln[i])
							tau[i] ~ dgamma(tau.a, tau.tb)
              }
          for (j in 1:n2) {

          		v[j] ~dnorm(0,tau.u)
							y_nonsampled[j] ~ dlnorm(mu.ln.nonsampled[j],tau.nonsampled[j])
							mu.ln.nonsampled[j] <- mu.b[1] + sum(mu.b[2:nvar]*x_nonsampled[j,]) +v[j]
							mu.nonsampled[j] <- exp(mu.ln.nonsampled[j])
							tau.nonsampled[j] ~ dgamma(tau.a, tau.tb)

					    }
					# prior
          for (k in 1:nvar){
					    b[k] ~ dnorm(mu.b[k],tau.b[k])
              }
				  tau.a ~ dgamma(tau.aa, tau.ab)
				  tau.tb ~ dgamma(tau.ba, tau.bb)
					tau.u ~ dgamma(tau.ua, tau.ub)
					a.var <- 1 / tau.u
			  }", file="LogNormal.txt")
      jags.m <- jags.model( file = "LogNormal.txt", data=dat, inits=inits, n.chains=1, n.adapt=500 )
      file.remove("LogNormal.txt")
      params <- c("mu","mu.nonsampled","a.var","b","tau.a","tau.tb","tau.u")
      samps <- coda.samples( jags.m, params, n.iter=iter.mcmc, thin=thin)
      samps1 <- window(samps, start=burn.in+1, end=iter.mcmc)
      result_samps=summary(samps1)
      a.var=result_samps$statistics[1]
      beta=result_samps$statistics[2:(nvar+1),1:2]
      for (i in 1:nvar){
        mu.b[i]  = beta[i,1]
        tau.b[i] = 1/(beta[i,2]^2)
      }
      tau.aa =  result_samps$statistics[(n+nvar+2),1]^2/result_samps$statistics[(n+nvar+2),2]^2
      tau.ab =  result_samps$statistics[(n+nvar+2),1]/result_samps$statistics[(n+nvar+2),2]^2

      tau.ba =  result_samps$statistics[(n+nvar+3),1]^2/result_samps$statistics[(n+nvar+3),2]^2
      tau.bb =  result_samps$statistics[(n+nvar+3),1]/result_samps$statistics[(n+nvar+3),2]^2

      tau.ua = result_samps$statistics[(4+nvar+n),1]^2/result_samps$statistics[(4+nvar+n),2]^2
      tau.ub = result_samps$statistics[(4+nvar+n),1]/result_samps$statistics[(4+nvar+n),2]^2

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

    mu=result_samps$statistics[(nvar+2):(1+nvar+n1),1:2]
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
