#' @title Small Area Estimation using Hierarchical Bayesian under Beta Distribution
#' @description This function is implemented to variable of interest \eqn{(y)} that assumed to be a Beta Distribution. The range of data must be \eqn{0<y<1}. The data proportion is supposed to be implemented with this function.
#' @param formula Formula that describe the fitted model
#' @param iter.update Number of updates with default \code{3}
#' @param iter.mcmc Number of total iterations per chain with default \code{2000}
#' @param thin Thinning rate, must be a positive integer with default \code{1}
#' @param burn.in Number of iterations to discard at the begining with default \code{1000}
#' @param data The data frame
#'
#' @return This function returns a list of the following objects:
#'    \item{Est}{A vector with the values of Small Area mean Estimates using Hierarchical bayesian method }
#'    \item{sd}{A vector with the values of Standard deviation of Small Area Mean Estimates using Hierarchical bayesian method}
#'    \item{refVar}{Estimated random effect variances}
#'    \item{coefficient}{A dataframe with the estimated model coefficient}
#'    \item{plot}{Trace, Dencity, Autocorrelation Function Plot of MCMC samples}
#'
#' @export Beta
#'
#' @examples
#'
#' \dontrun{
#' ##Load Dataset
#' data(dataBeta)
#'
#'
#' #Compute Fitted Model
#' #y ~ x1 +x2
#'
#'
#' ## For data without any nonsampled area
#' ##Load Dataset
#' data(dataBeta)
#' saeHBbeta <- Beta(formula = y ~ x1 +x2, data = dataBeta )
#' saeHBbeta$Est                                 #Small Area mean Estimates
#' saeHBbeta$sd                                  #Standard deviation of Small Area Mean Estimates
#' saeHBbeta$refVar                              #refVar
#' saeHBbeta$coefficient                         #coefficient
#' autocorr.plot(saeHBbeta$plot[[3]])            #ACF Plot
#' plot(saeHBbeta$plot[[3]])                     #Dencity and trace plot
#'
#' ## For data with nonsampled area
#' ##Load Dataset
#' data(dataBetaNs)
#' saeHBbetaNs <- Beta(formula = y ~ x1 +x2,  data = dataBetaNs )
#' saeHBbetaNs$Est                                 #Small Area mean Estimates
#' saeHBbetaNs$sd                                  #Standard deviation of Small Area Mean Estimates
#' saeHBbetaNs$refVar                              #refVar
#' saeHBbetaNs$coefficient                         #coefficient
#' autocorr.plot(saeHBbetaNs$plot[[3]])            #ACF Plot
#' plot(saeHBbetaNs$plot[[3]])                     #Dencity and trace plot
#'
#' }
#'
#'
#'
Beta <- function(formula, iter.update=3, iter.mcmc=2000, thin = 1, burn.in =1000, data){

  result <- list(Est = NA, sd = NA, refVar = NA, coefficient = NA,
                 plot=NA)

  formuladata <- model.frame(formula,data,na.action=NULL)
  if (any(is.na(formuladata[,-1])))
    stop("Auxiliary Variables contains NA values.")


  if (!any(is.na(formuladata[,1]))){
    formuladata <- as.matrix(na.omit(formuladata))

    if (any(formuladata[,1]==0) || any(formuladata[,1]==1)){
      stop("response variable must be 0 < " ,formula[2], " < 1")
    }
    x <- model.matrix(formula,data = as.data.frame(formuladata))
    x <- as.matrix(x)
    n <- nrow(formuladata)
    nvar <- ncol(formuladata)
    mu.b =rep(0,nvar)
    tau.b = rep(1,nvar)
    tau.ub = tau.ua = phi.a = phi.b = a.var = 1

    for(iter in 1:iter.update){

      dat <- list("n"=n,"nvar"=nvar,"y"=formuladata[,1],"x"=as.matrix(formuladata[,2:nvar]),"mu.b"=mu.b,
                  "tau.b"=tau.b,"tau.ua"=tau.ua, "tau.ub"=tau.ub, "phi.a"=phi.a, "phi.b"=phi.b)

      inits <- list(phi=1, u=rep(0,n), b = mu.b, tau.u = 1)



      cat("model{
					for (i in 1:n) {
							y[i] ~ dbeta(A[i],B[i])
							A[i] <- mu[i] * phi
              B[i] <- (1-mu[i]) * phi
              logit(mu[i]) <- b[1] + sum(b[2:nvar]*x[i,]) + u[i]
							u[i] ~ dnorm(0,tau.u)
					}

					for (k in 1:nvar){
					    b[k] ~ dnorm(mu.b[k],tau.b[k])
					}

					phi ~ dgamma(phi.a,phi.b)
					tau.u ~ dgamma(tau.ua,tau.ub)
					a.var <- 1 / tau.u

			}", file="saeHBbeta.txt")

      jags.m <- jags.model(file = "saeHBbeta.txt", data=dat, inits=inits, n.chains=1, n.adapt=500)

      params <- c("mu","a.var","b","phi","tau.u")
      samps  <- coda.samples( jags.m, params, n.iter=iter.mcmc, thin=thin)
      samps1 <- window(samps,start = burn.in+1, end =iter.mcmc )
      result_samps=summary(samps1)
      a.var=result_samps$statistics[1]
      beta=result_samps$statistics[2:(nvar+1),1:2]
      for (i in 1:nvar){
        mu.b[i]  = beta[i,1]
        tau.b[i] = 1/(beta[i,2]^2)
      }

      phi.a = result_samps$statistics[2+nvar+n,1]^2/result_samps$statistics[2+nvar+n,2]^2
      phi.a = result_samps$statistics[2+nvar+n,1]/result_samps$statistics[2+nvar+n,2]^2
      tau.ua = result_samps$statistics[3+nvar+n,1]^2/result_samps$statistics[3+nvar+n,2]^2
      tau.ub = result_samps$statistics[3+nvar+n,1]/result_samps$statistics[3+nvar+n,2]^2


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
    colnames(Estimation)=c("mean","sd")

    Quantiles <- as.data.frame(result_samps$quantiles[1:(2+nvar+n),])
    q_beta <- (Quantiles[2:(nvar+1),])
    rownames(q_beta) <- b.varnames
    beta <- cbind(beta,q_beta)

  }

  else {
    formuladata <- as.data.frame(formuladata)
    n <- nrow(formuladata)
    nvar <- ncol(formuladata)
    mu.b =rep(0,nvar)
    tau.b = rep(1,nvar)
    tau.ub = tau.ua = phi.a = phi.b = a.var = 1

    formuladata$idx <- rep(1:n)
    data_sampled <- na.omit(formuladata)

    if (any(data_sampled[,1]==0) || any(data_sampled[,1]==1)){
      stop("response variable must be 0 < " ,formula[2], " < 1")}

    data_nonsampled <- formuladata[-data_sampled$idx,]
    r=data_nonsampled$idx
    n1=nrow(data_sampled)
    n2=nrow(data_nonsampled)



    for(iter in 1:iter.update){
      dat <- list("n1" = n1,"n2"=n2,"nvar"=nvar,"y_sampled" = data_sampled[,1], "x_sampled"=as.matrix(data_sampled[,2:nvar]),
                  "x_nonsampled"=as.matrix(data_nonsampled[,2:nvar]),"mu.b"=mu.b,
                  "tau.b"=tau.b, "tau.ua"=tau.ua,"tau.ub"=tau.ub, "phi.a"=phi.a,"phi.b"=phi.b)

      inits <- list(phi=1,u = rep(0,n1),v = rep(0,n2), b = mu.b, tau.u = 1)

      cat("model {
					for (i in 1:n1) {
							y_sampled[i] ~ dbeta(A[i],B[i])
							A[i] <- mu[i] * phi
              B[i] <- (1-mu[i]) * phi
              logit(mu[i]) <- b[1] +  sum(b[2:nvar]*x_sampled[i,]) + u[i]
							u[i] ~ dnorm(0,tau.u)

					}

					for (j in 1:n2) {
					v[j]~dnorm(0,tau.u)
            logit(mu.nonsampled[j]) <- b[1] +sum(b[2:nvar]*x_nonsampled[j,])+v[j]
          }

					for (k in 1:nvar){
					    b[k] ~ dnorm(mu.b[k],tau.b[k])
					}

					phi ~ dgamma(phi.a,phi.b)
					tau.u ~ dgamma(tau.ua,tau.ub)
					a.var <- 1 / tau.u

			}", file="saeHBbeta.txt")

      jags.m <- jags.model( file = "saeHBbeta.txt", data=dat, inits=inits, n.chains=1, n.adapt=500 )
      params <- c("mu","mu.nonsampled","a.var","b", "phi", "tau.u")
      samps <- coda.samples( jags.m, params, n.iter=iter.mcmc, thin=thin)
      samps1 <- window(samps,start = burn.in+1, end =iter.mcmc )
      result_samps=summary(samps1)
      a.var=result_samps$statistics[1]
      beta=result_samps$statistics[2:(nvar+1),1:2]
      for (i in 1:nvar){
        mu.b[i]  = beta[i,1]
        tau.b[i] = 1/(beta[i,2]^2)
      }
      phi.a = result_samps$statistics[2+nvar+n,1]^2/result_samps$statistics[2+nvar+n,2]^2
      phi.a = result_samps$statistics[2+nvar+n,1]/result_samps$statistics[2+nvar+n,2]^2
      tau.ua = result_samps$statistics[3+nvar+n,1]^2/result_samps$statistics[3+nvar+n,2]^2
      tau.ub = result_samps$statistics[3+nvar+n,1]/result_samps$statistics[3+nvar+n,2]^2
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
    colnames(Estimation)=c("mean","sd")

    Quantiles <- as.data.frame(result_samps$quantiles[1:(2+nvar+n),])
    q_beta <- (Quantiles[2:(nvar+1),])
    rownames(q_beta) <- b.varnames
    beta <- cbind(beta,q_beta)
  }

  result$Est = Estimation$mean
  result$sd         = Estimation$sd
  result$refVar      = a.var
  result$coefficient = beta
  result$plot       = list(graphics.off() ,par(mar=c(2,2,2,2)),autocorr.plot(result_mcmc,col="brown2",lwd=2),plot(result_mcmc,col="brown2",lwd=2))
  return(result)


}
