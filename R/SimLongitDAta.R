#' @title The function that simulates data form longitudinal bivariate binary semicompeting risks data with baseline covariates
#' @description Given observed non-terminal and terminal event times, censoring indicators and possible left-truncation, 
#' this function returns the outcome data in the longitudinal bivariate binary representation according to given interval partition
#'  as proposed in Nevo et al.(2020+)
#' @param n.sample Desired sample size.
#' @param times A vector of increasing times. Normally of equal length
#' @param alpha.nt True value for \eqn{\alpha_1(t)} for each \eqn{t} in \code{times}.
#' @param alpha.t True value for \eqn{\alpha_2(t)} for each \eqn{t} in \code{times}.
#' @param alpha.or True value for \eqn{\alpha_\theta(t)} for each \eqn{t} in \code{times}.
#' @param beta.nt True value for \eqn{\beta_1}. 
#' @param beta.t True value for \eqn{\beta_2}. 
#' @param beta.or True value for \eqn{\beta_\theta}.
#' @param beta.y True value for \eqn{\beta_y}.
#' @param gamma.inter True value for \eqn{\gamma_[inter]}
#' @param X A matrix of time-fixed covariates to be used for simulating the data. Number of rows should be 
#' \code{n.sample} and 
#' number of columns should be equal to \code{length(beta.nt)}. If not specified, \code{X} is simulated 
#' as iid Gaussian random variables.
#' @return A list with the covariates used (\code{X}), at-risk indicators for each unit at each interval
#'  (\code{risk.NT} and \code{risk.T}) and outcome data
#'  at each interval (\code{YNT} and \code{YT}).
#'  
#' @examples
#' \dontrun{
#' # Simulate semicompeting risks data
#' set.seed(314)
#' times <- seq(1, 15, 1)
#' alpha.nt <- LongitSemiComp:::logit(dchisq(times,3, ncp =5)/2 + 0.025)
#' alpha.t <- LongitSemiComp:::logit(times*(0.075/10)  - 0.0005*(times/20)^2  + 0.05)
#' alpha.or <- 0.15 - times/10 + 0.75*(times/10)^2 + 0.3*(times/20)^3
#' plot(x = times, y= exp(alpha.or))
#' plot(x = times, y= LongitSemiComp:::expit(alpha.nt))
#' plot(x = times, y= LongitSemiComp:::expit(alpha.t))
#' beta.nt <- log(c(0.7, 3))
#' beta.t <- log(c(0.5, 1))
#' beta.or <- log(c(1, 1))
#' beta.y <- log(1.4)
#' my.data <- SimLongitData(n.sample = 200, times = times,  beta.y,  
#'                          alpha.nt, alpha.t, alpha.or, 
#'                          beta.nt, beta.t, beta.or)
#' longit.data <- my.data[-1]
#' X <- my.data[1]
#' }
#'
#' @author Daniel Nevo
#' @export
SimLongitData <- function(n.sample, times = 1:100,  beta.y,  alpha.nt, alpha.t, alpha.or, beta.nt, beta.t, beta.or, 
                          gamma.inter=NULL, X = NULL)
{
  if (length(alpha.nt) != length(times)) {stop("alpha.nt should be in the same length as times")}
  if (length(alpha.t) != length(times)) {stop("alpha.t should be in the same length as times")}
  if (length(alpha.or) != length(times)) {stop("alpha.or should be in the same length as times")}
  K <- length(times)
  p <- length(beta.nt)
  if (is.null(X)) {
  X <- matrix(nrow = n.sample, ncol = p, rnorm(n.sample*p))}
  if (!all(dim(X)==c(n.sample,p))) {stop("X should have n.sample rows and number of columns identical to length(beta.nt)")} 
  risk.NT <- risk.T <- YNT <- YT <- matrix(nrow = n.sample, ncol = K, 0)
  risk.NT[,1] <- 1
  risk.T[,1] <- 1
  p.nt.first <- expit(alpha.nt[1] + X%*%beta.nt)
  p.t.first <- expit(alpha.t[1] + X%*%beta.t)
  OR.first <- exp(alpha.or[1]+ X%*%beta.or)
  first.probs <- MargORtoJoint(p1marg = p.nt.first, p2marg = p.t.first, OR = OR.first)
  crude.data <- apply(X = first.probs, MARGIN = 1, FUN = sample, x = 1:4, size = 1, replace = T)
  YNT[,1] <- crude.data %in% c(2,4)
  YT[,1] <- crude.data %in% c(3,4)
  YNT[YNT[,1]==1, 2] <- 1
  YT[YT[,1]==1, 2] <- 1
  for(k in 2:K)
  {
    risk.NT[,k] <- (YNT[, k-1]==0 & YT[, k-1]==0)
    risk.T[,k] <- 1*(YT[, k-1]==0)
    at.risk.T.only <- risk.NT[, k]==0 & risk.T[, k]==1
    at.risk.both <- risk.NT[, k]==1 & risk.T[, k]==1
    # at risk for terminal event only
    if (sum(at.risk.T.only)>0)
    {
      if (is.null(gamma.inter)) {
        probs.T.only <- expit(alpha.t[k] + X[at.risk.T.only, ]%*%beta.t + beta.y)
      } else {
        probs.T.only <- expit(alpha.t[k] + X[at.risk.T.only, ]%*%beta.t + beta.y + X[at.risk.T.only, ]%*%gamma.inter)
        }
      YT[at.risk.T.only, k] <- rbinom(sum(at.risk.T.only), 1, probs.T.only)
    }
    #at risk for both events
    if  (sum(at.risk.both) > 0)
    {
      p.nt.both <- expit(alpha.nt[k] + X[at.risk.both, ]%*%beta.nt)
      p.t.both <- expit(alpha.t[k] + X[at.risk.both, ]%*%beta.t)
      OR.both <- exp(alpha.or[k]+ X[at.risk.both, ]%*%beta.or)
      probs.both <- MargORtoJoint(p1marg = p.nt.both, p2marg = p.t.both, OR = OR.both)
      crude.data.both <- apply(X = probs.both, MARGIN = 1, FUN = sample, x = 1:4, size = 1, replace = T)
      YNT[at.risk.both, k] <- crude.data.both %in% c(2, 4)
      YT[at.risk.both, k] <- crude.data.both %in% c(3, 4)
    }
    if(k < K)  {
      YNT[YNT[, k]==1, k + 1] <- 1
      YT[YT[, k]==1, k + 1] <- 1
    }
  }
  return(list(X = X, YNT = YNT, YT = YT, risk.NT = risk.NT, risk.T = risk.T))
}
