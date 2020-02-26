#' @title The function that simulates data form longitudinal bivariate binary semicompeting risks data with baseline covarites
#' @description Given observed non-terminal and terminal event times, censoring indicators and possible left-truncation, 
#' this function returns the outcome data in the longitudinal bivariate binary representation according to given interval partition
#'  as proposed in Nevo et al.
#' @param n.sample Desired sample size.
#' @param times A vector of increasing times. Normally of equal length
#' @param alpha.nt True value for \eqn{\alpha_1(t)} for each \eqn{t} in \code{times}.
#' @param alpha.t True value for \eqn{\alpha_2(t)} for each \eqn{t} in \code{times}.
#' @param alpha.or True value for \eqn{\alpha_\theta(t)} for each \eqn{t} in \code{times}.
#' @param beta.nt True value for \eqn{\beta_1}. 
#' @param beta.t True value for \eqn{\beta_2}. 
#' @param beta.or True value for \eqn{\beta_\theta}.
#' @param beta.y True value for \eqn{\beta_y}.
#' @param X A matrix of time-fixed covariates to be used for simulating the data. Number of rows should be 
#' \code{n.sample} and 
#' number of columns should be equal to \code{length(beta.nt)}. If not specfied, \code{X} is simulated 
#' as iid Gaussian random variables.
#' @return A list with the covariates used (\code{X}), at-risk indicators for each unit at each interval
#'  (\code{risk.NT} and \code{risk.T}) and outcome data
#'  at each interval (\code{YNT} and \code{YT}).
#'  @author Daniel Nevo
#' @export
SimLongitData <- function(n.sample, times = 1:100,  beta.y,  alpha.nt, alpha.t, alpha.or, beta.nt, beta.t, beta.or, X = NULL)
{
  if (length(alpha.nt) != length(times)) {stop("alpha.nt should be in the same length as times")}
  if (length(alpha.t) != length(times)) {stop("alpha.t should be in the same length as times")}
  if (length(alpha.or) != length(times)) {stop("alpha.or should be in the same length as times")}
  J <- length(times)
  p <- length(beta.nt)
  if (is.null(X)) {
  X <- matrix(nrow = n.sample, ncol = p, rnorm(n.sample*p))}
  if (!all(dim(X)==c(n.sample,p))) {stop("X should have n.sample rows and number of columns identical to length(beta.nt)")} 
  risk.NT <- risk.T <- YNT <- YT <- matrix(nrow = n.sample, ncol = J, 0)
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
  for(j in 2:J)
  {
    risk.NT[,j] <- (YNT[,j-1]==0 & YT[,j-1]==0)
    risk.T[,j] <- 1*(YT[,j-1]==0)
    at.risk.T.only <- risk.NT[,j]==0 & risk.T[,j]==1
    at.risk.both <- risk.NT[,j]==1 & risk.T[,j]==1
    # at risk for terminal event only
    if (sum(at.risk.T.only)>0)
    {
      probs.T.only <- expit(alpha.t[j] + X[at.risk.T.only, ]%*%beta.t + beta.y)
      YT[at.risk.T.only, j] <- rbinom(sum(at.risk.T.only), 1, probs.T.only)
    }
    #at risk for both events
    if  (sum(at.risk.both)>0)
    {
      p.nt.both <- expit(alpha.nt[j] + X[at.risk.both, ]%*%beta.nt)
      p.t.both <- expit(alpha.t[j] + X[at.risk.both, ]%*%beta.t)
      OR.both <- exp(alpha.or[j]+ X[at.risk.both, ]%*%beta.or)
      probs.both <- MargORtoJoint(p1marg = p.nt.both, p2marg = p.t.both, OR = OR.both)
      crude.data.both <- apply(X = probs.both, MARGIN = 1, FUN = sample, x = 1:4, size = 1, replace = T)
      YNT[at.risk.both, j] <- crude.data.both %in% c(2, 4)
      YT[at.risk.both, j] <- crude.data.both %in% c(3, 4)
    }
    if(j < J)  {
      YNT[YNT[, j]==1, j + 1] <- 1
      YT[YT[, j]==1, j + 1] <- 1
    }
  }
  return(list(X =X, YNT = YNT, YT = YT, risk.NT = risk.NT, risk.T = risk.T))
}
