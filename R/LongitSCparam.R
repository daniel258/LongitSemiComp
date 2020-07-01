#' @title The function to fit a longitudinal bivariate binary model for semi-competing risks data with time-fixed covariates
#' and using unstrucutred representation for the time-varying functions.
#' @description The function implements the proposed methodology in Nevo et al. (2020+) for time-fixed covariates under possible
#' right censoring and left truncation. Data should be first converted to longitudinal bivariate binary representation.
#' This could be done using \code{\link{TimesToLongit}}. The \code{LongitSCparam} function implements maximum likelihood to fit the model.
#' @param longit.data A list with entries named \code{risk.NT}, \code{risk.T},  \code{YNT}, \code{YT}. See details below.
#' @param times A vector of increasing times (for example, the interval partition points \eqn{\tau_1}... \eqn{\tau_K}).
#' This vector is used to construct the B-splines
#' @param formula.NT A formula of the form \code{ ~ x1 + x2} where \code{x1} and \code{x2} are covariates to be used for
#' the non-terminal event probability sub-model.
#' @param formula.T A formula of the form \code{ ~ x1 + x3} where \code{x1} and \code{x3} are covariates to be used for
#' for the non-terminal event probability sub-model.
#' @param formula.OR A formula of the form \code{ ~ x1 + x4} where \code{x1} and \code{x4} are covariates to be used for
#' for the odds ratio sub-model.
#' @param  formula.inter A formula of the form \code{ ~ x1 + x4} where \code{x1} and \code{x4} are covariates to be used for 
#' the interavtion between non-terminal event and covariates sub-model.
#' @param data A data.frame with the covariates specified in \code{formula.NT}, \code{formula.T} and \code{formula.OR}.
#' @param epsOR How close it the OR allowed to be to one before assuming it equals to one. Default is \code{10^(-10)}
#' @param init Initial values for the parameters.
#' @param maxit.optim For internal use of \code{optim}. Default is 50000.
#' @details Each of the matrices in \code{longit.data},  \code{risk.NT}, \code{risk.T},  \code{YNT}, \code{YT} have a row for each
#' unit and a column for each interval. Then, \code{risk.NT} and \code{risk.T} indicate whether the unit is at risk in each interval
#' for each of the events, respectively.  The matrices \code{YNT} and \code{YT} indicate whether the non-terminal and terminal
#' event, respectively, were observed by the end of each interval.
#' The function \code{\link{TimesToLongit}} can be used to obtain this representation of semicompeting risks time-to-event data.
#' @return The function returns an object of class \code{LongitSC} including estimates and confidence intervals for
#' the time-varying functions and coefficients.
#' @note For time-varying covariates use \code{\link{LongitSCtimeDep}} or \code{\link{LongitSCparamTimeDep}}.
#' @examples
#' \dontrun{
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
#' beta.or <- log(c(1.4, 1))
#' beta.y <- log(1.4)
#' my.data <- SimLongitData(n.sample = 2000, times = times,  beta.y,
#'                          alpha.nt, alpha.t, alpha.or,
#'                          beta.nt, beta.t, beta.or)
#' longit.data <- my.data[-1]
#' X <- as.data.frame(my.data[1])
#' LongitSC(longit.data = longit.data, times = times,  formula.NT =  ~ X.1 + X.2,
#'          formula.T =  ~ X.1 + X.2,
#'          formula.OR = ~ X.1 + X.2,
#'          data = X, epsOR = 10^(-10),
#'          knots = 5, lambda = 1)
#' }
#' @author Daniel Nevo
#' @export
LongitSCparam <- function(longit.data, times = NULL, formula.NT, formula.T, formula.OR = NULL, formula.inter = ~ 1, data,
                               epsOR = 10^(-10), init = NULL, maxit.optim = 50000)
{
  if (!is.null(formula.NT)) {
    XNTmat <- as.matrix(model.matrix(formula.NT, data = data)[, -1])
    pNT <- ncol(XNTmat)
  } else {
    XNTmat <- NULL
    pNT <- 0
  }
  if (!is.null(formula.T)) {
    XTmat <- as.matrix(model.matrix(formula.T, data = data)[, -1])
    pT <- ncol(XTmat)
  } else {
    XTmat <- NULL
    pT <- 0
  }
  if (!is.null(formula.OR)) {
    XORmat <- as.matrix(model.matrix(formula.OR, data = data)[, -1])
    pOR <- ncol(XORmat)
  } else {
    XORmat <- NULL
    pOR <- 0
  }
  if (!is.null(formula.inter)) {
    XinterMat <- as.matrix(model.matrix(formula.inter, data = data)) # we need the intercept for beta_y
    pInter <- ncol(XinterMat)
  } else {
    XinterMat <- NULL
    pInter <- 0
  }
  K <- ncol(longit.data$YNT)
  if (is.null(times)) times <- 1:K
  p <- pNT + pT + pOR + pInter
  n.params <- 3*K + p
  if (is.null(init))
  {
    init <- rep(-0.1,n.params)
  }
  res.opt <- optim(par = init, fn = ParamLogLik, gr = GradParamLogLik, hessian = T,
                   control = list(maxit = maxit.optim), method = "L-BFGS-B", epsOR = epsOR,
                   XNT = XNTmat, XT = XTmat, XOR = XORmat, XinterMat = XinterMat,
                   YNT = longit.data$YNT, YT = longit.data$YT,
                   riskNT = longit.data$risk.NT, riskT = longit.data$risk.T)
  fit <- list()
  fit$formula.NT <- formula.NT
  fit$formula.T <- formula.T
  fit$formula.OR <- formula.OR
  fit$formula.inter <- formula.inter
  fit$optim.conv <- res.opt$convergence
  fit$est <- res.opt$par
  fit$lik <- -res.opt$value
  fit$hess <- res.opt$hessian
  fit$v.hat <- solve(res.opt$hessian)
  fit$se <- sqrt(diag(fit$v.hat))
  if (any(is.nan(fit$se)))
  {
    fit$se[is.nan(fit$se)] <- 100
    warning("close to infinite SE for some parameters, set SE=100")
  }
  fit$time.int.NT <- expit(fit$est[1:K])
  fit$time.int.T <- expit(fit$est[(K + 1):(2*K)])
  fit$time.int.OR <- exp(fit$est[(2*K + 1):(3*K)])
  fit$coef.NT <- fit$est[(3*K + 1):(3*K + pNT)]
  fit$coef.T <- fit$est[(3*K + pNT + 1):(3*K + pNT + pT)]
  fit$coef.OR <- fit$est[(3*K + pNT + pT + 1):(3*K + pNT + pT + pOR)]
  fit$coef.inter <- fit$est[(3*K + pNT + pT + pOR + 1):(3*K + pNT + pT + pOR + pInter)]
  fit$se.coef.NT <- fit$se[(3*K + 1):(3*K + pNT)]
  fit$se.coef.T <- fit$se[(3*K + pNT + 1):(3*K + pNT + pT)]
  fit$se.coef.OR <- fit$se[(3*K + pNT + pT + 1):(3*K + pNT + pT + pOR)]
  fit$se.coef.inter <- fit$se[(3*K + pNT + pT + pOR + 1):(3*K + pNT + pT + pOR + pInter)]
  # calculate ci for the baseline prob.T1, prob.T2 and OR.
  # Using the appropriate transformation of the CI for B%*%eta
  sub.vhat.NT <- fit$v.hat[1:K, 1:K]
  sub.vhat.T <- fit$v.hat[(K + 1):(2*K), (K + 1):(2*K)]
  sub.vhat.OR <- fit$v.hat[(2*K + 1):(3*K), (2*K + 1):(3*K)]
  fit$ci.L.alpha.NT <- expit(fit$est[1:K] - qnorm(0.975)* sqrt(diag(sub.vhat.NT)))
  fit$ci.H.alpha.NT <- expit(fit$est[1:K] + qnorm(0.975)* sqrt(diag(sub.vhat.NT)))
  fit$ci.L.alpha.T <- expit(fit$est[(K + 1):(2*K)] - qnorm(0.975) * sqrt(diag(sub.vhat.T)))
  fit$ci.H.alpha.T <- expit(fit$est[(K + 1):(2*K)] + qnorm(0.975) * sqrt(diag(sub.vhat.T)))
  fit$ci.L.alpha.OR <- exp(fit$est[(2*K + 1):(3*K)] - qnorm(0.975) * sqrt(diag(sub.vhat.OR)))
  fit$ci.H.alpha.OR <- exp(fit$est[(2*K + 1):(3*K)] + qnorm(0.975) * sqrt(diag(sub.vhat.OR)))
  if (pNT==1) {
    names(fit$coef.NT) <- names(fit$se.coef.NT) <- all.vars(formula.NT)} else {
      names(fit$coef.NT) <- names(fit$se.coef.NT) <- colnames(XNTmat)}
  if (pT==1) {
    names(fit$coef.T) <- names(fit$se.coef.T) <- all.vars(formula.T)} else {
      names(fit$coef.T) <- names(fit$se.coef.T) <- colnames(XTmat)}
  if (pOR==1) {
    names(fit$coef.OR) <- names(fit$se.coef.OR) <- all.vars(formula.OR)} else {
      names(fit$coef.OR) <- names(fit$se.coef.OR) <- colnames(XORmat)}
  if (pInter==1) {
    names(fit$coef.inter) <- names(fit$se.coef.inter) <- all.vars(formula.inter)} else {
      names(fit$coef.inter) <- names(fit$se.coef.inter) <- paste0("NT-event:",colnames(XinterMat))}
  names(fit$coef.inter)[1] <- names(fit$se.coef.inter)[1] <- "NT-event"
  class(fit) <- "LongitSC"
  fit
}
