#' @title The function to fit a longitudinal bivariate binary model for semi-competing risks data with time-depending
#' covariates, and using P-splines for the time-varying functions.
#' @description The function implements the proposed methodology in Nevo et al. (2020+) for time-depending covariates under 
#' possible right censoring and left truncation. Data should be first converted to longitudinal bivariate binary representation,
#' similar to the counting-process representation. See details below.
#' The \code{LongitSCtimeDep} function uses B-splines representation the time-varying functions
#' and implements penalized maximum likelihood to fit the model.
#' @param data A data.frame or a list with columns named \code{ID}, \code{TM},   \code{YNT}, \code{YT} as well as all covariate 
#' names used in \code{formula.NT}, \code{formula.T} and \code{formula.OR}. See details below. Other names can be used for
#' \code{YNT}, \code{YT}, but then their names should be given in the formulas below.
#' @param times Vector of increasing times (for example, the interval partition points \eqn{\tau_1,}..., \eqn{\tau_K}).
#' This vector is used to construct the B-splines
#' @param formula.NT A formula of the form \code{YNT ~ x1 + x2} where \code{x1} and \code{x2} are covariates to be used for 
#' for the non terminal probability sub-model.
#' @param formula.T A formula of the form \code{YT ~ x1 + x3} where \code{x1} and \code{x3} are covariates to be used for 
#' for the terminal probability sub-model.
#' @param formula.OR A formula of the form \code{ ~ x1 + x4} where \code{x1} and \code{x4} are covariates to be used for 
#' for the odds ratio sub-model.
#' @param epsOR How close it the OR allowed to be to one before assuming it equals to one. Default is \code{10^(-10)}
#' @param knots Number of knots for the B-splines.
#' @param lambda Penalization level for the penalized maximum likelihood. Either a vector of three values or a single number to be used for all three time-varying 
#' functions.
#' @param init Initial values for the parameters.
#' @param maxit.optim For internal use of \code{optim}. Default is 50000
#' @details For \code{data}, the representation is similar to the counting process representation of time-to-event data.  
#' \code{ID} identify each person, where \code{TM} identifies the intervals in which this person is under followup.
#' \code{YNT} and \code{YT} indicate whether a non-terminal event and the terminal event occurred by the end of interval \code{TM}.
#' See examples in  \code{\link{SimLongitDataTimeDep}}.
#' @return The function returns an object of class \code{LongitSC} including estimates and confidence intervals for 
#' the time-varying functions and coefficients.
#' @note  For unrestricted baseline functions (no B-splines or penalization) use \code{\link{LongitSCparamTimeDep}}.
#' For time-fixed covariates use \code{\link{LongitSCtimeDep}}.
#' 
#' @examples
#' \dontrun{
#' # Simulate semicompeting risks data
#' set.seed(314)
#' times <- seq(1,15,1)
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
#' sim.data <- SimLongitDataTimeDep(n.sample = 1500, times = times,  beta.y,  cens.poten.rate = 0.5,
#'                                  alpha.nt, alpha.t, alpha.or, 
#'                                  beta.nt, beta.t, beta.or)
#' # Censoring rate
#' mean(sim.data$cens)
#' my.df <- sim.data$df.return
#' # Analysis
#' res <- LongitSCtimeDep(data = my.df, times = times,  formula.NT = YNT ~ X.1 + X.2, 
#'                        formula.T = YT  ~ X.1 + X.2, 
#'                        formula.OR = ~ X.1 + X.2, 
#'                        knots = 5, lambda = 1)
#'  res
#' }
#'
#' @author Daniel Nevo
#' @export
LongitSCtimeDep <- function(times = NULL, data, formula.NT, formula.T, 
                            formula.OR = NULL, epsOR = 10^(-10),
                     knots = NULL, lambda = 0, init = NULL, maxit.optim = 50000)
{
  if (is.null(knots)) knots <- 5
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
  YNT <- model.response(model.frame(formula.NT, data = data))
  YT <- model.response(model.frame(formula.T, data = data))
  ID <- data$ID
  TM <- data$TM
  if(length(lambda)==1) lambda <- rep(lambda, 3)
  if(length(lambda)!=3) stop("lambda should be of length 1 or 3")
  K <- length(unique(data$TM))
  if (is.null(times)) times <- sort(unique(data$TM))
  smooth.aux <- mgcv::smooth.construct.ps.smooth.spec(mgcv::s(times,bs="ps", k = knots), data = list(times = times),
                                                      knots = list(times = c(min(times),max(times))))
  S.penal <- smooth.aux$S[[1]]
  Bsplines <- smooth.aux$X
  Q <- ncol(Bsplines)
  p <- pNT + pT + pOR
  n.params <- 1 + 3*Q + p 
  if (is.null(init))
  {
    init <- rep(0.1,n.params)
  }
  res.opt <- tryCatch(optim(par = init, fn = PenalLogLikTimeDep, gr = GradPenalLogLikTimeDep, 
                   hessian = T, method = "L-BFGS-B", control = list(maxit = maxit.optim),
                   YT = YT, YNT = YNT, TM = TM, ID = ID, 
                   XNT = XNTmat,  XT = XTmat, XOR = XORmat,
                   TimeBase = Bsplines, TimePen = S.penal, 
                   lambda = lambda, epsOR = epsOR), error=function(e) {e})
  if(inherits(res.opt, "error")){
    res.opt <- optim(par = init, fn = PenalLogLikTimeDep, gr = GradPenalLogLikTimeDep, 
                              hessian = T, method = "BFGS", control = list(maxit = maxit.optim),
                              YT = YT, YNT = YNT, TM = TM, ID = ID, 
                              XNT = XNTmat,  XT = XTmat, XOR = XORmat,
                              TimeBase = Bsplines, TimePen = S.penal, 
                              lambda = lambda, epsOR = epsOR)
    }
  fit <- list()
  fit$formula.NT <- formula.NT
  fit$formula.T <- formula.T
  fit$formula.OR <- formula.OR
  fit$Bsplines <- Bsplines
  fit$optim.conv <- res.opt$convergence
  fit$est <- res.opt$par
  fit$penal.lik <- -res.opt$value 
  fit$lik <- PenalLogLikTimeDep(param = fit$est, 
                         YT = YT, YNT = YNT, TM = TM,  
                         XNT = XNTmat,  XT = XTmat, XOR = XORmat,
                         ID = ID,
                         TimeBase = Bsplines, epsOR = epsOR, 
                         TimePen = S.penal, lambda = rep(0, 3)) # used for aic
  fit$hess.penal <- res.opt$hessian
  hess.inv <- tryCatch(solve(res.opt$hessian), error=function(e) {e})
  err.hess.inv <- inherits(hess.inv, "error")
  if(err.hess.inv){
    hess.inv <- tryCatch(ginv(res.opt$hessian), error=function(e) {e})
    err.hess.inv <- inherits(hess.inv, "error")}
  if (!err.hess.inv)
  {
  fit$se.naive <- sqrt(diag(hess.inv))
  my.grad.sqrd <- GradPenalLogLikPersTimeDep(param = res.opt$par,
                                               YT = YT, YNT = YNT, TM = TM,  
                                               XNT = XNTmat,  XT = XTmat, XOR = XORmat,
                                               ID = ID,
                                               TimeBase = Bsplines, epsOR = epsOR, 
                                               TimePen = S.penal, lambda = rep(0,3))
  fit$v.hat <- hess.inv%*%my.grad.sqrd%*%hess.inv
  fit$se.rob <- sqrt(diag(fit$v.hat))
  hess.no.penal <- numDeriv::jacobian(func = GradPenalLogLikTimeDep, x = res.opt$par, 
                                       YT = YT, YNT = YNT, TM = TM,  
                                       XNT = XNTmat,  XT = XTmat, 
                                       XOR = XORmat,
                                       ID = ID,
                                       TimeBase = Bsplines, epsOR = epsOR, 
                                       TimePen = S.penal, lambda = rep(0,3))
  if(!identical(dim(hess.no.penal), dim(res.opt$hessian))) {
    fit$df <- 0 # Indicates a problem
  } else {
    fit$df <- sum(diag((hess.no.penal%*%hess.inv)))
  }
  fit$aic <- -2*fit$lik - 2*fit$df # lik is minus the log-likelihood without the peanlty
  fit$coef.longterm <-  fit$est[1]
  fit$time.int.NT <- expit(Bsplines%*%fit$est[2:(1 + Q)])
  fit$time.int.T <- expit(Bsplines%*%fit$est[(1 + Q + 1):(1 + 2*Q)])
  fit$time.int.OR <- exp(Bsplines%*%fit$est[(1 + 2*Q + 1):(1 + 3*Q)])
  fit$coef.NT <- fit$est[(1 + 3*Q + 1):(1 + 3*Q + pNT)]
  fit$coef.T <- fit$est[(1 + 3*Q + pNT + 1):(1 + 3*Q + pNT + pT)]
  fit$coef.OR <- fit$est[(1 + 3*Q + pNT + pT + 1):(1 + 3*Q + pNT + pT + pOR)]
  fit$se.longterm <- fit$se.rob[1]
  fit$se.rob.NT <- fit$se.rob[(1 + 3*Q + 1):(1 + 3*Q + pNT)]
  fit$se.rob.T <- fit$se.rob[(1 + 3*Q + pNT + 1):(1 + 3*Q + pNT + pT)]
  fit$se.rob.OR <- fit$se.rob[(1 + 3*Q + pNT + pT + 1):(1 + 3*Q + pNT + pT + pOR)]
  # calculate ci for the baseline prob.T1, prob.T2 and OR.
  # Using the appropoiate transformation of the CI for B%*%alpha
  sub.vhat.NT <- fit$v.hat[2:(1 + Q), 2:(1 + Q)]
  sub.vhat.T <- fit$v.hat[(1 + Q + 1):(1 + 2*Q), (1 + Q + 1):(1 + 2*Q)]
  sub.vhat.OR <- fit$v.hat[(1 + 2*Q + 1):(1 + 3*Q), (1 + 2*Q + 1):(1 + 3*Q)]
  fit$ci.L.alpha.NT <- expit(Bsplines%*%fit$est[2:(1 + Q)] - 
                           qnorm(0.975)*sqrt(diag(Bsplines%*%sub.vhat.NT%*%t(Bsplines))))
  fit$ci.H.alpha.NT <- expit(Bsplines%*%fit$est[2:(1 + Q)] +
                           qnorm(0.975)*sqrt(diag(Bsplines%*%sub.vhat.NT%*%t(Bsplines))))
  fit$ci.L.alpha.T <- expit(Bsplines%*%fit$est[(1 + Q + 1):(1 + 2*Q)] - 
                           qnorm(0.975)*sqrt(diag(Bsplines%*%sub.vhat.T%*%t(Bsplines))))
  fit$ci.H.alpha.T <- expit(Bsplines%*%fit$est[(1 + Q + 1):(1 + 2*Q)] +
                           qnorm(0.975)*sqrt(diag(Bsplines%*%sub.vhat.T%*%t(Bsplines))))
  fit$ci.L.alpha.OR <- exp(Bsplines%*%fit$est[(1 + 2*Q + 1):(1 + 3*Q)] - 
                          qnorm(0.975)*sqrt(diag(Bsplines%*%sub.vhat.OR%*%t(Bsplines))))
  fit$ci.H.alpha.OR <- exp(Bsplines%*%fit$est[(1 + 2*Q + 1):(1 + 3*Q)] +
                          qnorm(0.975)*sqrt(diag(Bsplines%*%sub.vhat.OR%*%t(Bsplines))))
  if (pNT==1) {
    names(fit$coef.NT) <- names(fit$se.rob.NT) <- all.vars(formula.NT)} else {
      names(fit$coef.NT) <- names(fit$se.rob.NT) <- colnames(XNTmat)}
  if (pT==1) {
    names(fit$coef.T) <- names(fit$se.rob.T) <- all.vars(formula.T)} else {
      names(fit$coef.T) <- names(fit$se.rob.T) <- colnames(XTmat)}
  if (pOR==1) {
    names(fit$coef.OR) <- names(fit$se.rob.OR) <- all.vars(formula.OR)} else {
      names(fit$coef.OR) <- names(fit$se.rob.OR) <- colnames(XORmat)}
  } else {
    fit$coef.longterm <-  fit$est[1]
    fit$time.int.NT <- expit(Bsplines%*%fit$est[2:(1 + Q)])
    fit$time.int.T <- expit(Bsplines%*%fit$est[(1 + Q + 1):(1 + 2*Q)])
    fit$time.int.OR <- exp(Bsplines%*%fit$est[(1 + 2*Q + 1):(1 + 3*Q)])
    fit$coef.NT <- fit$est[(1 + 3*Q + 1):(1 + 3*Q + pNT)]
    fit$coef.T <- fit$est[(1 + 3*Q + pNT + 1):(1 + 3*Q + pNT + pT)]
    fit$coef.OR <- fit$est[(1 + 3*Q + pNT + pT + 1):(1 + 3*Q + pNT + pT + pOR)]
    if (pNT==1) {
      names(fit$coef.NT) <- all.vars(formula.NT)} else {
        names(fit$coef.NT) <- colnames(XNTmat)}
    if (pT==1) {
      names(fit$coef.T) <-  all.vars(formula.T)} else {
        names(fit$coef.T) <- colnames(XTmat)}
    if (pOR==1) {
      names(fit$coef.OR) <- all.vars(formula.OR)} else {
        names(fit$coef.OR) <-  colnames(XORmat)}
  } 
  class(fit) <- "LongitSC"
  fit
}
