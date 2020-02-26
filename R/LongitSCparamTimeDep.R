#' @title The function to fit a longitudinal bivariate binary model for semi-competing risks data with time-dependening 
#' covariates, and using unrestricted (per-interval parameter) time-varying functions.
#' @description The function implements the proposed methodology in Nevo et al. (2020+) for time-fixed covariates under possible 
#' right-censoring and left truncation. Data should be first converted to longitudinal bivariate binary representation. 
#' This could be done using \code{TimesToLongit}. This function uses B-splines representation the time-varying functions
#' and implements penalized maximum likelihood to fit the model.
#' @param data A data.frame or a list with columns named \code{ID}, \code{TM}   \code{YNT}, \code{YT} as well as all covariate 
#' names used in \code{formula.NT}, \code{formula.T} and \code{formula.OR}. See details below. Other names can be used for
#'  \code{YNT}, \code{YT}, but then their names should be given in the formulas below.
#' @param times Vector of increasing times (for example, the interval partition points \eqn{\tau_1,}..., \eqn{\tau_K}).
#' This vector is used to construct the B-splines
#' @param formula.NT A formula of the form \code{YNT ~ x1 + x2} where \code{x1} and \code{x2} are covariates to be used for 
#' for the non terminal probability sub-model.
#' @param formula.T A formula of the form \code{YT ~ x1 + x3} where \code{x1} and \code{x3} are covariates to be used for 
#' for the terminal probability sub-model.
#' @param formula.OR A formula of the form \code{ ~ x1 + x4} where \code{x1} and \code{x4} are covariates to be used for 
#' for the odds ratio sub-model.
#' @param epsOR How close it the OR allowed to be to one before assuming it equals to one. Default is \code{10^(-10)}
#' @param init Initial values for the parameters.
#' @param maxit.optim For internal use of \code{optim}. Default is 50000
#' @details For \code{data}, the represenation is similiar for the counting process represenation for survival data. 
#' \code{ID} identify each person, where \code{TM} identifes the intervals in which this person is under followup.
#'  \code{YNT} and \code{YT} indicate whetehr a non-terminal event and the terminal event occured by the end of interval \code{TM}.
#'  See examples in  \code{\link{SimLongitDataTimeDep}} 
#' @return The function returns an object of class \code{LongitSC} including estimates and confidence intervals for 
#' the time-varying functions and coefficients.
#' @note  For unrestricted baseline functions (no B-splines or penalization) use \code{\link{LongitSCparamTimeDep}}.
#' For time-fixed covariates use \code{\link{LongitSCtimeDep}}.
#'  @author Daniel Nevo
#'   @export
LongitSCparamTimeDep <- function(times = NULL, data, formula.NT, formula.T, 
                            formula.OR = NULL,  epsOR = 10^(-10),
                             init = NULL, maxit.optim = 50000)
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
  YNT <- model.response(model.frame(formula.NT, data = data))
  YT <- model.response(model.frame(formula.T, data = data))
  ID <- data$ID
  TM <- data$TM
  J <- length(unique(data$TM))
  if (is.null(times)) times <- sort(unique(data$TM))
  p <- pNT + pT + pOR
  n.params <- 1 + 3*J + p 
  if (is.null(init))
  {
    init <- rep(0.1,n.params)
  }
  res.opt <- optim(par = init, fn = ParamLogLikTimeDep, gr = GradParamLogLikTimeDep, 
                   hessian = T, method = "L-BFGS-B", control = list(maxit = 50000),
                   YT = YT, YNT = YNT, TM = TM, ID = ID,
                   XNT = XNTmat,  XT = XTmat, XOR = XORmat,
                   epsOR = epsOR)
    fit <- list()
    fit$formula.NT <- formula.NT
    fit$formula.T <- formula.T
    fit$formula.OR <- formula.OR
    fit$optim.conv <- res.opt$convergence
    fit$est <- res.opt$par
    fit$penal.lik <- -res.opt$value 
    fit$lik <- ParamLogLikTimeDep(param = fit$est, 
                           YT = YT, YNT = YNT, TM = TM,  
                           XNT = XNTmat,  XT = XTmat, XOR = XORmat,
                           ID = ID, epsOR = epsOR) # used for aic
    fit$hess.penal <- res.opt$hessian
    fit$se.naive <- sqrt(diag(solve(res.opt$hessian)))
    my.grad.sqrd <- GradParamLogLikPersTimeDep(param = res.opt$par,
                                               YT = YT, YNT = YNT, TM = TM,  
                                               XNT = XNTmat,  XT = XTmat, XOR = XORmat,
                                               ID = ID, epsOR = epsOR)
    fit$v.hat <- solve(res.opt$hessian)%*%my.grad.sqrd%*%(solve(res.opt$hessian))
    fit$se.rob <- sqrt(diag(fit$v.hat))
    hess.no.penal <- numDeriv::jacobian(func = GradParamLogLikTimeDep, x = res.opt$par, 
                                       YT = YT, YNT = YNT, TM = TM,  
                                       XNT = XNTmat,  XT = XTmat, 
                                       XOR = XORmat,
                                       ID = ID, epsOR = epsOR)
    if(!identical(dim(hess.no.penal), dim(res.opt$hessian))) {
      fit$df <- 0 # Indicates a problem
    } else {
      fit$df <- sum(diag((hess.no.penal%*%solve(res.opt$hessian))))
    }
    fit$aic <- -2*fit$lik - 2*fit$df # lik is minus the log-likelihood without the peanlty
    fit$coef.longterm <-  fit$est[1]
    fit$time.int.NT <- expit(fit$est[2:(1 + J)])
    fit$time.int.T <- expit(fit$est[(1 + J + 1):(1 + 2*J)])
    fit$time.int.OR <- exp(fit$est[(1 + 2*J + 1):(1 + 3*J)])
    fit$coef.NT <- fit$est[(1 + 3*J + 1):(1 + 3*J + pNT)]
    fit$coef.T <- fit$est[(1 + 3*J + pNT + 1):(1 + 3*J + pNT + pT)]
    fit$coef.OR <- fit$est[(1 + 3*J + pNT + pT + 1):(1 + 3*J + pNT + pT + pOR)]
    fit$se.longterm <- fit$se.rob[1]
    fit$se.rob.NT <- fit$se.rob[(1 + 3*J + 1):(1 + 3*J + pNT)]
    fit$se.rob.T <- fit$se.rob[(1 + 3*J + pNT + 1):(1 + 3*J + pNT + pT)]
    fit$se.rob.OR <- fit$se.rob[(1 + 3*J + pNT + pT + 1):(1 + 3*J + pNT + pT + pOR)]
  # calculate ci for the baseline prob.T1, prob.T2 and OR.
  # Using the appropoiate transformation of the CI for B%*%alpha
  sub.vhat.NT <- fit$v.hat[2:(1 + J), 2:(1 + J)]
  sub.vhat.T <- fit$v.hat[(1 + J + 1):(1 + 2*J), (1 + J + 1):(1 + 2*J)]
  sub.vhat.OR <- fit$v.hat[(1 + 2*J + 1):(1 + 3*J), (1 + 2*J + 1):(1 + 3*J)]
  fit$ci.L.alpha.NT <- expit(fit$est[2:(1 + J)] - qnorm(0.975)*sqrt(diag(sub.vhat.NT)))
  fit$ci.H.alpha.NT <- expit(fit$est[2:(1 + J)] + qnorm(0.975)*sqrt(diag(sub.vhat.NT)))
  fit$ci.L.alpha.T <- expit(fit$est[(1 + J + 1):(1 + 2*J)] - 
                              qnorm(0.975)*sqrt(diag(sub.vhat.T)))
  fit$ci.H.alpha.T <- expit(fit$est[(1 + J + 1):(1 + 2*J)] +
                           qnorm(0.975)*sqrt(diag(sub.vhat.T)))
  fit$ci.L.alpha.OR <- exp(fit$est[(1 + 2*J + 1):(1 + 3*J)] - 
                          qnorm(0.975)*sqrt(diag(sub.vhat.OR)))
  fit$ci.H.alpha.OR <- exp(fit$est[(1 + 2*J + 1):(1 + 3*J)] +
                          qnorm(0.975)*sqrt(diag(sub.vhat.OR)))
  if (pNT==1) {
    names(fit$coef.NT) <- names(fit$se.rob.NT) <- all.vars(formula.NT)} else {
      names(fit$coef.NT) <- names(fit$se.rob.NT) <- colnames(XNTmat)}
  if (pT==1) {
    names(fit$coef.T) <- names(fit$se.rob.T) <- all.vars(formula.T)} else {
      names(fit$coef.T) <- names(fit$se.rob.T) <- colnames(XTmat)}
  if (pOR==1) {
    names(fit$coef.OR) <- names(fit$se.rob.OR) <- all.vars(formula.OR)} else {
      names(fit$coef.OR) <- names(fit$se.rob.OR) <- colnames(XORmat)}
  class(fit) <- "LongitSC"
  fit
}
