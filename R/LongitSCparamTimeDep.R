#' @title The function that simulates data form longitudinal bivariate binary semicompeting risks data with time-varying covarites
#' and unrestricted (per-interval parametr) time-varying functions.
#' @description Given observed non-terminal and terminal event times, and time-depending covaraites, this function implments
#' theta proposed longitudinal bivariate binary representation according to given interval partition  in Nevo et al.
#' @param n.sample Desired sample size.
#' @param times A vector of increasing times. Normally of equal length
#' @param alpha.nt True value for \eqn{\alpha_1(t)} for each \eqn{t} in \code{times}.
#' @param alpha.t True value for \eqn{\alpha_2(t)} for each \eqn{t} in \code{times}.
#' @param alpha.or True value for \eqn{\alpha_\theta(t)} for each \eqn{t} in \code{times}.
#' @param beta.nt True value for \eqn{\beta_1}. 
#' @param beta.t True value for \eqn{\beta_2}.
#' @param beta.or True value for \eqn{\beta_\theta}.
#' @param beta.y True value for \eqn{\beta_y}.
#' @param cens.poten.rate Potential cenosoring rate. At each time interval the probability of each alive observation to be censored
#' @return A list with two objects: \code{df.return} returns the data in a way similiat to counting process presentation, 
#' each unique \code{ID} has multiple rows, one for each interval. A time-fixed normally distributed random variable and 
#' a binary time-dependent covariate simulated as described in Nevo et al (\code{X}). The outcome data at each interval
#'  is given by \code{YNT} and \code{YT}. The second returned object is \code{cens}, a vector with per-person censoring indicator.
#'  This is not needed for the analysis as the data has a counting-process style representation, but it is useful for keeping
#'  track of the censoring rate when simulating data. 
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
