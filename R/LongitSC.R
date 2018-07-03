#' @export
LongitSC <- function(longit.data, times = NULL, formula.NT, formula.T, formula.OR, data, epsOr = 10^(-10),
                     knots = NULL, lambda = 0, init = NULL, maxit.optim = 50000)
{
 # m <- match.call()
  if (is.null(knots)) knots <- 5
  XNTmat <- as.matrix(model.matrix(formula.NT, data = data)[, -1])
  XTmat <- as.matrix(model.matrix(formula.T, data = data)[, -1])
  XORmat <- as.matrix(model.matrix(formula.OR, data = data)[, -1])
  if(length(lambda)==1) lambda <- rep(lambda, 3)
  if(length(lambda)!=3) stop("lambda should be of length 1 or 3")
  J <- ncol(longit.data$YNT)
  if (is.null(times)) times <- 1:J
  smooth.aux <- mgcv::smooth.construct.ps.smooth.spec(mgcv::s(times,bs="ps", k = knots), data = list(times = times),
                                                      knots = list(times = c(min(times),max(times))))
  S.penal <- smooth.aux$S[[1]]
  Bsplines <- smooth.aux$X
  Q <- ncol(Bsplines)
  pNT <- ncol(XNTmat)
  pT <- ncol(XTmat)
  pOR <- ncol(XORmat)
  p <- pNT + pT + pOR
  n.params <- 1 + 3*Q + p
  if (is.null(init))
  {
    init <- rep(0,n.params)
  }
  res.opt <- optim(par = init, fn = PenalLogLik, gr = GradPenalLogLik, hessian = T,
                   control = list(maxit = maxit.optim), method = "L-BFGS-B", epsOR = epsOr,
                   XNT = XNTmat, XT = XTmat, XOR = XORmat,
                   YNT = longit.data$YNT, YT = longit.data$YT,
                   riskNT = longit.data$risk.NT, riskT = longit.data$risk.T,
                   TimeBase = Bsplines,
                   TimePen = S.penal, lambda = lambda)
  fit <- list()
  fit$formula.NT <- formula.NT
  fit$formula.T <- formula.T
  fit$formula.OR <- formula.OR
  fit$Bsplines <- Bsplines
  fit$optim.conv <- res.opt$convergence
  fit$est <- res.opt$par
  fit$penal.lik <- -res.opt$value 
  fit$lik <- PenalLogLik(param = fit$est, XNT = XNTmat, XT = XTmat, XOR = XORmat,
                         YNT = longit.data$YNT, YT = longit.data$YT,
                         riskNT = longit.data$risk.NT, riskT = longit.data$risk.T,
                         TimeBase = Bsplines,
                         TimePen = S.penal, lambda = rep(0,3)) # used for aic
  fit$hess.penal <- res.opt$hessian
  fit$se.naive <- sqrt(diag(solve(res.opt$hessian)))
  my.grad.sqrd <- GradPenalLogLikPers(param = res.opt$par, epsOR = epsOr,
                                      XNT = XNTmat, XT = XTmat, XOR = XORmat,
                                      YNT = longit.data$YNT, YT = longit.data$YT,
                                      riskNT = longit.data$risk.NT, riskT = longit.data$risk.T,
                                      TimeBase = Bsplines,
                                      TimePen = S.penal, lambda = lambda)
  fit$v.hat <- solve(res.opt$hessian)%*%my.grad.sqrd%*%(solve(res.opt$hessian))
  fit$se.rob <- sqrt(diag(fit$v.hat))
  hess.no.penal <- numDeriv::hessian(func = PenalLogLik, x = res.opt$par, epsOR = epsOr,XNT = XNTmat, XT = XTmat, XOR = XORmat,
                                     YNT = longit.data$YNT, YT = longit.data$YT, riskNT = longit.data$risk.NT, riskT = longit.data$risk.T,
                                     TimeBase = Bsplines,  TimePen = S.penal, lambda = 0)
  # hess.penal <- numDeriv::hessian(func = PenalLogLik, x = res.opt$par, epsOR = epsOr,XNT = XNTmat, XT = XTmat, XOR = XORmat,
  #                                 YNT = longit.data$YNT, YT = longit.data$YT, riskNT = longit.data$risk.NT, riskT = longit.data$risk.T,
  #                                 TimeBase = Bsplines,  TimePen = S.penal, lambda = lambda)
  if(!identical(dim(hess.no.penal), dim(res.opt$hessian))) {
    fit$df <- 0 # Indicate a problem
  } else {
    fit$df <- sum(diag((hess.no.penal%*%solve(res.opt$hessian))))
  }
  fit$aic <- 2*fit$lik - 2*fit$df # Need to fix that and rerun simulations - this includes the penalty and I need the acutal lik
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
