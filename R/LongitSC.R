#' @export
LongitSC <- function(longit.data, times = NULL, formula.NT, formula.T, formula.OR = NULL, 
                     formula.inter.gamma = NULL,
                     data, epsOr = 10^(-10),
                     knots = NULL, lambda = 0, init = NULL, maxit.optim = 50000)
{
 # m <- match.call()
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
  if (!is.null(formula.inter.gamma)) {
    XinterMat <- as.matrix(model.matrix(formula.inter.gamma, data = data)[, -1])
    p.inter.gamma <- ncol(XinterMat)
  } else {
    XinterMat <- NULL
    p.inter.gamma <- 0
  }
  # XNTmat <- as.matrix(model.matrix(formula.NT, data = data)[, -1])
  # XTmat <- as.matrix(model.matrix(formula.T, data = data)[, -1])
  # XORmat <- as.matrix(model.matrix(formula.OR, data = data)[, -1])
  if(length(lambda)==1) lambda <- rep(lambda, 3)
  if(length(lambda)!=3) stop("lambda should be of length 1 or 3")
  J <- ncol(longit.data$YNT)
  if (is.null(times)) times <- 1:J
  smooth.aux <- mgcv::smooth.construct.ps.smooth.spec(mgcv::s(times,bs="ps", k = knots), data = list(times = times),
                                                      knots = list(times = c(min(times),max(times))))
  S.penal <- smooth.aux$S[[1]]
  Bsplines <- smooth.aux$X
  Q <- ncol(Bsplines)
  # pNT <- ncol(XNTmat)
  # pT <- ncol(XTmat)
  # pOR <- ncol(XORmat)
  
  p <- pNT + pT + pOR
  n.params <- 1 + 3*Q + p + p.inter.gamma
  if (is.null(init))
  {
    init <- rep(0,n.params)
  }
  if (pOR==0) {
    res.opt <- optim(par = init, fn = PenalLogLikNullModelOR, gr = GradPenalLogLikNullModelOR, hessian = T,
                     control = list(maxit = maxit.optim), method = "L-BFGS-B", epsOR = epsOr,
                     XNT = XNTmat, XT = XTmat, #XOR = XORmat, 
                     YNT = longit.data$YNT, YT = longit.data$YT,
                     riskNT = longit.data$risk.NT, riskT = longit.data$risk.T,
                     TimeBase = Bsplines,
                     TimePen = S.penal, lambda = lambda)
    fit <- list()
    fit$formula.NT <- formula.NT
    fit$formula.T <- formula.T
    fit$formula.OR <- formula.OR
    fit$formula.inter.gamma <- formula.inter.gamma
    fit$Bsplines <- Bsplines
    fit$optim.conv <- res.opt$convergence
    fit$est <- res.opt$par
    fit$penal.lik <- -res.opt$value 
    fit$lik <- PenalLogLikNullModelOR(param = fit$est, epsOR = epsOr, XNT = XNTmat, XT = XTmat, #XOR = XORmat,
                           YNT = longit.data$YNT, YT = longit.data$YT, 
                           riskNT = longit.data$risk.NT, riskT = longit.data$risk.T,
                           TimeBase = Bsplines,
                           TimePen = S.penal, lambda = rep(0,3)) # used for aic
    fit$hess.penal <- res.opt$hessian
    fit$se.naive <- sqrt(diag(solve(res.opt$hessian)))
    my.grad.sqrd <- GradPenalLogLikPersNullModelOR(param = res.opt$par, epsOR = epsOr,
                                        XNT = XNTmat, XT = XTmat, #XOR = XORmat,
                                        YNT = longit.data$YNT, YT = longit.data$YT,
                                        riskNT = longit.data$risk.NT, riskT = longit.data$risk.T,
                                        TimeBase = Bsplines,
                                        TimePen = S.penal, lambda = lambda)
    fit$v.hat <- solve(res.opt$hessian)%*%my.grad.sqrd%*%(solve(res.opt$hessian))
    fit$se.rob <- sqrt(diag(fit$v.hat))
    hess.no.penal <- numDeriv::hessian(func = PenalLogLikNullModelOR, x = res.opt$par, epsOR = epsOr,
                                       XNT = XNTmat, XT = XTmat,# XOR = XORmat,
                                       YNT = longit.data$YNT, YT = longit.data$YT, riskNT = longit.data$risk.NT, riskT = longit.data$risk.T,
                                       TimeBase = Bsplines,  TimePen = S.penal, lambda = 0)
    # hess.penal <- numDeriv::hessian(func = PenalLogLik, x = res.opt$par, epsOR = epsOr,XNT = XNTmat, XT = XTmat, XOR = XORmat,
    #                                 YNT = longit.data$YNT, YT = longit.data$YT, riskNT = longit.data$risk.NT, riskT = longit.data$risk.T,
    #                                 TimeBase = Bsplines,  TimePen = S.penal, lambda = lambda)
    if(!identical(dim(hess.no.penal), dim(res.opt$hessian))) {
      fit$df <- 0 # Indicates a problem
    } else {
      fit$df <- sum(diag((hess.no.penal%*%solve(res.opt$hessian))))
    }
    fit$aic <- -2*fit$lik - 2*fit$df # lik is minus the log-likelihood without the peanlty
    fit$coef.longterm <-  fit$est[1]
    fit$time.int.NT <- expit(Bsplines%*%fit$est[2:(1 + Q)])
    fit$time.int.T <- expit(Bsplines%*%fit$est[(1 + Q + 1):(1 + 2*Q)])
    fit$time.int.OR <- exp(Bsplines%*%fit$est[(1 + 2*Q + 1):(1 + 3*Q)])
    fit$coef.NT <- fit$est[(1 + 3*Q + 1):(1 + 3*Q + pNT)]
    fit$coef.T <- fit$est[(1 + 3*Q + pNT + 1):(1 + 3*Q + pNT + pT)]
    fit$coef.OR <- NULL
    #fit$coef.OR <- fit$est[(1 + 3*Q + pNT + pT + 1):(1 + 3*Q + pNT + pT + pOR)]
    fit$se.longterm <- fit$se.rob[1]
    fit$se.rob.NT <- fit$se.rob[(1 + 3*Q + 1):(1 + 3*Q + pNT)]
    fit$se.rob.T <- fit$se.rob[(1 + 3*Q + pNT + 1):(1 + 3*Q + pNT + pT)]
    fit$se.rob.OR <- NULL
    #fit$se.rob.OR <- fit$se.rob[(1 + 3*Q + pNT + pT + 1):(1 + 3*Q + pNT + pT + pOR)]
    # calculate ci for the baseline prob.T1, prob.T2 and OR.
    # Using the appropoiate transformation of the CI for B%*%alpha
  } else {
    if (p.inter.gamma>0)
    {
  res.opt <- optim(par = init, fn = PenalLogLikWithInter, gr = GradPenalLogLikWithInter, hessian = T,
                   control = list(maxit = maxit.optim), method = "L-BFGS-B", epsOR = epsOr,
                   XNT = XNTmat, XT = XTmat, XOR = XORmat, XinterMat = XinterMat,
                   YNT = longit.data$YNT, YT = longit.data$YT,
                   riskNT = longit.data$risk.NT, riskT = longit.data$risk.T,
                   TimeBase = Bsplines,
                   TimePen = S.penal, lambda = lambda)
  fit <- list()
  fit$formula.NT <- formula.NT
  fit$formula.T <- formula.T
  fit$formula.OR <- formula.OR
  fit$formula.inter.gamma <- formula.inter.gamma
  fit$Bsplines <- Bsplines
  fit$optim.conv <- res.opt$convergence
  fit$est <- res.opt$par
  fit$penal.lik <- -res.opt$value 
  fit$lik <- PenalLogLikWithInter(param = fit$est, epsOR = epsOr, XNT = XNTmat, XT = XTmat, XOR = XORmat,
                         XinterMat = XinterMat,
                         YNT = longit.data$YNT, YT = longit.data$YT, 
                         riskNT = longit.data$risk.NT, riskT = longit.data$risk.T,
                         TimeBase = Bsplines,
                         TimePen = S.penal, lambda = rep(0,3)) # used for aic (this is minus the log-lik!!)
  fit$hess.penal <- res.opt$hessian
  fit$se.naive <- sqrt(diag(solve(res.opt$hessian)))
  my.grad.sqrd <- GradPenalLogLikPersWithInter(param = res.opt$par, epsOR = epsOr,
                                      XNT = XNTmat, XT = XTmat, XOR = XORmat,
                                      XinterMat = XinterMat,
                                      YNT = longit.data$YNT, YT = longit.data$YT,
                                      riskNT = longit.data$risk.NT, riskT = longit.data$risk.T,
                                      TimeBase = Bsplines,
                                      TimePen = S.penal, lambda = lambda)
  fit$v.hat <- solve(res.opt$hessian)%*%my.grad.sqrd%*%(solve(res.opt$hessian))
  fit$se.rob <- sqrt(diag(fit$v.hat))
  hess.no.penal <- numDeriv::hessian(func = PenalLogLikWithInter, x = res.opt$par, epsOR = epsOr,
                                     XNT = XNTmat, XT = XTmat, XOR = XORmat,
                                     XinterMat = XinterMat,
                                     YNT = longit.data$YNT, YT = longit.data$YT, riskNT = longit.data$risk.NT, riskT = longit.data$risk.T,
                                     TimeBase = Bsplines,  TimePen = S.penal, lambda = 0)
  if(!identical(dim(hess.no.penal), dim(res.opt$hessian))) {
    fit$df <- 0 # Indicates a problem
  } else {
    fit$df <- sum(diag((hess.no.penal%*%solve(res.opt$hessian))))
  }
  fit$aic <- -2*fit$lik - 2*fit$df # lik is minus the log-likelihood without the peanlty
  fit$coef.longterm <-  fit$est[1]
  fit$time.int.NT <- expit(Bsplines%*%fit$est[2:(1 + Q)])
  fit$time.int.T <- expit(Bsplines%*%fit$est[(1 + Q + 1):(1 + 2*Q)])
  fit$time.int.OR <- exp(Bsplines%*%fit$est[(1 + 2*Q + 1):(1 + 3*Q)])
  fit$coef.NT <- fit$est[(1 + 3*Q + 1):(1 + 3*Q + pNT)]
  fit$coef.T <- fit$est[(1 + 3*Q + pNT + 1):(1 + 3*Q + pNT + pT)]
  fit$coef.OR <- fit$est[(1 + 3*Q + pNT + pT + 1):(1 + 3*Q + pNT + pT + pOR)]
  fit$coef.inter.longterm <- fit$est[(1 + 3*Q + pNT + pT + pOR + 1):(1 + 3*Q + pNT + pT + pOR + p.inter.gamma)]
  fit$se.longterm <- fit$se.rob[1]
  fit$se.rob.NT <- fit$se.rob[(1 + 3*Q + 1):(1 + 3*Q + pNT)]
  fit$se.rob.T <- fit$se.rob[(1 + 3*Q + pNT + 1):(1 + 3*Q + pNT + pT)]
  fit$se.rob.OR <- fit$se.rob[(1 + 3*Q + pNT + pT + 1):(1 + 3*Q + pNT + pT + pOR)]
  fit$se.rob.longterm.inter <- fit$se.rob[(1 + 3*Q + pNT + pT + pOR + 1):(1 + 3*Q + pNT + pT + pOR + p.inter.gamma)]
    }
    else {
    res.opt <- optim(par = init, fn = PenalLogLik, gr = GradPenalLogLik, hessian = T,
                     control = list(maxit = maxit.optim), method = "L-BFGS-B", epsOR = epsOr,
                     XNT = XNTmat, XT = XTmat, XOR = XORmat,# XinterMat = XinterMat,
                     YNT = longit.data$YNT, YT = longit.data$YT,
                     riskNT = longit.data$risk.NT, riskT = longit.data$risk.T,
                     TimeBase = Bsplines,
                     TimePen = S.penal, lambda = lambda)
    fit <- list()
    fit$formula.NT <- formula.NT
    fit$formula.T <- formula.T
    fit$formula.OR <- formula.OR
    fit$formula.inter.gamma <- formula.inter.gamma
    fit$Bsplines <- Bsplines
    fit$optim.conv <- res.opt$convergence
    fit$est <- res.opt$par
    fit$penal.lik <- -res.opt$value 
    fit$lik <- PenalLogLik(param = fit$est, epsOR = epsOr, XNT = XNTmat, XT = XTmat, XOR = XORmat,
                           #XinterMat = XinterMat,
                           YNT = longit.data$YNT, YT = longit.data$YT, 
                           riskNT = longit.data$risk.NT, riskT = longit.data$risk.T,
                           TimeBase = Bsplines,
                           TimePen = S.penal, lambda = rep(0,3)) # used for aic
    fit$hess.penal <- res.opt$hessian
    fit$se.naive <- sqrt(diag(solve(res.opt$hessian)))
    my.grad.sqrd <- GradPenalLogLikPers(param = res.opt$par, epsOR = epsOr,
                                        XNT = XNTmat, XT = XTmat, XOR = XORmat,
                                #        XinterMat = XinterMat,
                                        YNT = longit.data$YNT, YT = longit.data$YT,
                                        riskNT = longit.data$risk.NT, riskT = longit.data$risk.T,
                                        TimeBase = Bsplines,
                                        TimePen = S.penal, lambda = lambda)
    fit$v.hat <- solve(res.opt$hessian)%*%my.grad.sqrd%*%(solve(res.opt$hessian))
    fit$se.rob <- sqrt(diag(fit$v.hat))
    hess.no.penal <- numDeriv::hessian(func = PenalLogLik, x = res.opt$par, epsOR = epsOr,
                                       XNT = XNTmat, XT = XTmat, XOR = XORmat,
                               #        XinterMat = XinterMat,
                                       YNT = longit.data$YNT, YT = longit.data$YT, riskNT = longit.data$risk.NT, riskT = longit.data$risk.T,
                                       TimeBase = Bsplines,  TimePen = S.penal, lambda = 0)
    if(!identical(dim(hess.no.penal), dim(res.opt$hessian))) {
      fit$df <- 0 # Indicates a problem
    } else {
      fit$df <- sum(diag((hess.no.penal%*%solve(res.opt$hessian))))
    }
    fit$aic <- -2*fit$lik - 2*fit$df # lik is minus the log-likelihood without the peanlty
    fit$coef.longterm <-  fit$est[1]
    fit$time.int.NT <- expit(Bsplines%*%fit$est[2:(1 + Q)])
    fit$time.int.T <- expit(Bsplines%*%fit$est[(1 + Q + 1):(1 + 2*Q)])
    fit$time.int.OR <- exp(Bsplines%*%fit$est[(1 + 2*Q + 1):(1 + 3*Q)])
    fit$coef.NT <- fit$est[(1 + 3*Q + 1):(1 + 3*Q + pNT)]
    fit$coef.T <- fit$est[(1 + 3*Q + pNT + 1):(1 + 3*Q + pNT + pT)]
    fit$coef.OR <- fit$est[(1 + 3*Q + pNT + pT + 1):(1 + 3*Q + pNT + pT + pOR)]
 #   fit$coef.inter.longterm <- fit$est[(1 + 3*Q + pNT + pT + pOR + 1):(1 + 3*Q + pNT + pT + pOR + p.inter.gamma)]
    fit$se.longterm <- fit$se.rob[1]
    fit$se.rob.NT <- fit$se.rob[(1 + 3*Q + 1):(1 + 3*Q + pNT)]
    fit$se.rob.T <- fit$se.rob[(1 + 3*Q + pNT + 1):(1 + 3*Q + pNT + pT)]
    fit$se.rob.OR <- fit$se.rob[(1 + 3*Q + pNT + pT + 1):(1 + 3*Q + pNT + pT + pOR)]
#    fit$se.rob.longterm.inter <- fit$se.rob[(1 + 3*Q + pNT + pT + pOR + 1):(1 + 3*Q + pNT + pT + pOR + p.inter.gamma)]
  }}
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
  # fit$ci.base.NT <- c(ci.L.alpha.NT, ci.H.alpha.NT)
  # fit$ci.base.T <- c(ci.L.alpha.T, ci.H.alpha.T)
  # fit$ci.base.OR <- c(ci.L.alpha.OR, ci.H.alpha.OR)
  if (pNT==1) {
    names(fit$coef.NT) <- names(fit$se.rob.NT) <- all.vars(formula.NT)} else {
      names(fit$coef.NT) <- names(fit$se.rob.NT) <- colnames(XNTmat)}
  if (pT==1) {
    names(fit$coef.T) <- names(fit$se.rob.T) <- all.vars(formula.T)} else {
      names(fit$coef.T) <- names(fit$se.rob.T) <- colnames(XTmat)}
  if (pOR==1) {
    names(fit$coef.OR) <- names(fit$se.rob.OR) <- all.vars(formula.OR)} else {
      names(fit$coef.OR) <- names(fit$se.rob.OR) <- colnames(XORmat)}
  if (p.inter.gamma==1) {
    names(fit$coef.inter.longterm) <- names(fit$se.rob.longterm.inter) <- all.vars(formula.inter.gamma)} else {
      names(fit$coef.inter.longterm) <- names(fit$se.rob.longterm.inter) <- colnames(XinterMat)}
  class(fit) <- "LongitSC"
  fit
}
