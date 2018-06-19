LongitSC <- function(longit.data, times = NULL, formula.NT, formula.T, formula.OR, data, epsOr = 10^(-10),
                     knots = NULL, lambda = 0, init = NULL, maxit.optim = 50000)
{
  m <- match.call()
  if (is.null(knots)) m$knots <- 5
  XNTmat <- model.matrix(formula.NT, data = data)[, -1] %>% as.matrix
  XTmat <- model.matrix(formula.T, data = data)[, -1] %>% as.matrix
  XORmat <- model.matrix(formula.OR, data = data)[, -1] %>% as.matrix
  if(length(lambda)==1) m$lambda <- rep(lambda, 3)
  m$J <- ncol(longit.data$YNT)
  if (is.null(times)) times <- 1:m$J
  smooth.aux <- mgcv::smooth.construct.ps.smooth.spec(mgcv::s(times,bs="ps", k = m$knots), data = list(times = times),
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
  fit$Bsplines <- Bsplines
  fit$optim.conv <- res.opt$convergence
  fit$est <- res.opt$par
  fit$se.naive <- sqrt(diag(solve(res.opt$hessian)))
  fit$formula.NT <- formula.NT
  fit$formula.T <- formula.T
  fit$formula.OR <- formula.OR
  my.grad.sqrd <- GradPenalLogLikPers(param = res.opt$par, epsOR = epsOr,
                                       XNT = XNTmat, XT = XTmat, XOR = XORmat,
                                       YNT = longit.data$YNT, YT = longit.data$YT,
                                       riskNT = longit.data$risk.NT, riskT = longit.data$risk.T,
                                       TimeBase = Bsplines,
                                       TimePen = S.penal, lambda = lambda)
  fit$v.hat <- solve(res.opt$hessian)%*%my.grad.sqrd%*%(solve(res.opt$hessian))
  fit$se.rob <- sqrt(diag(fit$v.hat))
  fit$coef.longterm <-  fit$est[1]
  fit$time.int.NT <- expit(Bsplines%*%fit$est[2:(1 + Q)])
  fit$time.int.T <- expit(Bsplines%*%fit$est[(1 + Q + 1):(1 + 2*Q)])
  fit$time.int.OR <- exp(Bsplines%*%fit$est[(1 + 2*Q + 1):(1 + 3*Q)])
  fit$coef.NT <- fit$est[(1 + 3*Q + 1):(1 + 3*Q + pNT)]
  fit$coef.T <- fit$est[(1 + 3*Q + pNT + 1):(1 + 3*Q + pNT + pT)]
  fit$coef.OR <- fit$est[(1 + 3*Q + pNT + pT + 1):(1 + 3*Q + pNT + pT + pOR)]
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
