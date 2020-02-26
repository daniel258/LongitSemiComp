#' @title The function that simulates data form longitudinal bivariate binary semicompeting risks data with time-varying covarites
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
#' @param cens.poten.rate Potential cenosoring rate. At each time interval the probability of each alive observation to be censored
#' @return A list with two objects: \code{df.return} returns the data in a way similiat to counting process presentation, 
#' each unique \code{ID} has multiple rows, one for each interval. A time-fixed normally distributed random variable and 
#' a binary time-dependent covariate simulated as described in Nevo et al (\code{X}). The outcome data at each interval
#'  is given by \code{YNT} and \code{YT}. The second returned object is \code{cens}, a vector with per-person censoring indicator.
#'  This is not needed for the analysis as the data has a counting-process style representation, but it is useful for keeping
#'  track of the censoring rate when simulating data. 
#'  @author Daniel Nevo
#' @export
SimLongitDataTimeDep <- function(n.sample, times = 1:100,  beta.y,  alpha.nt, alpha.t, alpha.or, beta.nt, beta.t, beta.or,
                                 cens.poten.rate = 0) # cens.poten.rate is not really the censrate
{
  if (length(alpha.nt) != length(times)) {stop("alpha.nt should be in the same length as times")}
  if (length(alpha.t) != length(times)) {stop("alpha.t should be in the same length as times")}
  if (length(alpha.or) != length(times)) {stop("alpha.or should be in the same length as times")}
  J <- length(times)
  p <- length(beta.nt)
  X.time.fixed <- as.matrix(rnorm(n.sample))
  X.time.dep <- matrix(nrow = n.sample, ncol = J, 0) 
  X.time.dep[, 1] <- rbinom(n.sample, 1, 0.6)  # At baseline 0.6 are married
  risk.NT <- risk.T <- YNT <- YT <- matrix(nrow = n.sample, ncol = J, 0)
  risk.NT[,1] <- 1
  risk.T[,1] <- 1
  Xfirst <- cbind(X.time.fixed, X.time.dep[,1])
  p.nt.first <- expit(alpha.nt[1] + Xfirst%*%beta.nt)
  p.t.first <- expit(alpha.t[1] + Xfirst%*%beta.t)
  OR.first <- exp(alpha.or[1]+ Xfirst%*%beta.or)
  first.probs <- MargORtoJoint(p1marg = p.nt.first, p2marg = p.t.first, OR = OR.first)
  crude.data <- apply(X = first.probs, MARGIN = 1, FUN = sample, x = 1:4, size = 1, replace = T)
  YNT[,1] <- crude.data %in% c(2,4)
  YT[,1] <- crude.data %in% c(3,4)
  YNT[YNT[,1]==1, 2] <- 1
  YT[YT[,1]==1, 2] <- 1
  for(j in 2:J)
  {
    Xnow <- X.time.dep[ ,j - 1]
    Xnow[Xnow==1] <- rbinom(sum(Xnow), 1, 0.9) 
    X.time.dep[ ,j] <- Xnow
    Xtemp <- cbind(X.time.fixed, Xnow)
    risk.NT[,j] <- (YNT[, j - 1]==0 & YT[, j - 1]==0)
    risk.T[,j] <- 1*(YT[, j - 1]==0)
    at.risk.T.only <- risk.NT[, j]==0 & risk.T[, j]==1
    at.risk.both <- risk.NT[, j]==1 & risk.T[, j]==1
    # at risk for terminal event only
    if (sum(at.risk.T.only)>0)
    {
      probs.T.only <- expit(alpha.t[j] + Xtemp[at.risk.T.only, ]%*%beta.t + beta.y)
      YT[at.risk.T.only, j] <- rbinom(sum(at.risk.T.only), 1, probs.T.only) 
    }
    #at risk for both events
    if  (sum(at.risk.both)>0)
    {
      p.nt.both <- expit(alpha.nt[j] + Xtemp[at.risk.both, ]%*%beta.nt)
      p.t.both <- expit(alpha.t[j] + Xtemp[at.risk.both, ]%*%beta.t)
      OR.both <- exp(alpha.or[j]+ Xtemp[at.risk.both, ]%*%beta.or)
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
  # Add censoring
  C <- sample(x = 2:J, replace = T, size = n.sample)
  somecens <- rbinom(n.sample, 1, cens.poten.rate)
  cens <- rep(0, n.sample)
  for (i in 1:n.sample)
  {
  if (somecens[i]==1) {
  YNT[i, C[i]:J] <- YT[i, C[i]:J] <- NA
  risk.NT[i, C[i]:J] <- risk.T[i, C[i]:J] <- 0
  if (risk.T[i, C[i] - 1]==0) {cens[i] = 1} # cens=1 if observation i was actually censored
  }}
  obs.per.person <- apply(risk.T,1, function(x) sum(x==1, na.rm = T))
  ID <- rep(1:n.sample, times = obs.per.person)
  Xcln <- matrix(ncol = 2, nrow = sum(obs.per.person))
  TM <- YNTcln <- YTcln <- vector(length = sum(obs.per.person))
  temp.ind <- 1
  for (i in 1:n.sample)
  {
    nobs.i <- obs.per.person[i]
    indicesY <- temp.ind:(temp.ind + nobs.i - 1)
    YNTcln[indicesY] <- YNT[i, 1:nobs.i]
    YTcln[indicesY] <- YT[i, 1:nobs.i]
    TM[indicesY] <- 1:nobs.i
    Xcln[indicesY, 1] <- rep(X.time.fixed[i], nobs.i)
    Xcln[indicesY, 2] <- X.time.dep[i, 1:nobs.i]
    temp.ind <- temp.ind + nobs.i
  }
  df.return <- data.frame(X = Xcln, YNT = YNTcln, YT = YTcln, ID = ID, TM = TM) 
  return(list(df.return = df.return, cens = cens))
}

