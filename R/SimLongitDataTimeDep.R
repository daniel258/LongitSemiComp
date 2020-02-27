#' @title The function that simulates data form longitudinal bivariate binary semicompeting risks data with time-depending covariates
#' @description Given observed non-terminal and terminal event times, censoring indicators and possible left-truncation, 
#' this function returns the outcome data in the longitudinal bivariate binary representation according to given interval partition
#'  as proposed in Nevo et al. (2020+)
#' @param n.sample Desired sample size.
#' @param times A vector of increasing times. Normally of equal length
#' @param alpha.nt True value for \eqn{\alpha_1(t)} for each \eqn{t} in \code{times}.
#' @param alpha.t True value for \eqn{\alpha_2(t)} for each \eqn{t} in \code{times}.
#' @param alpha.or True value for \eqn{\alpha_\theta(t)} for each \eqn{t} in \code{times}.
#' @param beta.nt True value for \eqn{\beta_1}. 
#' @param beta.t True value for \eqn{\beta_2}.
#' @param beta.or True value for \eqn{\beta_\theta}.
#' @param beta.y True value for \eqn{\beta_y}.
#' @param cens.poten.rate Potential censoring rate. At each time interval the probability of each alive observation to be censored
#' @return A list with two objects: \code{df.return} returns the data in a way similar to counting process presentation, 
#' each unique \code{ID} has multiple rows, one for each interval. A time-fixed normally distributed random variable and 
#' a binary time-dependent covariate simulated as described in Nevo et al (\code{X}). The outcome data at each interval
#' is given by \code{YNT} and \code{YT}. The second returned object is \code{cens}, a vector with per-person censoring indicator.
#' This is not needed for the analysis as the data has a counting-process style representation, but it is useful for keeping
#' track of the censoring rate when simulating data. 
#'@examples
#'\dontrun{
#' set.seed(314)
#' times <- seq(1,14,1)
#' alpha.nt <- LongitSemiComp:::logit(times*0.005  + 0.005*(times-2)^2 - 
#'                                    (0.0002*(times + 1)^3) + 0.005)
#' alpha.t <- LongitSemiComp:::logit(times*0.0075  + 0.001*(times^2)  + 0.03)
#' alpha.or <- 0.9 + 0.175*times - 0.02*times^2 #+ 0.3*(times/20)^3
#' alpha.or[times >= 13] <- 0
#' plot(x = times, y= exp(alpha.or))
#' plot(x = times, y= LongitSemiComp:::expit(alpha.nt))
#' plot(x = times, y= LongitSemiComp:::expit(alpha.t))
#' beta.nt <- log(c(0.7, 3))
#' beta.t <- log(c(0.5, 1))
#' beta.or <- log(c(1, 1))
#' beta.y <- log(1.4)
#' my.data <- SimLongitDataTimeDep(n.sample = 2000, times = times,  beta.y,  
#'                                 alpha.nt, alpha.t, alpha.or, 
#'                                 beta.nt, beta.t, beta.or, cens.poten.rate = 0.5)
#' df.data <- my.data$df.return
#' mean(my.data$cens)
#'}
#'
#' @author Daniel Nevo
#' @export
SimLongitDataTimeDep <- function(n.sample, times = 1:100,  beta.y,  alpha.nt, alpha.t, alpha.or, beta.nt, beta.t, beta.or,
                                 cens.poten.rate = 0) # cens.poten.rate is not really the censrate
{
  if (length(alpha.nt) != length(times)) {stop("alpha.nt should be in the same length as times")}
  if (length(alpha.t) != length(times)) {stop("alpha.t should be in the same length as times")}
  if (length(alpha.or) != length(times)) {stop("alpha.or should be in the same length as times")}
  K <- length(times)
  p <- length(beta.nt)
  X.time.fixed <- as.matrix(rnorm(n.sample))
  X.time.dep <- matrix(nrow = n.sample, ncol = K, 0) 
  X.time.dep[, 1] <- rbinom(n.sample, 1, 0.6)  # At baseline Pr(X(t)=1)=0.6
  risk.NT <- risk.T <- YNT <- YT <- matrix(nrow = n.sample, ncol = K, 0)
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
  for(k in 2:K)
  {
    Xnow <- X.time.dep[ ,k - 1]
    Xnow[Xnow==1] <- rbinom(sum(Xnow), 1, 0.9) 
    X.time.dep[ ,k] <- Xnow
    Xtemp <- cbind(X.time.fixed, Xnow)
    risk.NT[,k] <- (YNT[, k - 1]==0 & YT[, k - 1]==0)
    risk.T[,k] <- 1*(YT[, k - 1]==0)
    at.risk.T.only <- risk.NT[, k]==0 & risk.T[, k]==1
    at.risk.both <- risk.NT[, k]==1 & risk.T[, k]==1
    # at risk for terminal event only
    if (sum(at.risk.T.only)>0)
    {
      probs.T.only <- expit(alpha.t[k] + Xtemp[at.risk.T.only, ]%*%beta.t + beta.y)
      YT[at.risk.T.only, k] <- rbinom(sum(at.risk.T.only), 1, probs.T.only) 
    }
    #at risk for both events
    if  (sum(at.risk.both) > 0)
    {
      p.nt.both <- expit(alpha.nt[k] + Xtemp[at.risk.both, ]%*%beta.nt)
      p.t.both <- expit(alpha.t[k] + Xtemp[at.risk.both, ]%*%beta.t)
      OR.both <- exp(alpha.or[k]+ Xtemp[at.risk.both, ]%*%beta.or)
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
  # Add censoring
  C <- sample(x = 2:K, replace = T, size = n.sample)
  somecens <- rbinom(n.sample, 1, cens.poten.rate)
  cens <- rep(0, n.sample)
  for (i in 1:n.sample)
  {
  if (somecens[i]==1) {
  YNT[i, C[i]:K] <- YT[i, C[i]:K] <- NA
  risk.NT[i, C[i]:K] <- risk.T[i, C[i]:K] <- 0
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

