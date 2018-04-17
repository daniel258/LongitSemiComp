MargORtoJoint <- function(p1marg, p2marg, OR)
{
  p12 <- ifelse(OR > 0.999 & OR < 1.001, p1marg*p2marg,
     (1 + (p1marg + p2marg)*(OR - 1) - sqrt((1 + (p1marg + p2marg)*(OR - 1))^2 - 4*OR*(OR -1)*p1marg*p2marg))/(2*(OR - 1)))
  p1 <- p1marg - p12
  p2 <- p2marg - p12
  p0 <- 1 - p12 - p1 - p2
  prob <- cbind(p0, p1, p2,p12)
   if (any(prob<0) | any(prob>1)) {
     stop(paste("Problems with probabilities. prob =", prob))  }
  return(prob)
}
# J - number of time points
SimLongitDataParam <- function(n.sample, J,  gamma,  alpha.nt, alpha.t, alpha.or, beta.nt, beta.t, beta.or)
{
 p <- length(beta.nt)
 X <- matrix(nr = n.sample, nc = p, rnorm(n.sample*p))
 risk.NT <- risk.T <- YNT <- YT <- matrix(nr = n.sample, nc = J, 0)
 p.nt <- expit(alpha.nt[1] + X%*%beta.nt)
 p.t <- expit(alpha.t[1] + X%*%beta.t)
 OR <- exp(alpha.or[1]+ X%*%beta.or)
 risk.NT[,1] <- 1
 risk.T[,1] <- 1
 first.probs <- MargORtoJoint(p1marg = p.nt, p2marg = p.t, OR = OR)
 crude.data <- apply(X = first.probs, MARGIN = 1, FUN = sample, x = 1:4, size = 1, replace = T)
 YNT[,1] <- crude.data %in% c(2, 4)
 YT[,1] <- crude.data %in% c(3, 4)
 YNT[YNT[,1]==1, 2] <- 1
 YT[YT[,1]==1, 2] <- 1
 for(j in 2:J)
 {
   risk.NT[,j] <- 1*(YNT[,j-1]==0 & YT[,j-1]==0)
   risk.T[,j] <- 1*(YT[,j-1]==0)
   at.risk.T.only <- risk.NT[,j]==0 & risk.T[,j]==1
  at.risk.both <- risk.NT[,j]==1 & risk.T[,j]==1
  # at risk for terminal event only
  if (sum(at.risk.T.only)>0)
  {
  probs.T.only <- expit(alpha.t[j] + X[at.risk.T.only, ]%*%beta.t + gamma)
  YT[at.risk.T.only, j] <- rbinom(sum(at.risk.T.only), 1, probs.T.only)
  if (j < J)
    {
  YT[YT[,j]==1, j + 1] <- 1
  }  }
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
  if(j < J)  {
  YNT[YNT[, j]==1, j + 1] <- 1
  YT[YT[, j]==1, j + 1] <- 1
  }}}
 return(list(X =X, YNT = YNT, YT = YT, risk.NT = risk.NT, risk.T = risk.T))
}

SimLongitDataSmooth <- function(n.sample, times = 1:100,  gamma,  alpha.nt, alpha.t, alpha.or, beta.nt, beta.t, beta.or)
{
  if (length(alpha.nt) != length(times)) {stop("alpha.nt should be in the same length as times")}
  if (length(alpha.t) != length(times)) {stop("alpha.t should be in the same length as times")}
  if (length(alpha.or) != length(times)) {stop("alpha.or should be in the same length as times")}
  J <- length(times)
  p <- length(beta.nt)
  X <- matrix(nr = n.sample, nc = p, rnorm(n.sample*p))
  risk.NT <- risk.T <- YNT <- YT <- matrix(nr = n.sample, nc = J, 0)
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
      probs.T.only <- expit(alpha.t[j] + X[at.risk.T.only, ]%*%beta.t + gamma)
      YT[at.risk.T.only, j] <- rbinom(sum(at.risk.T.only), 1, probs.T.only)
      # if (j < J)
      # {
        # YT[YT[,j]==1, j + 1] <- 1
      # }
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
      # if(j < J)  {
      #   YNT[YNT[, j]==1, j + 1] <- 1
      #   YT[YT[, j]==1, j + 1] <- 1
      # }
    }
    if(j < J)  {
      YNT[YNT[, j]==1, j + 1] <- 1
      YT[YT[, j]==1, j + 1] <- 1
    }
    }
  return(list(X =X, YNT = YNT, YT = YT, risk.NT = risk.NT, risk.T = risk.T))
}


SimLongitDataSmoothNoX <- function(n.sample, times = 1:100,  gamma,  alpha.nt, alpha.t, alpha.or)
{
  J <- length(times)
  YNT <- YT <- matrix(nr = n.sample, nc = J, 0)
  p.nt.first <- expit(alpha.nt[1])
  p.t.first <- expit(alpha.t[1])
  OR.first <- exp(alpha.or[1])
  first.probs <- MargORtoJoint(p1marg = p.nt.first, p2marg = p.t.first, OR = OR.first)
  crude.data <- apply(X = first.probs, MARGIN = 1, FUN = sample, x = 1:4, size = n.sample, replace = T)
  YNT[,1] <- crude.data %in% c(2, 4)
  YT[,1] <- crude.data %in% c(3, 4)
  YNT[YNT[,1]==1, 2] <- 1
  YT[YT[,1]==1, 2] <- 1
  for(j in 2:J)
  {
    at.risk.T.only <- YNT[,j-1]==1 & YT[,j-1]==0
    at.risk.both <- YNT[,j-1]==0 & YT[,j-1]==0
    # at risk for terminal event only
    if (sum(at.risk.T.only)>0)
    {
      probs.T.only <- expit(alpha.t[j] + gamma)
      YT[at.risk.T.only, j] <- rbinom(sum(at.risk.T.only), 1, probs.T.only)
      if (j < J)
      {
        YT[YT[,j]==1, j + 1] <- 1
      }  }
    #at risk for both events
    if  (sum(at.risk.both)>0)
    {
      p.nt.both <- expit(alpha.nt[j] )
      p.t.both <- expit(alpha.t[j])
      OR.both <- exp(alpha.or[j])
      probs.both <- MargORtoJoint(p1marg = p.nt.both, p2marg = p.t.both, OR = OR.both)
      crude.data.both <- apply(X = probs.both, MARGIN = 1, FUN = sample, x = 1:4, size = sum(at.risk.both), replace = T)
      YNT[at.risk.both, j] <- crude.data.both %in% c(2, 4)
      YT[at.risk.both, j] <- crude.data.both %in% c(3, 4)
      if(j < J)  {
        YNT[YNT[, j]==1, j + 1] <- 1
        YT[YT[, j]==1, j + 1] <- 1
      }}}
  return(list(YNT = YNT, YT = YT))
}
