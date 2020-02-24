###### Auxiliary functions #####
############################################################
#### logit and expit####
logit <- function(p) {log(p/(1-p))}
expit <- function(x) {exp(x)/(1+exp(x))}
############################################################
##### Functions that find status at the end of the study ###
############################################################
# Under time-fixed (one row per person).
# This function returns the status fo a single unit
FindStatus <- function(risk.NT, risk.T, YNT, YT)
{
  J <- length(risk.T)
  if(any(YT[!is.na(YT)]==1)) {
    death <- 1
  } else {
    death <- 0
  }
  if(any(YNT[!is.na(YNT)]==1)) 
  {
    non.terminal <- 1
  } else {
    non.terminal <- 0
  }
  if (risk.T[J]==1 & YT[J]==0) {
    CensEnd <- 1
  } else {
    CensEnd <- 0
  }
  if (death==0 & is.na(YT[J])) {
    CensMid <- 1
  } else {
    CensMid <- 0
  }
  if(non.terminal==0 & death==1) {
    status <- "Death without non-terminal event"
  } else if  (non.terminal==1 & death==1) {
    status <- "Death with non-terminal event"
  } else if (non.terminal==0 & CensMid==1) {
    status <- "Censored mid study before non-terminal event or death"
  } else if (non.terminal==1 & CensMid==1) {
    status <- "Censored mid study after non-terminal event before Death"
  } else if (non.terminal==0 & CensEnd==1) {
    status <- "Censored at end of the study without non-terminal event"
  } else if (non.terminal==1 & CensEnd==1) {
    status <- "Censored at end of the study with non-terminal event"
  }
  return(status)
}
#################################################################
# Under multiple rows per-person (similiar to counting process) #
# This function return a status table for the entire sample     #
#################################################################
FindStatusTimeDep <- function(ID, YNT, YT)
{
  status <- vector(length = length(unique(ID)))
  i = 0
  for (iID in unique(ID))
  {
    i <- i + 1
    iYNT <- YNT[ID==iID]
    iYT <- YT[ID==iID]
    if(any(iYT==1)) { death <- 1
    } else {
      death <- 0
    }
    if(any(iYNT==1)) { non.terminal <- 1
    } else {
      non.terminal <- 0
    }
    if(non.terminal==0 & death==1) {
      status[i] <- "Death without non terminal event"
    } else if  (non.terminal==1 & death==1) {
      status[i] <- "Death with non terminal event"
    } else if (non.terminal==0 & death==0) {
      status[i] <- "Censored without non terminal event or death"
    } else if (non.terminal==1 & death==0) {
      status[i] <- "Censored with non terminal event"
    }}
  return(table(status))
}
############################################################################
# Given marginal probabilities and the odds ratio, return cell probabilities
############################################################################
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
