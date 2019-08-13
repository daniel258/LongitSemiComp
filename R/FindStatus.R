#' @export
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
    AD <- 1
  } else {
    AD <- 0
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
  if(AD==0 & death==1) {
    status <- "Death without AD"
  } else if  (AD==1 & death==1) {
    status <- "Death with AD"
  } else if (AD==0 & CensMid==1) {
    status <- "Censored MidStudy before AD or Death"
  } else if (AD==1 & CensMid==1) {
    status <- "Censored MidStudy after AD before Death"
  } else if (AD==0 & CensEnd==1) {
    status <- "Censored at age 100 without AD"
  } else if (AD==1 & CensEnd==1) {
    status <- "Censored at age 100 with AD"
  }
  #c(AD, death, CensEnd, CensMid)
  return(status)
}