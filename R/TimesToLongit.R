#' @title Transfomring semicompeting risks outcome data to longitudinal bivariate binary representation
#' @description Given observed non-terminal and terminal event times, censoring indicators and possible left-truncation, 
#' this function returns the outcome data in the longitudinal bivariate binary representation according to given interval partition
#'  as proposed in Nevo et al.
#' @param T1 Observed non-terminal event time.
#' @param T2 Observed terminal event time.
#' @param delta1 Non-terminal event indicator (1: no-terminal event 0: censored)
#' @param delta2 Terminal event indicator (1: no-terminal event 0: censored)
#' @param inter.vec Partition points creating the interval.
#' @param TruncData Logical. Is the data include left truncation?
#' @param TruncTime Truncation time. Only used if TruncData==T
#' @return A list with at-risk indicators for each unit in each interval (risk.NT and risk.T) and outcome data
#'  at each interval (YNT and YT).
#' @examples
#' \dontrun{
#' data(scrData)
#' plot(scrData$time1, y = scrData$time2)
#' longit.data <- TimesToLongit(T1 = scrData$time1, T2 = scrData$time2, delta1 = scrData$event1, 
#'                              delta2 = scrData$event2, inter.vec = seq(0, 60, 5))
#' longit.data$YNT[1:3,]
#' longit.data$YT[1:3,]
#' longit.data$risk.NT[1:3,]
#' longit.data$risk.T[1:3,]
#' }
# I want the following two lines in the namespace, can be moved to other files if this script is deleted
#'@useDynLib LongitSemiComp 
#'@importFrom Rcpp evalCpp 
#'@import MASS stats
#' @export
TimesToLongit <- function(T1,T2, delta1, delta2, inter.vec, TruncData = F, TruncTime = NULL) # TruncTime is the time of entry to the study
{
  J <- length(inter.vec)-1
  YNT <- YT <- matrix(data = 0, nrow = length(T1), ncol = J)
  risk.T <- risk.NT <- matrix(data = 1, nrow = length(T1), ncol = J)
  if (any(T1 < TruncTime) | any(T2 < TruncTime)) {stop("Data problem: for some observation event time is before truncation time")}
  for (i in 1:length(T1))
  {
    if(delta1[i]==1)
    {
      YNT[i , inter.vec[-1]>=T1[i]] <- 1
    }
    if(delta2[i]==1)
    {
      YT[i, inter.vec[-1]>=T2[i]] <- 1
    }
  indexT1 <- findInterval(T1[i], inter.vec, rightmost.closed = T, left.open = T) + 1
  indexT2 <- findInterval(T2[i], inter.vec, rightmost.closed = T, left.open = T) + 1
  if (indexT1<=J) {risk.NT[i, indexT1:J] <- 0}
  if (indexT2<=J)
      {
      risk.NT[i, indexT2:J] <- 0
      risk.T[i, indexT2:J] <- 0
      }
  if (TruncData == T)
  {
  indexTrunc <- findInterval(TruncTime[i], inter.vec, rightmost.closed = T, left.open = T) - 1
  if (indexTrunc > 1)
      {
    risk.NT[i, 1:indexTrunc] <- 0
    risk.T[i, 1:indexTrunc] <- 0
      }
  }}
  return(list(YNT = YNT, YT = YT, risk.NT = risk.NT, risk.T = risk.T))
}