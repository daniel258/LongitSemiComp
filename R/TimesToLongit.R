# Daniel Nevo
#inter.vec - interval used for the grouping of the events into longit
#'@useDynLib LongitSemiComp
#' @export
TimesToLongit <- function(T1,T2, delta1, delta2, inter.vec, TruncData = F, TruncTime = NULL) # TruncTime is the time of entry to the study
{
  J <- length(inter.vec)-1
  YNT <- YT <- matrix(data = 0, nr = length(T1), nc = J)
  risk.T <- risk.NT <- matrix(data = 1, nr = length(T1), nc = J)
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
#' @export
TimesToInter <- function(T1,T2, delta1, delta2, inter.vec)
{
  YNT <- YT <- rep(Inf, length(T1))
  for (i in 1:length(T1))
  {
    if(delta1[i]==1)
    {
      YNT[i] <- findInterval(T1[i], inter.vec, rightmost.closed = T)
    }
    if(delta2[i]==1)
    {
      YT[i] <- findInterval(T2[i], inter.vec, rightmost.closed = T)
    }
  }
  return(list(YNT=YNT, YT=YT))
}
# Rcpp::cppFunction('int FirstOne(NumericVector x) {
#   int index = 0;
#   int length = x.size();
#   for (int i = 0; i < length; i++) {
#             if(x[i]==1) {
#             index = i + 1;
#             break;
#             }}
#             return index;
#             }')
#' @export
LongitToTimes <- function(YNT, YT, times) # should add risk set indicator at some point, right now it assumes no censoring or truncation
{
  n <- nrow(YNT)
  J <- length(times)
  indices.nt <- apply(YNT, 1, FirstOne)
  indices.t <- apply(YT, 1, FirstOne)
  delta1 <- delta2 <- T1 <- T2 <- rep(0, n)
  delta1[indices.nt>0] <- 1
  delta2[indices.t>0] <- 1
  for (i in 1:n)
  {
  if (i %% 10000==0) {cat("i =" , i, "\n")}
  T1[i] <- ifelse(indices.nt[i]>0, times[indices.nt[i]] ,times[J])
  T2[i] <- ifelse(indices.t[i]>0, times[indices.t[i]] ,times[J])
    #if(T1[i] >T2[i] & T1[i]!=times[J]) {
  if(T1[i] >T2[i]) {
    T1[i] <- T2[i]
  #  warning(paste("T1 time was changes to T2 for observation", i, "T1 = ", T1, " T2 = ", T2 ))
    }}
  return(list(T1 = T1, T2 = T2, delta1 = delta1, delta2 = delta2))
}
# inter.vec <- c(0,4,10)
# T1 <- 2
# T2 <- 7
# YNT <- YT <- rep(0, length(inter.vec)-1)
# YNT[inter.vec[-1]>=T1] <- 1
# YT[inter.vec[-1]>=T2] <- 1
# YNT
# YT
