#' @export
FindStatusTimeDep <- function(ID, YNT, YT)
{
  status <- vector(length = length(unique(ID)))
  i = 0
  for (iID in ID)
  {
    i <- i + 1
    iYNT <- YNT[ID==iID]
    iYT <- YT[ID==iID]
    if(any(iYT==1)) { death <- 1
  } else {
    death <- 0
  }
  if(any(iYNT==1)) { AD <- 1
  } else {
    AD <- 0
  }
  if(AD==0 & death==1) {
    status[i] <- "Death without AD"
  } else if  (AD==1 & death==1) {
    status[i] <- "Death with AD"
  } else if (AD==0 & death==0) {
    status[i] <- "Censored without AD or Death"
  } else if (AD==1 & death==0) {
    status[i] <- "Censored with AD"
  }}
  return(table(status))
}