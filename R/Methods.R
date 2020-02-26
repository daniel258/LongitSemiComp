#' @export
 print.LongitSC <- function(x, digits = max(options()$digits - 4, 3),...)
 {
   cat("formula.NT:\n")
   print(x$formula.NT)
   cat("\n")
   cat("formula.T:\n")
   print(x$formula.T)
   cat("\n")
   cat("formula.OR:\n")
   print(x$formula.OR)
   cat("\n")
   cat("\n------------------------------------------------------------------------------")
   cat("\n Non-terminal event :\n")
   coef.NT <- x$coef.NT
   sd.err.NT <- x$se.rob.NT
   zvalue.NT <- coef.NT/sd.err.NT
   pvalue.NT <- 2*pnorm(abs(zvalue.NT),lower.tail = F)
   coef.table.NT <- cbind(coef.NT, sd.err.NT, zvalue.NT, pvalue.NT)
   dimnames(coef.table.NT) <- list(names(coef.NT), c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))
   printCoefmat(coef.table.NT, digits = digits)
   cat("\n")
   cat("\n------------------------------------------------------------------------------")
   cat("\n Terminal event :\n")
   coef.T <- x$coef.T
  coef.T.names <- x$T.varnames
   sd.err.T <- x$se.rob.T
   zvalue.T <- coef.T/sd.err.T
   pvalue.T <- 2*pnorm(abs(zvalue.T),lower.tail = F)
   coef.table.T <- cbind(coef.T, sd.err.T, zvalue.T, pvalue.T)
   dimnames(coef.table.T) <- list(names(coef.T), c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))
   printCoefmat(coef.table.T, digits = digits)
    cat("\n------------------------------------------------------------------------------")
    cat("\n Odds ratio between non-terminal and terminal events :\n")
   coef.OR <- x$coef.OR
   coef.OR.names <- x$OR.varnames
   sd.err.OR <- x$se.rob.OR
   zvalue.OR <- coef.OR/sd.err.OR
   pvalue.OR <- 2*pnorm(abs(zvalue.OR),lower.tail = F)
   coef.table.OR <- cbind(coef.OR, sd.err.OR, zvalue.OR, pvalue.OR)
   dimnames(coef.table.OR) <- list(names(coef.OR), c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))
   printCoefmat(coef.table.OR, digits = digits)
 }
 #' @export
 coef.LongitSC <- function(object, ...)
 {
   coef.NT <- object$coef.NT
   coef.T <- object$coef.T
   coef.OR <- object$coef.OR
  return(list(coef.NT = coef.NT, coef.T = coef.T, coef.OR = coef.OR))
 }