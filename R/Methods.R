#' @export
 print.LongitSC <- function(fit, digits = max(options()$digits - 4, 3),...)
 {
   cat("formula.NT:\n")
   print(fit$formula.NT)
   cat("\n")
   cat("formula.T:\n")
   print(fit$formula.T)
   cat("\n")
   cat("formula.OR:\n")
   print(fit$formula.OR)
   cat("\n")
   cat("\n------------------------------------------------------------------------------")
   cat("\n Non-terminal event :\n")
   coef.NT <- fit$coef.NT
   sd.err.NT <- fit$se.rob.NT
   zvalue.NT <- coef.NT/sd.err.NT
   pvalue.NT <- 2*pnorm(abs(zvalue.NT),lower.tail = F)
   coef.table.NT <- cbind(coef.NT, sd.err.NT, zvalue.NT, pvalue.NT)
   dimnames(coef.table.NT) <- list(names(coef.NT), c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))
   printCoefmat(coef.table.NT, digits = digits)
   cat("\n")
   cat("\n------------------------------------------------------------------------------")
   cat("\n Terminal event :\n")
   coef.T <- fit$coef.T
  coef.T.names <- fit$T.varnames
   sd.err.T <- fit$se.rob.T
   zvalue.T <- coef.T/sd.err.T
   pvalue.T <- 2*pnorm(abs(zvalue.T),lower.tail = F)
   coef.table.T <- cbind(coef.T, sd.err.T, zvalue.T, pvalue.T)
   dimnames(coef.table.T) <- list(names(coef.T), c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))
   printCoefmat(coef.table.T, digits = digits)
    cat("\n------------------------------------------------------------------------------")
    cat("\n Odds ratio between non-terminal and terminal events :\n")
   coef.OR <- fit$coef.OR
   coef.OR.names <- fit$OR.varnames
   sd.err.OR <- fit$se.rob.OR
   zvalue.OR <- coef.OR/sd.err.OR
   pvalue.OR <- 2*pnorm(abs(zvalue.OR),lower.tail = F)
   coef.table.OR <- cbind(coef.OR, sd.err.OR, zvalue.OR, pvalue.OR)
   dimnames(coef.table.OR) <- list(names(coef.OR), c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))
   printCoefmat(coef.table.OR, digits = digits)
 }
 #' @export
 coef.LongitSC <- function(fit, digits = max(options()$digits - 4, 3),...)
 {
   coef.NT <- fit$coef.NT
   coef.T <- fit$coef.T
   coef.OR <- fit$coef.OR
  return(list(coef.NT = coef.NT, coef.T = coef.T, coef.OR = coef.OR))
 }