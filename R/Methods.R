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
   cat("formula.inter:\n")
   print(x$formula.inter)
   cat("\n")
   cat("\n------------------------------------------------------------------------------")
   cat("\n Non-terminal event:\n")
   coef.NT <- x$coef.NT
   exp.coef.NT <- exp(coef.NT)
   if (!is.null(x$se.coef.rob.NT))
   {
   sd.err.NT <- x$se.coef.rob.NT
   } 
    else {
       stop("No robust SE found")#, using standard SE for parameters associated with the non-terminal event")
       #sd.err.NT <- x$se.NT
    }
   zvalue.NT <- coef.NT/sd.err.NT
   pvalue.NT <- 2*pnorm(abs(zvalue.NT),lower.tail = F)
   coef.table.NT <- cbind(coef.NT, sd.err.NT, exp.coef.NT, zvalue.NT, pvalue.NT)
   dimnames(coef.table.NT) <- list(names(coef.NT), c("Estimate", "Std. Error", "Exp(Estimate)", "z value", "Pr(>|z|)"))
   printCoefmat(coef.table.NT, digits = digits)
   cat("\n")
   cat("\n------------------------------------------------------------------------------")
   cat("\n Terminal event :\n")
   coef.T <- x$coef.T
   coef.inter <- x$coef.inter
   exp.coef.T <- exp(coef.T)
   exp.coef.inter <- exp(coef.inter)
   #coef.T.names <- x$T.varnames
   if (!is.null(x$se.coef.rob.T))
   {
      sd.err.T <- x$se.coef.rob.T
   } else {
      stop("No robust SE found")#, using standard SE for parameters associated with the terminal event")
      #sd.err.inter <- x$se.inter
   } 
   if (!is.null(x$se.coef.rob.inter))
   {
      sd.err.inter <- x$se.coef.rob.inter
   } else {
      stop("No robust SE found)#, using standard SE for parameters associated with the terminal event")
      #sd.err.T <- x$se.T
   } 
   zvalue.T <- coef.T/sd.err.T
   zvalue.inter <- coef.inter/sd.err.inter
   pvalue.T <- 2*pnorm(abs(zvalue.T),lower.tail = F)
   pvalue.inter <- 2*pnorm(abs(zvalue.inter),lower.tail = F)
   coef.table.T <- cbind(c(coef.T, coef.inter), c(sd.err.T, sd.err.inter), c(exp.coef.T, exp.coef.inter),  
                                                  c(zvalue.T, zvalue.inter), c(pvalue.T, pvalue.inter))
   dimnames(coef.table.T) <- list(c(names(coef.T), names(coef.inter)), 
                                  c("Estimate", "Std. Error", "Exp(Estimate)", "z value", "Pr(>|z|)"))
   printCoefmat(coef.table.T, digits = digits)
   cat("\n------------------------------------------------------------------------------")
   cat("\n Odds ratio between non-terminal and terminal events :\n")
   coef.OR <- x$coef.OR
   exp.coef.OR <- exp(coef.OR)
   #coef.OR.names <- x$OR.varnames
   if (!is.null(x$se.coef.rob.OR))
   {
      sd.err.OR <- x$se.coef.rob.OR
   } else {
      stop("No robust SE found")#, using standard SE for parameters associated with the odds ratio")
      #sd.err.OR <- x$se.OR
   }   
   zvalue.OR <- coef.OR/sd.err.OR
   pvalue.OR <- 2*pnorm(abs(zvalue.OR),lower.tail = F)
   coef.table.OR <- cbind(coef.OR, sd.err.OR, exp.coef.OR, zvalue.OR, pvalue.OR)
   dimnames(coef.table.OR) <- list(names(coef.OR), c("Estimate", "Std. Error", "Exp(Estimate)", "z value", "Pr(>|z|)"))
   printCoefmat(coef.table.OR, digits = digits)
   cat("\n------------------------------------------------------------------------------")
   cat("\n Long term dependence parameter between non-terminal and terminal events :\n")
   # coef.longterm <- x$coef.longterm
   # exp.coef.longterm <- exp(coef.longterm)
   # sd.err.longterm <- x$se.longterm
   # zvalue.longterm <- coef.longterm/sd.err.longterm
   # pvalue.longterm <- 2*pnorm(abs(zvalue.longterm),lower.tail = F)
   # coef.table.longterm <- cbind(coef.longterm, sd.err.longterm, exp.coef.longterm, zvalue.longterm, pvalue.longterm)
   # dimnames(coef.table.longterm) <- list("Long-term dependence", c("Estimate", "Std. Error", "Exp(Estimate)", "z value", "Pr(>|z|)"))
   # printCoefmat(coef.table.longterm, digits = digits)
 }
 #' @export
 coef.LongitSC <- function(object, ...)
 {
   coef.NT <- object$coef.NT
   coef.T <- object$coef.T
   coef.inter <- object$coef.inter
   coef.OR <- object$coef.OR
  return(list(coef.NT = coef.NT, coef.T = coef.T, coef.inter, coef.OR = coef.OR))
 }