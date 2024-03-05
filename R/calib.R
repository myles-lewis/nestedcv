#' Calibration of binary `nestcv.glmnet` fitted object
 
#' Calibration of probabilities by Platt scaling 
#' 
#' @param fit a `nestcv.glmnet` fitted object
#' @importFrom stats predict
#' @return updated `nestcv.glmnet` fitted object with calibrated probabilities 
#'  in the output dataframe of`nestcv.glmnet` fitted object 
#'  under the column of "calib.p"
#' @export
calib <- function(fit){
  logitP <- fit$output$predyp
  #convert logit to p value
  pred <- data.frame("pvalue" = exp(logitP)/(1 + exp(logitP)))
  rownames(pred) <- rownames(fit$output)
  
  pred$bin <- fit$output$testy[match(rownames(pred),
                                     rownames(fit$output))]
  
  calib.model <- stats::glm(bin ~ pvalue, pred, family = binomial)
  
  recalib <- stats::predict(calib.model, pred, type = "response")
  fit$output$calib.p <-   recalib[match(rownames(fit$output),
                                        names(recalib))]
  fit

}


#' Plotting calibrated p values from Platt scaling with original predicted p values
#' 
#' A plot created using ggplot2. 
#' Basis of code built from https://github.com/etlundquist/eRic/
#' 
#' @param fit a `nestcv.glmnet` fitted object
#' @param nbin the number of bins to create the plot
#' @importFrom ggplot2 ggplot geom_line geom_abline coord_cartesian theme
#' scale_color_manual scale_x_continuous scale_y_continuous 
#' @export
calib_plot <- function(fit, nbin){
  output <- fit$output
  output$orig.p <- exp(output$predyp)/(1 + exp(output$predyp))
  #converted logit to p values
  raw.bins <- cut(output$orig.p, nbin, include.lowest = TRUE)
  raw.xval <- tapply(output$orig.p, raw.bins, mean)
  raw.yval <- tapply(as.numeric(as.vector(output$testy)), raw.bins, mean) 
  raw.cali <- data.frame(method = rep('Original', nbin), x = raw.xval, y = raw.yval)
  
  cal.bins <- cut(output$calib.p, nbin, include.lowest = TRUE)
  cal.xval <- tapply(output$calib.p, cal.bins, mean)
  cal.yval <- tapply(as.numeric(as.vector(output$testy)), cal.bins, mean)
  cal.cali <- data.frame(method = rep('Calibrated', nbin), x = cal.xval, y = cal.yval)
  
  comb <- rbind(raw.cali, cal.cali)
  
  ggplot(comb) +
    geom_line(aes(x = x, y = y, color = as.factor(method)), size = 0.75) + 
    geom_abline(slope = 1, intercept = 0, linetype = 2, size = 0.50) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    scale_color_manual(name = '', values = c('red', 'blue')) + 
    scale_x_continuous(breaks = seq(0, 1, 0.25)) + 
    scale_y_continuous(breaks = seq(0, 1, 0.25)) + 
    labs(x = "Predicted probability", 
         y = "Actual class proportion",
         title = "Platt scaling")+
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, size = 14)) +
    theme(axis.title.x = element_text(size = 12)) +
    theme(axis.title.y = element_text(size = 12)) +
    theme(axis.text.x = element_text(size = 12, colour = "black")) +
    theme(axis.text.y = element_text(size = 12, colour = "black")) +
    theme(legend.title = element_text(size = 12)) +
    theme(legend.text = element_text(size = 12))+
    theme(legend.position = c(0.75, 0.1))
}