rm(list=ls())

library(nestedcv)
library(gbm)

# Models from Myles -  Loads etaMod, tocMod and rtxMod objects
load("STRAP_final_models_140524.RData")

var_stability(etaMod)
plot_var_stability(etaMod)
plot_var_stability(tocMod)
plot_var_stability(rtxMod)

var_stability.2 <- function(x,
                            percent = TRUE,
                            level = 1,
                            sort = TRUE, returnFULL = FALSE, ranking = ranking, ...) {

  if(inherits(x, "nestcv.glmnet")){ 
    m <- nestedcv::cv_coef(x, level)
    vdir <- sign(rowMeans(m))
    based_On_freq <- FALSE
  }else if(inherits(x, "nestcv.train")){ 
    message("#  why percentages can be used only for glmnet?")
    message("#  why nestcv.train filters feature based on $freq>0 but glmnet does it using $Final!= 0")
    
    percent <- FALSE 
    based_On_freq <- TRUE  
    m_fv <- x$final_vars
    m <- nestedcv::cv_varImp(x)

    warning("These genes are not in the list of final variables: \n", paste0(sort(rownames(m)[!rownames(m) %in% m_fv], decreasing = F), collapse = "\n"))
    m <- m[rownames(m) %in% m_fv,]
    vdir <- nestedcv::var_direction(x)
    }else{ "unkown class "}
  

  if (percent) {
    # cm <- Rfast::colMaxs(m, value = TRUE)
    cm <- colSums(abs(m))
    if (any(cm == 0)) {
      m <- m[, cm != 0]
      cm <- cm[cm != 0]
    }
    m <- t(t(abs(m)) / cm * 100)
  }

  mm <- rowMeans(m)
  msd <- apply(m, 1, sd)
  msem <- msd/sqrt(ncol(m))
  freq <- apply(m, 1, function(i) sum(i!=0))
  df <- data.frame(mean = mm, sd = msd, sem = msem, frequency = freq, check.names = FALSE)
  
  df$sign <- vdir[rownames(df)]
  df$direction <- factor(df$sign, levels = c(-1, 1),
                         labels = c("negative", "positive"))
  df$final <- "no"
  df$final[m[, "Final"] != 0] <- "yes"
  df$final <- factor(df$final)
  if (!sort){
  df[order(abs(df$mean), decreasing = TRUE), ]
  }

  if(returnFULL){
    mfull <- as.data.frame(m)
    
    if(based_On_freq){
    message("# in github code it is  df <- df[df$freq > 0, ]. It works well but better to use $frequency")
    mfull <- mfull[rownames(df)[df$frequency > 0],] 
    }else{
    mfull <- mfull[mfull$Final!= 0,]
    }

    if(ranking=="reverse_rank"){
      mfull <- apply(mfull, 2, function(x) {rank(-x, ties.method = "average")})
    }else if(ranking=="rank"){
      mfull <- apply(mfull, 2, function(x) {rank(x, ties.method = "average")})
    }else{
      print("No ranking applied")
    }
    
    fac_levels <- rownames(mfull)[order(rowMeans(mfull), decreasing = T)]
    mfull <- reshape::melt(as.matrix(mfull))
    colnames(mfull) <- c("Gene","Fold","VarImpScore")
    mfull$Gene <- factor(mfull$Gene, levels = fac_levels)
    return(mfull)
  }else{
  return(df)
  }
}


var_imp_dotplot <- function(mymodel, ranking=match.arg(c("reverse_rank","rank","no_ranking"))){
library(ggplot2)
library(ggbeeswarm)
library(Hmisc)
set.seed(12345)

mdf <- var_stability.2(x = mymodel, returnFULL = T, percent = T, level = 1, sort = T, ranking = ranking)

p1 <- ggplot(data = mdf, mapping = aes(x = reorder(Gene, rev(Gene)), y = VarImpScore, fill = Gene, col= Gene)) +
    geom_beeswarm(varwidth = T, side = 0, method = "center", dodge.width = 0.5, cex=0.8, 
    width = 0.5,
    data    = mdf,
    shape = 21, stroke = 0.5, size = 1.5
  ) +
  geom_beeswarm(varwidth = T, side = 0, method = "center", dodge.width = 0.5, cex=0.8, 
                width = 0.5,
                data    = mdf,
                shape = 20, stroke = 0.5, size = 0.2, color = "black"
  ) +
  stat_summary(fun = mean, geom = 'point', size = 4, shape = 5, col = "black") +
  scale_y_continuous(breaks = scales::extended_breaks(n = 10)) 
  
if(ranking=="no_ranking" | ranking=="rank"){
    p1 <- p1 + scale_x_discrete(limits = c(rev(levels(mdf$Gene))," ")) + coord_flip()
  }else if(ranking=="reverse_rank"){
    p1 <- p1 + scale_x_discrete(limits = c((levels(mdf$Gene))," ")) + coord_flip()
  }else{
    stop("unkonwn ranking argument")
  }

p1 <- p1 +theme_classic(base_size = 16) + 
  geom_vline(
    linewidth = 0.2,
    data = data.frame(x = seq(from = 1, to = length(levels(mdf$Gene)), by = 1) + 0.5),
    mapping = aes(xintercept = x),  color = "black", linetype = "longdash") + 
  guides(color="none", fill="none") + ylab("Variable importance") + xlab("") +
  theme(axis.ticks = element_line(color = "black"), axis.text = element_text(color = "black"))
  
return(p1)
}


pdf("Figure5d_VarImp_Plot.pdf",  height = 9, width = 6) 
  var_imp_dotplot(mymodel = etaMod, ranking = "rank") + ggtitle("Etanercept") + ylab("Rank")
  var_imp_dotplot(mymodel = tocMod, ranking = "rank") + ggtitle("Tocilizumab") + ylab("Rank")
  var_imp_dotplot(mymodel = rtxMod, ranking = "rank") + ggtitle("Rituximab") + ylab("Rank")
  
  var_imp_dotplot(mymodel = etaMod, ranking = "reverse_rank") + ggtitle("Etanercept") + ylab("Rank")
  var_imp_dotplot(mymodel = tocMod, ranking = "reverse_rank") + ggtitle("Tocilizumab") + ylab("Rank")
  var_imp_dotplot(mymodel = rtxMod, ranking = "reverse_rank") + ggtitle("Rituximab") + ylab("Rank")
dev.off()
