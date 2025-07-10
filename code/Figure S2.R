
#library(dyplr)
library(ggplot2)
library(gridExtra)
library(matrixStats)
library(plotROC)
library(grid)
library(tidyverse)
library(cowplot)


simulated_data <- readRDS("simulated_data_random_num_variants.rds")
delta=c(0.02, 0.05, 0.1, 0.2)
rep_use_index=c(3,2,1,3)
p_delta=list()
p_delta_burden=list()

AUC=list(
  "0.02"=c("Burden(0.812)", "SKAT-O(0.876)", "MIRAGE(0.976)", "CMC(0.899)", "ASUM(0.762)", "ACAT(0.653)"),
  "0.05"=c( "Burden(0.888)", "SKAT-O(0.918)", "MIRAGE(0.984)", "CMC(0.899)", "ASUM(0.806)", "ACAT(0.756)"),
  "0.1"=c("Burden(0.811)", "SKAT-O(0.878)", "MIRAGE(0.961)", "CMC(0.879)", "ASUM(0.730)", "ACAT(0.703)"),
  "0.2"=c( "Burden(0.858)", "SKAT-O(0.915)", "MIRAGE(0.974)", "CMC(0.909)", "ASUM(0.825)", "ACAT(0.745)"))


AUC_burden=list(
  "0.02"=c("Burden(0.812)", "Burden-adj(0.839)", "Burden-combine(0.837)"),
  "0.05"=c("Burden(0.888)", "Burden-adj(0.910)", "Burden-combine(0.914)"),
  "0.1"=c("Burden(0.811)", "Burden-adj(0.843)", "Burden-combine(0.836)"),
  "0.2"=c( "Burden(0.858)", "Burden-adj(0.920)", "Burden-combine(0.923)"))

# Define reusable theme
plot_roc <- function(title_size = 14, axis_y_size = 14, legend_text_size = 11) {
  list(
    theme_classic(),
    theme(
      plot.title = element_text(hjust = 0.5, size = title_size),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = axis_y_size),
      axis.title.y = element_text(size = 16),
      axis.title.x = element_text(size = 16),
      legend.title = element_blank(),
      legend.justification = c(1, 0),
      legend.position = c(1, 0),
      legend.text = element_text(size = legend_text_size)
    ),
    geom_abline(intercept = 0, slope = 1, col = "black")
  )
}

# Loop over delta values
for (i in seq_along(rep_use_index)) {
  
 
  Ui <- simulated_data[[i]]$Ui
  mirage_pvalue <- simulated_data[[i]]$mirage_pvalue
  mirage_BF <- simulated_data[[i]]$mirage_BF
  fisher_pvalue <- simulated_data[[i]]$fisher_pvalue
  skato_pvalue <- simulated_data[[i]]$skato_pvalue
  cmc_pvalue <- simulated_data[[i]]$cmc_pvalue
  asum_pvalue <- simulated_data[[i]]$asum_pvalue
  acat.pvalue <- simulated_data[[i]]$acat.pvalue
  fisher.separate.pvalue <- simulated_data[[i]]$fisher.separate.pvalue
  
  
  # Calculate adjusted and combined Fisher p-values
  fisher.adj.pvalue <- apply(fisher.separate.pvalue, 1, function(p) min(p.adjust(p, method = "bonferroni")))
  fisher.combine.pvalue <- apply(fisher.separate.pvalue, 1, function(p) pchisq(-sum(log(p)), df = length(p), lower.tail = TRUE))
  
  num_run <- length(mirage_pvalue)
  
  ### ROC for 6 methods
  method <- rep(AUC[[as.character(delta[i])]], each = num_run)
  roc_single_run <- data.frame(
    D = Ui,
    m = c(fisher_pvalue, skato_pvalue, -mirage_BF, cmc_pvalue, asum_pvalue, acat.pvalue),
    method = method
  )
  
  p1.basic <- ggplot(roc_single_run, aes(d = D, m = m, color = method), size = 2) +
    geom_roc(increasing = FALSE, n.cuts = 0) +
    style_roc(theme = theme_grey)
  
  p_delta[[i]] <- p1.basic +
    ggtitle(bquote(delta == .(delta[i]))) +
    plot_roc()
  
  ### ROC for 3 burden methods
  method_burden <- rep(AUC_burden[[as.character(delta[i])]], each = num_run)
  roc_single_run_burden <- data.frame(
    D = Ui,
    m = c(fisher_pvalue, fisher.adj.pvalue, -fisher.combine.pvalue),
    method = method_burden
  )
  
  p.basic_burden <- ggplot(roc_single_run_burden, aes(d = D, m = m, color = method), size = 2) +
    geom_roc(increasing = FALSE, n.cuts = 0, labels = FALSE) +
    style_roc(theme = theme_grey)
  
  p_delta_burden[[i]] <- p.basic_burden +
    ggtitle(bquote(delta == .(delta[i]))) +
    plot_roc(axis_y_size = 14, legend_text_size = 13)
}

names(p_delta)=paste0("delta_", delta)
names(p_delta_burden)=paste0("delta_", delta)



p_delta
p_delta_burden
