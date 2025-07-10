library(ggplot2)
library(plotROC)
library(gridExtra)

delta_all <- c("0.02","0.05","0.1","0.2")

load(paste0("AUC_summary.rdata"))

selected_rep <- 5
num_run <- 1000
p1_delta <- list()
for (delta in delta_all) {
  
  AUC_selected <- AUC_summary_all[[paste0("delta-",delta)]][selected_rep,]
  
  method_labels <- sprintf("%s(%.3f)", names(AUC_selected), round(AUC_selected, 3))
  method <- rep(method_labels, each = num_run)
  
  load(paste0("7methods_delta",delta, "_gammamean335_5rep.RData"))
  
  Ui <- Gene.Risk.Status[[selected_rep]]
  mirage_BF <- MIRAGE.BF[[selected_rep]]
  fisher_pvalue <- Fisher.pvalue[[selected_rep]]
  skato_pvalue <- SKATO.pvalue[[selected_rep]]
  cmc_pvalue <- CMC.pvalue[[selected_rep]]
  asum_pvalue <- ASUM.pvalue[[selected_rep]]
  acat.pvalue <- ACAT.pvalue[[selected_rep]]
  
  methods_order <- c("ACAT", "ASUM", "Burden", "CMC", "MIRAGE","SKAT-O")
  p_values_list <- list(
    acat.pvalue,
    asum_pvalue,
    fisher_pvalue,
    cmc_pvalue,
    -mirage_BF,  # Note: MIRAGE uses BF values, not p-values
    skato_pvalue
  )
  names(p_values_list) <- methods_order
  
  # Build ROC data frame
  roc_single_run <- data.frame(
    D = Ui,
    m = do.call(c, p_values_list),
    method = method
  )
  
  p1.basic <- ggplot(roc_single_run, aes(d = D, m = m, color = method)) + 
    geom_roc(increasing = FALSE, n.cuts = 0) +
    style_roc(theme = theme_grey)
  
  p1_delta[[paste0("delta-",delta)]] <- p1.basic + 
    ggtitle(bquote(delta == .(delta)))+
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      axis.title.y = element_text(size = 16),
      axis.title.x = element_text(size = 16),
      legend.title = element_blank(),
      legend.justification = c(1, 0),
      legend.position = c(1, 0),
      legend.text = element_text(size = 12)
    ) +
    geom_abline(intercept = 0, slope = 1, col = "black")
  
}

#all <- grid.arrange(p1_delta[["delta-0.02"]], p1_delta[["delta-0.1"]], nrow=1) 
#ggsave(all,filename = paste0("plots/fig2.roc.delta00201_rep",selected_rep,".pdf"),width = 8,height = 4,dpi = 600)
#all <- grid.arrange(p1_delta[["delta-0.05"]], p1_delta[["delta-0.2"]], nrow=1) 
#ggsave(all,filename = paste0("plots/fig2.roc.delta00502_rep",selected_rep,".pdf"),width = 8,height = 4,dpi = 600)
