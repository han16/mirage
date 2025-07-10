library(tidyverse)

# Parameters
gene.prior <- 0.1
num.run <- 20
bayesian.fdr <- 0.05
gammamean.values <- 3:6
gamma.bar.values <- 3:6  # bar(gamma) values used within each Gammamean

# Store summary results for each Gammamean
summary_list <- list()

for (gammamean in gammamean.values) {
  actual.fdr.range <- matrix(nrow = length(gamma.bar.values), ncol = num.run)
  
  for (g in seq_along(gamma.bar.values)) {
    gamma.bar <- gamma.bar.values[g]
    
    file.name <- paste0("Mixed_Gene_Gammamean", gammamean, "_useGammabar", gamma.bar, ".delta0.1_replicate100.RData")
    load(file.name)
    
    for (run in 1:num.run) {
      BF <- all.BF.gene[run, ]
      RiskStatus <- all.Ui[run, ]
      post.prob <- gene.prior * BF / (1 - gene.prior + gene.prior * BF)
      gene.post <- tibble(post.prob, RiskStatus)
      
      tau.grid <- seq(0, 0.999, by = 0.001)
      fdr.df <- map_dfr(tau.grid, function(tau) {
        pred <- gene.post$post.prob > tau
        num.pred <- sum(pred)
        false.disc <- sum((1 - gene.post$post.prob)[pred])
        bar.fdr <- ifelse(num.pred > 0, false.disc / num.pred, NA)
        tibble(tau = tau, bar.fdr = bar.fdr)
      }) %>% drop_na()
      
      bayes.fdr.grid <- seq(0, 0.5, length.out = 100)
      actual.fdr <- map_dbl(bayes.fdr.grid, function(fdr.target) {
        tau.sel <- max(fdr.df %>% filter(bar.fdr > fdr.target & bar.fdr < fdr.target + 0.1) %>% pull(tau), na.rm = TRUE)
        selected <- gene.post %>% filter(post.prob > tau.sel)
        TP <- sum(selected$RiskStatus == 1)
        total <- nrow(selected)
        if (total > 0) 1 - TP / total else NA
      })
      
      actual.fdr.range[g, run] <- max(actual.fdr[bayes.fdr.grid < (bayesian.fdr - 1e-4)], na.rm = TRUE)
    }
  }
  
  # Summarize results for this Gammamean
  fdr.mean <- rowMeans(actual.fdr.range, na.rm = TRUE)
  fdr.sd <- apply(actual.fdr.range, 1, sd, na.rm = TRUE)
  summary_df <- tibble(
    Gammamean = gammamean,
    bar.gamma = gamma.bar.values,
    mean = fdr.mean,
    sd = fdr.sd
  )
  
  summary_list[[as.character(gammamean)]] <- summary_df
}

# Combine all summaries
#combined_summary <- bind_rows(summary_list)

# Optional: Plot all in one figure using facet
#ggplot(combined_summary, aes(x = bar.gamma, y = mean)) +
#  geom_point() +
#  geom_pointrange(aes(ymin = mean - sd, ymax = mean + sd)) +
#  geom_hline(yintercept = bayesian.fdr, linetype = "dashed", color = "red") +
#  facet_wrap(~ Gammamean, labeller = label_bquote("True " ~ bar(gamma) == .(Gammamean))) +
#  theme_classic() +
#  xlab(expression(paste(bar(gamma)))) +
#  ylab("False Discovery Proportion") +
#  theme(strip.text = element_text(size = 10), axis.title.y = element_text(size = 8))

