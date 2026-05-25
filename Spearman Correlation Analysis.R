
lnc_exp <- load.csv("D:/Thesis/DGE/DEG_siglncRNAs.csv")


cor_results <- apply(lnc_exp, 1, function(x) {
  +     res <- cor.test(x, gsh_score, method = "spearman", exact = FALSE)
  +     c(cor = res$estimate, p.value = res$p.value)
  + )

cor_df <- as.data.frame(t(cor_results))
cor_df$gene_id <- rownames(cor_df)
cor_df$cor <- as.numeric(cor_df$cor)
cor_df$p.value <- as.numeric(cor_df$p.value)
cor_df$FDR <- p.adjust(cor_df$p.value, method = "BH")

sig_cor <- subset(cor_df, abs(cor) >= 0.25 & FDR < 0.05)

summary(abs(cor_df$cor))
hist(cor_df$cor, breaks = 50, main = "Spearman correlation distribution")

sig_cor_annot <- merge(sig_cor, sig_lnc[, c("gene_id", "gene_name", "gene_type")])
write.csv(write.csv(sig_cor_annot,"./glutathione related lncRNAs.csv"))

## Plotting for all 47 genes
library(ggplot2)
library(dplyr)
library(ggrepel)

top_genes <- volcano_df %>%
  arrange(FDR) %>%
  slice(1:10)

ggplot(volcano_df, aes(x = cor.rho, y = negLogFDR)) +
  geom_point(aes(color = direction), alpha = 0.8, size = 2) +
  geom_text_repel(
    data = top_genes,
    aes(label = gene_name),
    size = 3,
    max.overlaps = 20
  ) +
  scale_color_manual(values = c("Positive" = "#D55E00",
                                "Negative" = "#0072B2")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  labs(
    x = "Spearman correlation coefficient (ρ)",
    y = "-log10(FDR)",
    title = "Volcano plot of glutathione-related lncRNAs"
  ) +
  theme_classic()


