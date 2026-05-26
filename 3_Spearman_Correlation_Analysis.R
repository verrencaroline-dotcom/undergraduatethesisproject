######################## Spearman Correlation analysis between lncRNA expression and glutathione metabolism ssGSEA score
library(dplyr)
library(ggplot2)
library(ggrepel)


# Import significant lncRNA expression data
lnc_exp_raw <- read.csv(
  "D:/Thesis/DGE/DEG_siglncRNAs.csv",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


# Set gene IDs as row names
rownames(lnc_exp_raw) <- sub("\\..*", "", lnc_exp_raw$gene_id)


# Remove annotation or statistical columns that are not expression values
non_expression_cols <- c(
  "gene_id", "gene_name", "gene_type",
  "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B"
)

lnc_exp <- lnc_exp_raw[, !(colnames(lnc_exp_raw) %in% non_expression_cols)]


# Convert expression data into a numeric matrix
lnc_exp <- as.matrix(lnc_exp)
mode(lnc_exp) <- "numeric"


# Remove genes with missing values or zero variance
lnc_exp <- lnc_exp[
  apply(lnc_exp, 1, function(x) all(!is.na(x)) && var(x) > 0),
]


# Prepare glutathione metabolism ssGSEA score
gsh_score <- as.numeric(ssgsea_scores["Glutathione_Metabolism", ])
names(gsh_score) <- colnames(ssgsea_scores)


# Match samples between the lncRNA expression matrix and glutathione ssGSEA scores
common_samples <- intersect(colnames(lnc_exp), names(gsh_score))

lnc_exp <- lnc_exp[, common_samples]
gsh_score <- gsh_score[common_samples]

# Perform Spearman correlation analysis for each lncRNA
cor_results <- apply(lnc_exp, 1, function(x) {
  res <- cor.test(
    x,
    gsh_score,
    method = "spearman",
    exact = FALSE
  )
  
  c(
    cor = as.numeric(res$estimate),
    p.value = res$p.value
  )
})


# Convert the correlation result matrix into a data frame
cor_df <- as.data.frame(t(cor_results))


# Add gene IDs as a separate column
cor_df$gene_id <- rownames(cor_df)


# Make sure correlation coefficient and p-value columns are numeric
cor_df$cor <- as.numeric(cor_df$cor)
cor_df$p.value <- as.numeric(cor_df$p.value)


# Adjust p-values using the Benjamini-Hochberg method
cor_df$FDR <- p.adjust(cor_df$p.value, method = "BH")


# Filter significantly correlated lncRNAs
sig_cor <- subset(
  cor_df,
  abs(cor) >= 0.25 & FDR < 0.05
)


# Summarize the distribution of absolute correlation coefficients
summary(abs(cor_df$cor))


# Plot the overall distribution of Spearman correlation coefficients
hist(
  cor_df$cor,
  breaks = 50,
  main = "Spearman Correlation Distribution",
  xlab = "Spearman correlation coefficient",
  col = "grey",
  border = "black"
)


# Annotate significantly correlated lncRNAs with gene name and gene type
sig_lnc$gene_id <- sub("\\..*", "", sig_lnc$gene_id)

sig_cor_annot <- merge(
  sig_cor,
  sig_lnc[, c("gene_id", "gene_name", "gene_type")],
  by = "gene_id",
  all.x = TRUE
)


# Export the significantly correlated glutathione-related lncRNAs
write.csv(
  sig_cor_annot,
  file = "./glutathione_related_lncRNAs.csv",
  row.names = FALSE
)


########### Volcano-style plot for glutathione-related lncRNA correlations
# Create plotting data frame
volcano_df <- sig_cor_annot %>%
  mutate(
    cor.rho = cor,
    negLogFDR = -log10(FDR),
    direction = ifelse(cor.rho > 0, "Positive", "Negative")
  )


# If any FDR values are 0, replace them with a very small number
volcano_df$FDR[volcano_df$FDR == 0] <- 1e-300
volcano_df$negLogFDR <- -log10(volcano_df$FDR)


# Select the top 10 most significant lncRNAs based on FDR
top_genes <- volcano_df %>%
  arrange(FDR) %>%
  slice(1:10)


# Create volcano-style correlation plot
ggplot(volcano_df, aes(x = cor.rho, y = negLogFDR)) +
  
  # Plot each lncRNA as one point
  geom_point(aes(color = direction), alpha = 0.8, size = 2) +
  
  # Label the top 10 most significant lncRNAs
  geom_text_repel(
    data = top_genes,
    aes(label = gene_name),
    size = 3,
    max.overlaps = 20
  ) +
  
  # Define colours for positive and negative correlations
  scale_color_manual(
    values = c(
      "Positive" = "#D55E00",
      "Negative" = "#0072B2"
    )
  ) +
  
  # Add a horizontal dashed line showing the FDR = 0.05 threshold
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    color = "grey50"
  ) +
  
  # Add axis labels and title
  labs(
    x = "Spearman correlation coefficient (rho)",
    y = "-log10(FDR)",
    title = "Volcano Plot of Glutathione-Related lncRNAs"
  ) +
  
  # Use a clean classic theme
  theme_classic()