################################################### DOWNSTREAM ANALYSIS

############## Nomogram Construction
############## Nomogram Construction

library(survival)
library(rms)
library(Hmisc)

# Remove samples with missing stage
surv_data_nomogram <- surv_data[!is.na(surv_data$Stage), ]

# Rename lncRNA expression columns for clearer nomogram labels
surv_data_nomogram$lncRNA_X <- surv_data_nomogram$ENSG00000229325
surv_data_nomogram$lncRNA_Y <- surv_data_nomogram$ENSG00000281691

# Convert Stage to factor
surv_data_nomogram$Stage_factor <- factor(
  surv_data_nomogram$Stage,
  levels = c("I", "IIA", "IIB", "III", "IIIA", "IIIB", "IV")
)

# Optional: set nicer labels for the nomogram
label(surv_data_nomogram$lncRNA_X) <- "lncRNA-X"
label(surv_data_nomogram$lncRNA_Y) <- "lncRNA-Y"
label(surv_data_nomogram$age) <- "Age"
label(surv_data_nomogram$Stage_factor) <- "Stage"

options(contrasts = c("contr.treatment", "contr.treatment"))

# Set datadist after creating/renaming variables
dd <- datadist(surv_data_nomogram)
options(datadist = "dd")

# Fit Cox model using the renamed variables
cox_rms <- cph(
  Surv(OS.time, OS) ~ lncRNA_X + lncRNA_Y + age + Stage_factor,
  data = surv_data_nomogram,
  x = TRUE,
  y = TRUE,
  surv = TRUE
)

# Survival function
surv <- Survival(cox_rms)

# Construct nomogram
nom <- nomogram(
  cox_rms,
  fun = list(
    function(x) surv(365, x),
    function(x) surv(3 * 365, x),
    function(x) surv(5 * 365, x)
  ),
  funlabel = c("1-year OS", "3-year OS", "5-year OS")
)

# Plot nomogram
plot(nom)

############## Functional Enrichment Analysis for both lncRNAs
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

lncRNAs <- c("ENSG00000229325", "ENSG00000281691")

get_coexpressed <- function(exp_mat, lncRNA, cutoff = 0.3){
   if(!lncRNA %in% rownames(exp_mat)){
  stop(paste("lncRNA", lncRNA, "not found in expression matrix"))
  }
   cor_vals <- apply(exp_mat, 1, function(x) cor(x, exp_mat[lncRNA, ], method = "spearman"))
   coexp_genes <- names(cor_vals[cor_vals > cutoff & names(cor_vals) != lncRNA])
   return(coexp_genes)
}

coexp_list <- lapply(lncRNAs, function(l) get_coexpressed(exp_mat = exp_mat, lncRNA = l, cutoff = 0.3))
names(coexp_list) <- lncRNAs
coexp_combined <- unique(unlist(coexp_list))

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_map <- getBM(
  attributes = c("ensembl_gene_id", "entrezgene_id"),
  filters = "ensembl_gene_id",
  values = coexp_combined,
  mart = mart
  )

entrez_ids <- na.omit(gene_map$entrezgene_id)
ego <- enrichGO(
  gene = entrez_ids,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  readable = TRUE
  )
ekegg <- enrichKEGG(
  gene = entrez_ids,
  organism = "hsa",
  pvalueCutoff = 0.05
  )

dotplot(ego, showCategory = 10) + ggtitle("GO Biological Process Enrichment")
dotplot(ekegg, showCategory = 10) + ggtitle("KEGG Pathway Enrichment")

##Functional Enrichment Analysis for individual lncRNAs
FEA_ACAP2_IT1 <- run_enrichment(
  gene_list_ensembl = coexp_list[["ENSG00000229325"]],
  title_name = "lncRNA-X"
  )

FEA_RBM5_AS1 <- run_enrichment(
  gene_list_ensembl = coexp_list[["ENSG00000281691"]],
  title_name = "lncRNA-Y"
  )

########################## Immune Infiltration Analysis

gtf_data <- import("D:/Thesis/gencode.v22.annotation.gtf", format = "gtf")
gtf_df <- as.data.frame(gtf_data)

annot <- unique(gtf_df[, c("gene_id", "gene_name", "gene_type")])
rm(gtf_df)

immune_list <- list(
  B_cells = c("MS4A1","CD79A","CD79B","CD37","CD19","HLA-DRA","HLA-DRB1"),
  CD8_T_cells = c("CD8A","CD8B","GZMB","GZMK","PRF1","NKG7","IFNG"),
  CD4_T_cells = c("CD4","IL7R","CCR7","ICOS"),
  Treg = c("FOXP3","IL2RA","CTLA4","IKZF2","TNFRSF18"),
  NK_cells = c("NKG7","KLRD1","GNLY","PRF1","GZMB"),
  Macrophages = c("CD68","CD163","CSF1R","LYZ","FCGR3A"),
  Dendritic_cells = c("ITGAX","FCER1A","CLEC10A","CD1C"),
  Neutrophils = c("S100A8","S100A9","CXCR2","FCGR3B","MPO")
)
immune_ids_list <- lapply(names(immune_list), function(cell) {
  genes <- immune_list[[cell]]
  immune_annot <- annot[annot$gene_name %in% genes, ]
  immune_ids <- immune_annot$gene_id
  immune_ids <- gsub("\\..*", "", immune_ids)   
  immune_ids
})
names(immune_ids_list) <- names(immune_list)

immune_ids_list <- lapply(immune_ids_list, function(ids) {
  intersect(ids, rownames(exp_mat))
})
sapply(immune_ids_list, length)

immune_list$CD4_T_cells <- c(
  "CD4", "IL7R", "CCR7", "ICOS", "LTB", "MAL", "TRBC1", "TRBC2"
)
immune_list$Dendritic_cells <- c(
  "ITGAX", "FCER1A", "CLEC10A", "CD1C", "BATF3", "IRF7", "LILRA4", "IFI30"
)
immune_ids_list <- lapply(names(immune_list), function(cell) {
  genes <- immune_list[[cell]]
  immune_annot <- annot[annot$gene_name %in% genes, ]
  immune_ids <- immune_annot$gene_id
  immune_ids <- gsub("\\..*", "", immune_ids)
  intersect(immune_ids, rownames(exp_mat))
})

names(immune_ids_list) <- names(immune_list)
sapply(immune_ids_list, length)


library(GSVA)
library(GSEABase)
library(dplyr)

ssgsea_param <- ssgseaParam(expr = exp_mat, geneSets = immune_ids_list)
ssgsea_scores_immune <- gsva(ssgsea_param)

lnc_RBM5  <- exp_mat["ENSG00000281691", ]   
lnc_ACAP2 <- exp_mat["ENSG00000229325", ]  

cor_res <- apply(ssgsea_scores_immune, 1, function(cell) {
  c(
    RBM5_cor = cor(cell, lnc_RBM5, method = "spearman"),
    RBM5_p   = cor.test(cell, lnc_RBM5, method = "spearman")$p.value,
    ACAP2_cor = cor(cell, lnc_ACAP2, method = "spearman"),
    ACAP2_p   = cor.test(cell, lnc_ACAP2, method = "spearman")$p.value
  )
})
cor_res <- as.data.frame(t(cor_res))
cor_res$Immune_Cell <- rownames(cor_res)
cor_res$RBM5_FDR  <- p.adjust(cor_res$RBM5_p, method = "BH")
cor_res$ACAP2_FDR <- p.adjust(cor_res$ACAP2_p, method = "BH")


library(pheatmap)
mat <- cor_res[, c("RBM5_cor", "ACAP2_cor")]
rownames(mat) <- cor_res$Immune_Cell
pheatmap(mat,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         main = "Correlation of lncRNAs with Immune Cells")

library(ggplot2)

group_RBM5 <- ifelse(lnc_RBM5 > median(lnc_RBM5), "High", "Low")
group_ACAP2 <- ifelse(lnc_ACAP2 > median(lnc_ACAP2), "High", "Low")
immune_cells <- rownames(ssgsea_scores_immune)

immune_cells <- rownames(ssgsea_scores_immune)

for(cell in immune_cells){
  # Calculate Spearman correlation
  cor_test <- cor.test(lnc_RBM5, ssgsea_scores_immune[cell, ], method="spearman")
  R <- round(cor_test$estimate, 3)
  p <- signif(cor_test$p.value, 3)
  # Make scatter plot with annotation
  df <- data.frame(
    lncRNA = lnc_RBM5,
    Immune = ssgsea_scores_immune[cell, ]
  )
  
  p1 <- ggplot(df, aes(x=lncRNA, y=Immune)) +
    geom_point(alpha=0.6) +
    geom_smooth(method="lm", se=FALSE, color="blue") +
    theme_minimal() +
    labs(
      title = paste("RBM5-AS1 vs", cell),
      x = "RBM5-AS1 expression",
      y = paste(cell, "ssGSEA score")
    ) 
  
  annotate("text", 
           x = max(df$lncRNA)*0.7, 
           y = max(df$Immune), 
           label = paste0("R=", R, ", p=", p), 
           hjust = 0)
  print(p1)
}

for(cell in immune_cells){
  cor_test <- cor.test(lnc_ACAP2, ssgsea_scores_immune[cell, ], method="spearman")
  R <- round(cor_test$estimate, 3)
  p <- signif(cor_test$p.value, 3)
  df <- data.frame(
    lncRNA = lnc_ACAP2,
    Immune = ssgsea_scores_immune[cell, ]
  )
  
  p2 <- ggplot(df, aes(x=lncRNA, y=Immune)) +
    geom_point(alpha=0.6) +
    geom_smooth(method="lm", se=FALSE, color="red") +
    theme_minimal() +
    labs(
      title = paste("ACAP2-IT1 vs", cell),
      x = "ACAP2-IT1 expression",
      y = paste(cell, "ssGSEA score")
    ) 
  
  annotate("text", 
           x = max(df$lncRNA)*0.7, 
           y = max(df$Immune), 
           label = paste0("R=", R, ", p=", p), 
           hjust = 0)
  print(p2)
}

for(cell in immune_cells){
  df <- data.frame(
    Group = group_RBM5,
    Immune = ssgsea_scores_immune[cell, ]
  )
  
  # Wilcoxon test for High vs Low
  wilcox_res <- wilcox.test(Immune ~ Group, data=df)
  p <- signif(wilcox_res$p.value, 3)
  
  p3 <- ggplot(df, aes(x=Group, y=Immune, fill=Group)) +
    geom_boxplot() +
    theme_minimal() +
    labs(
      title = paste("RBM5-AS1 High vs Low:", cell),
      x = "RBM5-AS1 group",
      y = paste(cell, "ssGSEA score")
    ) +
    
    scale_fill_manual(values=c("lightblue","blue")) +
    annotate("text",
             x=1.5,
             y=max(df$Immune),
             label=paste0("p=", p),
             vjust=-0.5)
  
  print(p3)
}

immune_long <- melt(ssgsea_scores_immune)
colnames(immune_long) <- c("Immune_Cell", "Sample", "Score")

immune_long$RBM5_group <- rep(group_RBM5, each = nrow(ssgsea_scores_immune))
immune_long$ACAP2_group <- rep(group_ACAP2, each = nrow(ssgsea_scores_immune))
ggplot(immune_long, aes(x=Immune_Cell, y=Score, fill=RBM5_group)) +
  geom_boxplot(outlier.size = 0.8) +
  scale_fill_manual(values = c("Low"="#66c2a5", "High"="#fc8d62")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  labs(title="Immune infiltration by RBM5-AS1 expression",
       x="Immune Cell Type",
       y="ssGSEA score",
       fill="RBM5-AS1")

ggplot(immune_long, aes(x=Immune_Cell, y=Score, fill=ACAP2_group)) +
  geom_boxplot(outlier.size = 0.8) +
  scale_fill_manual(values = c("Low"="#8da0cb", "High"="#e78ac3")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  labs(title="Immune infiltration by ACAP2-IT1 expression",
       x="Immune Cell Type",
       y="ssGSEA score",
       fill="ACAP2-IT1")

library(ggpubr)

ggplot(immune_long, aes(x=RBM5_group, y=Score, fill=RBM5_group)) +
  geom_boxplot(outlier.size=0.8) +
  facet_wrap(~Immune_Cell, scales = "free_y") +
  stat_compare_means(method = "wilcox.test") +
  scale_fill_manual(values = c("Low"="#66c2a5", "High"="#fc8d62")) +
  theme_minimal() +
  labs(title="Immune infiltration by RBM5-AS1 (Wilcoxon test)",
       x="RBM5-AS1 group", y="ssGSEA score")

######### Immune Infiltration Separating M1 and M2 

macrophage_gene_sets <- list(
  M1_Macrophages = c("NOS2", "TNF", "IL1B", "IL12A", "IL6", "CXCL10"),
  M2_Macrophages = c("CD163", "MRC1", "IL10", "ARG1", "TGFB1", "CCL22")
)

macrophage_ids_list <- lapply(names(macrophage_gene_sets), function(cell) {
  genes <- macrophage_gene_sets[[cell]]
  macrophage_annot <- annot[annot$gene_name %in% genes, ]
  macrophage_ids<- macrophage_annot$gene_id
  macrophage_ids <- gsub("\\..*", "", macrophage_ids)   
  macrophage_ids
})
names(macrophage_ids_list) <- names(macrophage_gene_sets)
macrophage_ids_list <- lapply(macrophage_ids_list, function(ids) {
  intersect(ids, rownames(exp_mat))
})
sapply(macrophage_ids_list, length)

ssgsea_param <- ssgseaParam(expr = exp_mat, geneSets = macrophage_ids_list)
ssgsea_scores_macro <- gsva(ssgsea_param)

lnc_RBM5  <- exp_mat["ENSG00000281691", ]   
lnc_ACAP2 <- exp_mat["ENSG00000229325", ]  

cor_res <- apply(ssgsea_scores_macro, 1, function(cell) {
  c(
    RBM5_cor = cor(cell, lnc_RBM5, method = "spearman"),
    RBM5_p   = cor.test(cell, lnc_RBM5, method = "spearman")$p.value,
    ACAP2_cor = cor(cell, lnc_ACAP2, method = "spearman"),
    ACAP2_p   = cor.test(cell, lnc_ACAP2, method = "spearman")$p.value
  )
})
cor_res <- as.data.frame(t(cor_res))
cor_res$Immune_Cell <- rownames(cor_res)
cor_res$RBM5_FDR  <- p.adjust(cor_res$RBM5_p, method = "BH")
cor_res$ACAP2_FDR <- p.adjust(cor_res$ACAP2_p, method = "BH")

library(pheatmap)
mat <- cor_res[, c("RBM5_cor", "ACAP2_cor")]
rownames(mat) <- cor_res$Immune_Cell
pheatmap(mat,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         main = "Correlation of lncRNAs with Macrophage Cells")

######### Immune Infiltration Separating N1 and N2 (similar to above just different gene sets)
neutrophil_gene_sets <- list(
  N1_genes <- c("TNF", "ICAM1", "FAS", "CCL3", "CCL5", "CXCL9", "CXCL10", "IL12A", "IL12B", "STAT1"),
  N2_genes <- c("ARG1", "MMP9", "VEGFA", "CXCR4", "IL10", "TGFB1", "S100A8", "S100A9", "PDGFB")
)

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
neutrophil_symbols <- unique(unlist(neutrophil_gene_sets))
mapping <- getBM(
  attributes = c("hgnc_symbol", "ensembl_gene_id"),
  filters = "hgnc_symbol",
  values = neutrophil_symbols,
  mart = mart
)
neutrophil_ids_list <- lapply(neutrophil_gene_sets, function(genes) {
  ids <- mapping$ensembl_gene_id[mapping$hgnc_symbol %in% genes]
  intersect(ids, rownames(exp_mat))
})
sapply(neutrophil_ids_list, length)

ssgsea_param <- ssgseaParam(expr = exp_mat, geneSets = neutrophil_ids_list)
ssgsea_scores_neutro <- gsva(ssgsea_param)

cor_res <- apply(ssgsea_scores_neutro, 1, function(cell) {
  c(
    RBM5_cor = cor(cell, lnc_RBM5, method = "spearman"),
    RBM5_p   = cor.test(cell, lnc_RBM5, method = "spearman")$p.value,
    ACAP2_cor = cor(cell, lnc_ACAP2, method = "spearman"),
    ACAP2_p   = cor.test(cell, lnc_ACAP2, method = "spearman")$p.value
  )
})
cor_res <- as.data.frame(t(cor_res))
cor_res$Immune_Cell <- rownames(cor_res)
cor_res$RBM5_FDR  <- p.adjust(cor_res$RBM5_p, method = "BH")
cor_res$ACAP2_FDR <- p.adjust(cor_res$ACAP2_p, method = "BH")

library(pheatmap)
mat <- cor_res[, c("RBM5_cor", "ACAP2_cor")]
rownames(mat) <- cor_res$Immune_Cell
pheatmap(mat,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         main = "Correlation of lncRNAs with Neutrophil Cells")


############## Drug Sensitivity Analysis -> starts with expression data from cell expression analysis 

library(dplyr)
library(tidyr)
library(readr)
library(tidyverse)
library(readxl)
library(ggplot2)

gdsc_ic50 <- read_excel("D:/Thesis/Downstream Analysis/GDSC2_fitted_dose_response_27Oct23.xlsx")
gdsc_meta <- read_excel("D:/Thesis/Downstream Analysis/Cell_Lines_Details.xlsx")

expr_acap_meta <- expr_acap_meta %>%
  rename(CELL_LINE_NAME = CellLineName)
drug_expr_acap <- merge(expr_acap_meta, gdsc_ic50, by = "CELL_LINE_NAME")
drug_expr_acap <- drug_expr_acap %>%
  rename('ACAP2_IT1' = `ACAP2-IT1 (100874306)`)

cor_results_acap <- drug_expr_acap %>%
  group_by(DRUG_NAME) %>%
  summarise(
    cor_ACAP2_IC50 = cor(ACAP2_IT1, LN_IC50, use = "complete.obs"),
    p_value = cor.test(ACAP2_IT1, LN_IC50)$p.value
    ) %>%
  arrange(p_value)

top_drugs_acap <- cor_results_acap %>% slice_min(p_value, n = 10)

lncRNA_label <- "lncRNA-X"

ggplot(top_drugs_acap, aes(x = reorder(DRUG_NAME, cor_ACAP2_IC50), y = cor_ACAP2_IC50)) +
  geom_bar(stat = "identity", fill = "red") +
  coord_flip() +
  labs(
    title = paste("Top 10 drugs correlated with", lncRNA_label, "expression"),
    x = "Drug",
    y = paste("Pearson correlation (", lncRNA_label, " vs LN_IC50)", sep = "")
  ) +
  theme_minimal()

##### for RBM5 do exact same as above but replace with RBM5-AS1