################################################### DOWNSTREAM ANALYSIS

## The downstream analyses include:
# 1. Nomogram construction
# 2. Functional enrichment analysis
# 3. Immune infiltration analysis
# 4. Macrophage subtype analysis
# 5. Neutrophil subtype analysis
# 6. Drug sensitivity analysis

library(survival)
library(survminer)
library(rms)
library(Hmisc)
library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(rtracklayer)
library(GSVA)
library(GSEABase)
library(pheatmap)
library(ggpubr)
library(data.table)
library(readxl)
library(tidyr)

############## Nomogram Construction

## Remove samples with missing pathological stage information
surv_data_nomogram <- surv_data[!is.na(surv_data$Stage), ]

## Rename lncRNA expression columns for clearer nomogram labels
surv_data_nomogram$lncRNA_X <- surv_data_nomogram[[gene_id_X]]
surv_data_nomogram$lncRNA_Y <- surv_data_nomogram[[gene_id_Y]]

## Convert pathological stage into a factor
surv_data_nomogram$Stage_factor <- factor(
  surv_data_nomogram$Stage,
  levels = c("I", "IIA", "IIB", "III", "IIIA", "IIIB", "IV")
)

## Add labels to variables so the nomogram is easier to interpret
label(surv_data_nomogram$lncRNA_X) <- "lncRNA-X"
label(surv_data_nomogram$lncRNA_Y) <- "lncRNA-Y"
label(surv_data_nomogram$age) <- "Age"
label(surv_data_nomogram$Stage_factor) <- "Stage"

## Set contrast options for categorical variables
options(contrasts = c("contr.treatment", "contr.treatment"))

## Set datadist after creating/renaming variables
dd <- datadist(surv_data_nomogram)
options(datadist = "dd")

## Fit Cox proportional hazards model using rms::cph()
cox_rms <- cph(
  Surv(OS.time, OS) ~ lncRNA_X + lncRNA_Y + age + Stage_factor,
  data = surv_data_nomogram,
  x = TRUE,
  y = TRUE,
  surv = TRUE
)

## Create survival function from the fitted Cox model
surv <- Survival(cox_rms)

## Construct nomogram for 1-year, 3-year, and 5-year overall survival
nom <- nomogram(
  cox_rms,
  fun = list(
    function(x) surv(365, x),
    function(x) surv(3 * 365, x),
    function(x) surv(5 * 365, x)
  ),
  funlabel = c("1-year OS", "3-year OS", "5-year OS")
)
plot(nom)

############## Functional Enrichment Analysis 
## Define candidate lncRNAs for co-expression analysis
lncRNAs <- c("lncRNA-X" = gene_id_X,
             "lncRNA-Y" = gene_id_Y)

## Function to identify genes co-expressed with a selected lncRNA
get_coexpressed <- function(exp_mat, lncRNA, cutoff = 0.3){
  
## Check whether the selected lncRNA exists in the expression matrix
   if(!lncRNA %in% rownames(exp_mat)){
  stop(paste("lncRNA", lncRNA, "not found in expression matrix"))
   }
  
## Calculate Spearman correlation between the selected lncRNA
   cor_vals <- apply(exp_mat, 1, function(x) cor(x, exp_mat[lncRNA, ], method = "spearman"))
  
## Select positively co-expressed genes above the correlation cut-off
   coexp_genes <- names(cor_vals[cor_vals > cutoff & names(cor_vals) != lncRNA])
   return(coexp_genes)
}

## Identify co-expressed genes for each lncRNA
coexp_list <- lapply(lncRNAs, function(l) get_coexpressed(exp_mat = exp_mat, lncRNA = l, cutoff = 0.3))

## Name the co-expression lists using thesis-facing labels
names(coexp_list) <- lncRNAs

## Combine co-expressed genes from lncRNA-X and lncRNA-Y
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
## Import GTF annotation file
gtf_data <- import("D:/Thesis/gencode.v22.annotation.gtf", format = "gtf")
gtf_df <- as.data.frame(gtf_data)

## Create a simplified annotation table
annot <- unique(gtf_df[, c("gene_id", "gene_name", "gene_type")])
rm(gtf_df)

## Define immune cell marker gene sets
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

## Function to convert marker gene symbols into Ensembl IDs
immune_ids_list <- lapply(names(immune_list), function(cell) {
  genes <- immune_list[[cell]]
  immune_annot <- annot[annot$gene_name %in% genes, ]
  immune_ids <- immune_annot$gene_id
  immune_ids <- gsub("\\..*", "", immune_ids)   
  immune_ids
})
names(immune_ids_list) <- names(immune_list)

## Convert immune marker genes into Ensembl IDs
immune_ids_list <- lapply(immune_ids_list, function(ids) {
  intersect(ids, rownames(exp_mat))
})

## Check how many marker genes were matched for each immune cell type
sapply(immune_ids_list, length)

## Run ssGSEA to estimate immune infiltration scores
ssgsea_param <- ssgseaParam(expr = exp_mat, geneSets = immune_ids_list)
ssgsea_scores_immune <- gsva(ssgsea_param)

## Extract lncRNA-X and lncRNA-Y expression vectors
lnc_X  <- exp_mat[gene_id_X, ]   
lnc_Y <- exp_mat[gene_id_Y, ]  

## Function to correlate lncRNA expression with immune ssGSEA scores
cor_res <- apply(ssgsea_scores_immune, 1, function(cell) {
  c(
    lncX_cor = cor(cell, lnc_X, method = "spearman"),
    lncX_p   = cor.test(cell, lnc_X, method = "spearman")$p.value,
    lncY_cor = cor(cell, lnc_Y, method = "spearman"),
    lncY_p   = cor.test(cell, lnc_Y, method = "spearman")$p.value
  )
})
cor_res <- as.data.frame(t(cor_res))
cor_res$Immune_Cell <- rownames(cor_res)

## Adjust p-values using Benjamini-Hochberg correction
cor_res$lncX_FDR  <- p.adjust(cor_res$lncX_p, method = "BH")
cor_res$lncY_FDR <- p.adjust(cor_res$lncY_p, method = "BH")

## Create correlation matrix for heatmap
mat <- cor_res[, c("lncX_cor", "lncY_cor")]
rownames(mat) <- cor_res$Immune_Cell

## Plot heatmap of lncRNA-immune correlations
pheatmap(mat,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         main = "Correlation of lncRNAs with Immune Cells")

## Create high- and low-expression groups based on median lncRNA expression
group_lncX <- ifelse(lnc_X > median(lnc_X), "High", "Low")
group_lncY <- ifelse(lnc_Y > median(lnc_Y), "High", "Low")

## Extract immune cell names
immune_cells <- rownames(ssgsea_scores_immune)

## Scatter plots for lncRNA-X against immune cell ssGSEA scores
for(cell in immune_cells){
  cor_test <- cor.test(lnc_x, ssgsea_scores_immune[cell, ], method="spearman")
  R <- round(cor_test$estimate, 3)
  p <- signif(cor_test$p.value, 3)
  df <- data.frame(
    lncRNA = lnc_x,
    Immune = ssgsea_scores_immune[cell, ]
  )
  
  p1 <- ggplot(df, aes(x=lncRNA, y=Immune)) +
    geom_point(alpha=0.6) +
    geom_smooth(method="lm", se=FALSE, color="blue") +
    theme_minimal() +
    labs(
      title = paste("lncRNA-X vs", cell),
      x = "lncRNA-X expression",
      y = paste(cell, "ssGSEA score")
    ) 
  
  annotate("text", 
           x = max(df$lncRNA)*0.7, 
           y = max(df$Immune), 
           label = paste0("R=", R, ", p=", p), 
           hjust = 0)
  print(p1)
}

## Scatter plots for lncRNA-Y against immune cell ssGSEA scores
for(cell in immune_cells){
  cor_test <- cor.test(lnc_Y, ssgsea_scores_immune[cell, ], method="spearman")
  R <- round(cor_test$estimate, 3)
  p <- signif(cor_test$p.value, 3)
  df <- data.frame(
    lncRNA = lnc_Y,
    Immune = ssgsea_scores_immune[cell, ]
  )
  
  p2 <- ggplot(df, aes(x=lncRNA, y=Immune)) +
    geom_point(alpha=0.6) +
    geom_smooth(method="lm", se=FALSE, color="red") +
    theme_minimal() +
    labs(
      title = paste("lncRNA-Y vs", cell),
      x = "lncRNA-Y expression",
      y = paste(cell, "ssGSEA score")
    ) 
  
  annotate("text", 
           x = max(df$lncRNA)*0.7, 
           y = max(df$Immune), 
           label = paste0("R=", R, ", p=", p), 
           hjust = 0)
  print(p2)
}
######### 5. MACROPHAGE SUBTYPE ANALYSIS: M1 AND M2

## Define M1 and M2 macrophage marker gene sets
macrophage_gene_sets <- list(
  M1_Macrophages = c("NOS2", "TNF", "IL1B", "IL12A", "IL6", "CXCL10"),
  M2_Macrophages = c("CD163", "MRC1", "IL10", "ARG1", "TGFB1", "CCL22")
)

## Convert macrophage marker genes into Ensembl IDs
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

## Calculate ssGSEA scores for macrophage subtypes
ssgsea_param <- ssgseaParam(expr = exp_mat, geneSets = macrophage_ids_list)
ssgsea_scores_macro <- gsva(ssgsea_param)

## Correlate lncRNAs with macrophage subtype scores
cor_res <- apply(ssgsea_scores_macro, 1, function(cell) {
  c(
    LncX_cor = cor(cell, lnc_X, method = "spearman"),
    LncX_p   = cor.test(cell, lnc_X, method = "spearman")$p.value,
    LncY_cor = cor(cell, lnc_Y, method = "spearman"),
    LncY_p   = cor.test(cell, lnc_Y, method = "spearman")$p.value
  )
})

## Create macrophage correlation heatmap
mat <- cor_res[, c("LncX_cor", "LncY_cor")]
rownames(mat) <- cor_res$Immune_Cell
pheatmap(mat,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         main = "Correlation of lncRNAs with Macrophage Cells")

######### Immune Infiltration Separating N1 and N2 (similar to above just different gene sets)
## Define N1 and N2 neutrophil marker gene sets
neutrophil_gene_sets <- list(
  N1_genes <- c("TNF", "ICAM1", "FAS", "CCL3", "CCL5", "CXCL9", "CXCL10", "IL12A", "IL12B", "STAT1"),
  N2_genes <- c("ARG1", "MMP9", "VEGFA", "CXCR4", "IL10", "TGFB1", "S100A8", "S100A9", "PDGFB")
)

## Convert neutrophil marker genes into Ensembl IDs
neutrophil_ids_list <- lapply(neutrophil_gene_sets, function(genes) {
  ids <- mapping$ensembl_gene_id[mapping$hgnc_symbol %in% genes]
  intersect(ids, rownames(exp_mat))
})
sapply(neutrophil_ids_list, length)

## Calculate ssGSEA scores for neutrophil subtypes
ssgsea_param <- ssgseaParam(expr = exp_mat, geneSets = neutrophil_ids_list)
ssgsea_scores_neutro <- gsva(ssgsea_param)

## Correlate lncRNAs with neutrophil subtype scores
cor_res <- apply(ssgsea_scores_neutro, 1, function(cell) {
  c(
    LncX_cor = cor(cell, lnc_X, method = "spearman"),
    LncX_p   = cor.test(cell, lnc_X, method = "spearman")$p.value,
    LncY_cor = cor(cell, lnc_Y, method = "spearman"),
    LncY_p   = cor.test(cell, lnc_Y, method = "spearman")$p.value
  )
})

## Create neutrophil correlation heatmap
mat <- cor_res[, c("LncX_cor", "LncY_cor")]
rownames(mat) <- cor_res$Immune_Cell
pheatmap(mat,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         main = "Correlation of lncRNAs with Neutrophil Cells")


############## Drug Sensitivity Analysis 

## Import GDSC IC50 data
gdsc_ic50 <- read_excel("D:/Thesis/Downstream Analysis/GDSC2_fitted_dose_response_27Oct23.xlsx")

## Import cell line metadata
gdsc_meta <- read_excel("D:/Thesis/Downstream Analysis/Cell_Lines_Details.xlsx")

## Prepare lncRNA-X cell line expression data
expr_X_meta <- expr_X_meta %>%
  rename(CELL_LINE_NAME = CellLineName)

## Merge lncRNA-X expression with GDSC drug sensitivity data
drug_expr_X <- merge(expr_acap_meta, gdsc_ic50, by = "CELL_LINE_NAME")
drug_expr_X <- drug_expr_X %>%
  rename('lncRNA-X' = `gene_id_X`)

## Function to perform drug sensitivity correlation analysis
cor_results_lncRNAX <- drug_expr_X %>%
  group_by(DRUG_NAME) %>%
  summarise(
    cor_ACAP2_IC50 = cor(lncRNA_X, LN_IC50, use = "complete.obs"),
    p_value = cor.test(lncRNA_X, LN_IC50)$p.value
    ) %>%
  arrange(p_value)

## Select top 10 drugs most significantly correlated with lncRNA-X
top_drugs_acap <- cor_results_lncRNAX %>% slice_min(p_value, n = 10)

## Plot top drug correlations for lncRNA-X
ggplot(top_drugs_acap, aes(x = reorder(DRUG_NAME, cor_ACAP2_IC50), y = cor_ACAP2_IC50)) +
  geom_bar(stat = "identity", fill = "red") +
  coord_flip() +
  labs(
    title = paste("Top 10 drugs correlated with lncRNAX", lncRNA_label, "expression"),
    x = "Drug",
    y = paste("Pearson correlation (", lncRNA_label, " vs LN_IC50)", sep = "")
  ) +
  theme_minimal()

##### for lncRNA-Y do exact same as above but replace with lncRNA-Y