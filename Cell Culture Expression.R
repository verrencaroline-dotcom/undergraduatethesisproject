######################################### Cell Line Expression Analysis 
library(data.table)
library(dplyr)
library(ggplot2)

exprACAP <- fread(
  "D:/Thesis/Downstream Analysis/OmicsExpressionTPMLogp1HumanAllGenes.csv",
  select = c("ModelID", "ACAP2-IT1 (100874306)")
  )
exprRBM <- fread(
  "D:/Thesis/Downstream Analysis/OmicsExpressionTPMLogp1HumanAllGenes.csv",
  select = c("ModelID", "RBM5-AS1 (100775107)")
  )

meta <- fread("D:/Thesis/Downstream Analysis/model.csv")

breast_lines <- meta[OncotreeLineage == "Breast"]

expr_acap_meta <- merge(exprACAP, meta, by = "ModelID")
expracapmetaorder <- expr_acap_meta[order(-expr_acap_meta$`ACAP2-IT1 (100874306)`), ]
expr_rbm_meta <- merge(exprRBM, meta, by = "ModelID")
exprrbmmetaorder <- expr_rbm_meta[order(-expr_rbm_meta$`RBM5-AS1 (100775107)`), ]

ggplot(expr_acap_meta, aes(x = reorder(CellLineName, `ACAP2-IT1 (100874306)`), 
                           +                            y = `ACAP2-IT1 (100874306)`)) +
  geom_bar(stat = "identity", fill = "red") +
  coord_flip() +  # flips axes so cell lines are on y-axis (easier to read)
  labs(title = "ACAP2-IT1 Expression Across 71 Breast Cancer Cell Lines",
             x = "Cell Line",
             y = "Expression (FPKM or TPM)") +
  theme_minimal(base_size = 12) +
  theme(axis.text.y = element_text(size = 7))
ggplot(expr_rbm_meta, aes(x = reorder(CellLineName, `RBM5-AS1 (100775107)`), 
                          y = `RBM5-AS1 (100775107)`)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +  # flips axes so cell lines are on y-axis (easier to read)
  labs(title = "RBM5 AS1 Expression Across 71 Breast Cancer Cell Lines",
             x = "Cell Line",
             y = "Expression (FPKM or TPM)") +
  theme_minimal(base_size = 12) +
  theme(axis.text.y = element_text(size = 7))
