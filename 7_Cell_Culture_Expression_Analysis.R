######################################### Cell Line Expression Analysis 
library(data.table)
library(dplyr)
library(ggplot2)

# Import lncRNA-X expression data from the DepMap expression file
expr_X <- fread(
  "D:/Thesis/Downstream Analysis/OmicsExpressionTPMLogp1HumanAllGenes.csv",
  select = c("ModelID", "gene id")
)


# Import lncRNA-Y expression data from the DepMap expression file
expr_Y <- fread(
  "D:/Thesis/Downstream Analysis/OmicsExpressionTPMLogp1HumanAllGenes.csv",
  select = c("ModelID", "gene id")
)

# Import cell line metadata
meta <- fread("D:/Thesis/Downstream Analysis/model.csv")


# Filter metadata to include only breast cancer cell lines
breast_lines <- meta %>%
  filter(OncotreeLineage == "Breast")


# Merge lncRNA-X expression data with breast cancer cell line metadata
expr_X_meta <- merge(
  expr_X,
  breast_lines,
  by = "ModelID"
)


# Rank breast cancer cell lines from highest to lowest lncRNA-X expression
expr_X_meta_order <- expr_X_meta %>%
  arrange(desc(lncRNA_X))


# Merge lncRNA-Y expression data with breast cancer cell line metadata
expr_Y_meta <- merge(
  expr_Y,
  breast_lines,
  by = "ModelID"
)


# Rank breast cancer cell lines from highest to lowest lncRNA-Y expression
expr_Y_meta_order <- expr_Y_meta %>%
  arrange(desc(lncRNA_Y))


# Plot lncRNA-X expression across breast cancer cell lines
ggplot(
  expr_X_meta,
  aes(
    x = reorder(CellLineName, lncRNA_X),
    y = lncRNA_X
  )
) +
  geom_bar(stat = "identity", fill = "red") +
  
  # Flip axes so that cell line names are easier to read
  coord_flip() +
  
  # Add plot title and axis labels
  labs(
    title = "lncRNA-X Expression Across Breast Cancer Cell Lines",
    x = "Cell Line",
    y = "Expression, log2(TPM + 1)"
  ) +
  
  # Use a clean minimal theme
  theme_minimal(base_size = 12) +
  
  # Reduce y-axis text size because many cell lines are shown
  theme(
    axis.text.y = element_text(size = 7)
  )


# Plot lncRNA-Y expression across breast cancer cell lines
ggplot(
  expr_Y_meta,
  aes(
    x = reorder(CellLineName, lncRNA_Y),
    y = lncRNA_Y
  )
) +
  geom_bar(stat = "identity", fill = "steelblue") +
  
  # Flip axes so that cell line names are easier to read
  coord_flip() +
  
  # Add plot title and axis labels
  labs(
    title = "lncRNA-Y Expression Across Breast Cancer Cell Lines",
    x = "Cell Line",
    y = "Expression, log2(TPM + 1)"
  ) +
  
  # Use a clean minimal theme
  theme_minimal(base_size = 12) +
  
  # Reduce y-axis text size because many cell lines are shown
  theme(
    axis.text.y = element_text(size = 7)
  )
