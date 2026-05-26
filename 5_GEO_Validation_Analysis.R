########################## KM Validation with GEO data 
library(GEOquery)
library(data.table)
library(dplyr)
library(survival)
library(survminer)

## Import GEO expression data
expr_geo <- "./GSE96058/GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz"
expr_geo <- fread(expr_geo, data.table = FALSE)

## Set the first column as gene names
rownames(expr_geo) <- expr_geo$V1

## Remove the first column after assigning it as row names
expr_geo$V1 <- NULL

## Transpose expression matrix
expr_geo_t <- as.data.frame(t(expr_geo))

## Add sample ID column for merging with clinical data
expr_geo_t$geo_accession <- rownames(expr_geo_t)

## Import GEO clinical metadata from the series matrix file
clin_geo <- getGEO(filename = "GSE96058-GPL11154_series_matrix.txt.gz")

## Extract phenotype/clinical data
clin_geo <- pData(clin_geo) 

## Convert overall survival time into numeric format
clin_geo$OS.time <- as.numeric(clin_geo$`overall survival days:ch1`)

## Convert overall survival status into binary format
clin_geo$OS.status <- ifelse(clin_geo$`overall survival event:ch1` == "dead", 1, 0)
clin_geo$geo_accession <- clin_geo$title

## Merge GEO clinical data with expression data
merged_data_geo <- merge(clin_geo, expr_geo_t, by = "geo_accession")

## Create a clean data frame for Kaplan-Meier analysis
geo_Y_df <- merged_data_geo[, c("OS.time", "OS.status", "gene_id_Y")]
colnames(geo_rbm_df) <- c("time", "status", "expr")

## Determine the optimal expression cut-off for lncRNA-Y
cut_Y_geo <- surv_cutpoint(geo_Y_df,
              time = "time",
              event = "status",
              variables = "expr")
## Extract the optimal cut-off value
cut_Y_value <- cut_Y_geo$cutpoint$cutpoint[1]

## Divide samples into high- and low-expression groups based on the optimal cut-off
geo_Y_df$group <- ifelse(geo_Y_df$expr > 0.3076267, "High", "Low")

## Fit Kaplan-Meier survival model
fit_Y_geo <- survfit(Surv(time, status) ~ group, data = geo_Y_df)

## Plot Kaplan-Meier survival curve for lncRNA-Y in GEO dataset
ggsurvplot(fit_Y_geo,
           data = geo_Y_df,
           pval = TRUE,
           risk.table = TRUE,
           title = "lncRNA-Y (GEO Optimal Cutoff)")
