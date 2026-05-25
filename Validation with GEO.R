########################## KM Validation with GEO data 
library(GEOquery)
library(data.table)

## prepping expression data
expr_geo <- "./GSE96058/GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz"
expr_geo <- fread(expr_geo, data.table = FALSE)

rownames(expr_geo) <- expr_geo$V1
expr_geo$V1 <- NULL
expr_geo_t <- as.data.frame(t(expr_geo))
expr_geo_t$geo_accession <- rownames(expr_geo_t)

## prepping clinical data 
clin_geo <- getGEO(filename = "GSE96058-GPL11154_series_matrix.txt.gz")
clin_geo <- pData(clin_geo) 
clin_geo$OS.time <- as.numeric(clin_geo$`overall survival days:ch1`)
clin_geo$OS.status <- ifelse(clin_geo$`overall survival event:ch1` == "dead", 1, 0)
clin_geo$geo_accession <- clin_geo$title

## full merged data 
merged_data_geo <- merge(clin_geo, expr_geo_t, by = "geo_accession")

##RBM5-AS1 assessment KM 
library(survival)
library(survminer)

geo_rbm_df <- merged_data_geo[, c("OS.time", "OS.status", "RBM5-AS1")]
colnames(geo_rbm_df) <- c("time", "status", "expr")
cut_rbm5_geo <- surv_cutpoint(geo_rbm_df,
              time = "time",
              event = "status",
              variables = "expr")
cut_rbm5_geo

geo_rbm_df$group <- ifelse(geo_rbm_df$expr > 0.3076267, "High", "Low")
fit_rbm5_geo <- survfit(Surv(time, status) ~ group, data = geo_rbm_df)
ggsurvplot(fit_rbm5_geo,
           data = geo_rbm_df,
           pval = TRUE,
           risk.table = TRUE,
           title = "RBM5-AS1 (GEO Optimal Cutoff)")

## ACAP2-IT1 assessment KM -> not possible as it isn't available in GEO  

