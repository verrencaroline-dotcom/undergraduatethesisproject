###################### Prognostic Model Construction and Validation

# The workflow includes:
# 1. Preparing clinical and survival data
# 2. Merging clinical data with lncRNA expression data
# 3. Performing univariate Cox regression
# 4. Performing multivariate Cox regression
# 5. Generating Kaplan-Meier survival curves
# 6. Constructing a combined two-lncRNA risk score
# 7. Evaluating predictive performance using time-dependent ROC curves
# 8. Visualizing survival status across risk scores


##########################1. Preparing clinical and survival data
library(survminer)
library(survival)
library(dplyr)
library(stringr)
library(ggplot2)
library(timeROC)


# Import clinical data
clinical <- read.delim(
  "clinical.tsv",
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# Select only the clinical covariates needed for prognostic analysis
clinical_covariates <- clinical[, c(
  "sample",
  "age_at_earliest_diagnosis_in_years.diagnoses.xena_derived",
  "ajcc_pathologic_t.diagnoses",
  "ajcc_pathologic_n.diagnoses",
  "ajcc_pathologic_m.diagnoses",
  "ajcc_pathologic_stage.diagnoses"
)]


# Rename columns into simpler names for easier downstream analysis
colnames(clinical_covariates) <- c(
  "sample", "age", "T_stage", "N_stage", "M_stage", "Stage"
)


# Remove the final "A" from sample names if present
clinical_covariates$sample <- gsub("A$", "", clinical_covariates$sample)


# Merge survival data with selected clinical covariates
clinical_full <- merge(
  survival_brac,
  clinical_covariates,
  by = "sample"
)


# Remove survival outcomes that are not used in this analysis
clinical_full <- clinical_full %>%
  select(-any_of(c(
    "DSS", "DSS.time",
    "DFI", "DFI.time",
    "PFI", "PFI.time"
  )))


# Keep only samples with complete overall survival, stage, and age information
clinical_full <- clinical_full %>%
  filter(
    !is.na(OS),
    !is.na(OS.time),
    !is.na(Stage),
    !is.na(age)
  )


# Clean and standardize pathological stage values
clinical_full$Stage <- toupper(clinical_full$Stage)
clinical_full$Stage <- gsub("STAGE ", "", clinical_full$Stage)
clinical_full$Stage <- gsub("[^A-Z0-9]", "", clinical_full$Stage)


# Convert stage into an ordered factor
clinical_full$Stage <- factor(
  clinical_full$Stage,
  levels = c("I", "IA", "IB",
             "II", "IIA", "IIB",
             "III", "IIIA", "IIIB", "IIIC",
             "IV"),
  ordered = TRUE
)


# Convert age, survival time, and survival status into numeric variables
clinical_full$age <- as.numeric(clinical_full$age)
clinical_full$OS.time <- as.numeric(clinical_full$OS.time)
clinical_full$OS <- as.numeric(clinical_full$OS)


# Remove samples with missing values introduced during numeric conversion
clinical_full <- clinical_full %>%
  filter(
    !is.na(age),
    !is.na(OS.time),
    !is.na(OS)
  )


# Remove duplicated samples if any are present
clinical_full <- clinical_full[!duplicated(clinical_full$sample), ]


################################# 2. Subset significant lncRNAs and merge with survival data

# Remove version numbers from gene IDs in the expression matrix
rownames(exp) <- sub("\\..*", "", rownames(exp))


# Extract unique significant glutathione-related lncRNA IDs
siglncRNAnames <- unique(sig_cor_annot$gene_id)
siglncRNAnames <- sub("\\..*", "", siglncRNAnames)


# Keep only lncRNAs that are present in the expression matrix
siglncRNAnames <- intersect(siglncRNAnames, rownames(exp))


# Subset expression matrix to include only significant glutathione-related lncRNAs
exp_lnc <- exp[siglncRNAnames, , drop = FALSE]


# Transpose expression matrix
exp_lnc_t <- as.data.frame(t(exp_lnc), check.names = FALSE)


# Add sample IDs as a column for merging
exp_lnc_t$sample <- rownames(exp_lnc_t)


# Merge clinical/survival information with lncRNA expression data
surv_data <- merge(
  clinical_full,
  exp_lnc_t,
  by = "sample"
)


#############################3. Univariate Cox regression analysis

# Perform univariate Cox regression for each significant lncRNA
uni_results <- lapply(siglncRNAnames, function(gene) {
 
   cox_formula <- as.formula(
    paste0("Surv(OS.time, OS) ~ `", gene, "`")
  )
  
  cox_model <- coxph(cox_formula, data = surv_data)
  cox_summary <- summary(cox_model)
  
  data.frame(
    gene_id = gene,
    HR = exp(coef(cox_model)),
    p.value = cox_summary$coefficients[1, "Pr(>|z|)"]
  )
})


# Combine all univariate Cox results into one data frame
uni_table <- bind_rows(uni_results)


# Adjust p-values using the Benjamini-Hochberg method
uni_table$FDR <- p.adjust(uni_table$p.value, method = "BH")


# Sort genes by p-value
uni_table <- uni_table[order(uni_table$p.value), ]


# View the most significant prognostic lncRNAs
head(uni_table)


# Select significant lncRNAs from univariate Cox regression
sig_uni <- uni_table %>%
  filter(p.value < 0.05)


# Clean gene IDs in annotation table before merging
sig_lnc$gene_id <- sub("\\..*", "", sig_lnc$gene_id)


# Annotate significant univariate Cox lncRNAs with gene name and gene type
sig_uni_annot <- merge(
  sig_uni,
  sig_lnc[, c("gene_id", "gene_name", "gene_type")],
  by = "gene_id",
  all.x = TRUE
)


# Export significant univariate Cox results
write.csv(
  sig_uni_annot,
  file = "./univariate_glutathione_lncRNAs.csv",
  row.names = FALSE
)


#################################4. Multivariate Cox regression analysis
# Extract lncRNAs that were significant in univariate Cox analysis
sig_uni_genes <- sig_uni_annot$gene_id


# Keep only genes available in the merged survival-expression dataset
sig_uni_genes <- intersect(sig_uni_genes, colnames(surv_data))


# Create multivariate Cox dataset
multivar_data <- surv_data[, c("OS.time", "OS", sig_uni_genes)]


# Construct Cox formula using all significant lncRNAs
cox_formula <- as.formula(
  paste(
    "Surv(OS.time, OS) ~",
    paste(paste0("`", sig_uni_genes, "`"), collapse = " + ")
  )
)


# Fit multivariate Cox regression model
cox_multi <- coxph(cox_formula, data = multivar_data)


# Extract multivariate Cox regression summary
cox_multi_summary <- summary(cox_multi)


# Create multivariate Cox result table
multivar_results <- data.frame(
  gene_id = gsub("`", "", rownames(cox_multi_summary$coefficients)),
  HR = cox_multi_summary$conf.int[, "exp(coef)"],
  Lower95 = cox_multi_summary$conf.int[, "lower .95"],
  Upper95 = cox_multi_summary$conf.int[, "upper .95"],
  p.value = cox_multi_summary$coefficients[, "Pr(>|z|)"]
)


# Adjust p-values for multiple testing
multivar_results$FDR <- p.adjust(multivar_results$p.value, method = "BH")


# Annotate multivariate Cox results with gene name and gene type
multivar_results_annot <- merge(
  multivar_results,
  sig_lnc[, c("gene_id", "gene_name", "gene_type")],
  by = "gene_id",
  all.x = TRUE
)


# Identify independently prognostic lncRNAs
multivar_results_annot$independent <- with(
  multivar_results_annot,
  ifelse(p.value < 0.05 & FDR < 0.05, TRUE, FALSE)
)


# Export multivariate Cox results
write.csv(
  multivar_results_annot,
  file = "./multivariate_glutathione_lncRNAs.csv",
  row.names = FALSE
)


####################### 5. Kaplan-Meier curve for each individual lncRNA
# Make sure exp_mat has clean gene IDs
rownames(exp_mat) <- sub("\\..*", "", rownames(exp_mat))

# Function to generate Kaplan-Meier curve for one lncRNA
make_gene_km <- function(gene_id, gene_label) {
  
  # Extract expression values for the selected lncRNA
  expr_gene <- exp_mat[gene_id, ]
  
  # Convert expression to numeric while keeping sample names
  expr_gene <- as.numeric(expr_gene)
  names(expr_gene) <- colnames(exp_mat)
  
  # Create expression data frame
  df_KM <- data.frame(
    sample = names(expr_gene),
    expression = expr_gene
  )
  
  # Merge expression with clinical survival data
  df_KM <- merge(
    df_KM,
    clinical_full,
    by = "sample"
  )
  
  # Remove missing values
  df_KM <- df_KM %>%
    filter(
      !is.na(expression),
      !is.na(OS.time),
      !is.na(OS)
    )
  
  # Determine optimal expression cut-off using surv_cutpoint
  cut_gene <- surv_cutpoint(
    df_KM,
    time = "OS.time",
    event = "OS",
    variables = "expression"
  )
  
  # Extract cut-off value
  cut_value <- cut_gene$cutpoint[1, "cutpoint"]
  
  # Divide patients into high- and low-expression groups
  df_KM$group_opt <- ifelse(
    df_KM$expression > cut_value,
    "High",
    "Low"
  )
  
  # Fit Kaplan-Meier survival model
  fit_gene <- survfit(
    Surv(OS.time, OS) ~ group_opt,
    data = df_KM
  )
  
  # Plot Kaplan-Meier survival curve
  km_plot <- ggsurvplot(
    fit_gene,
    data = df_KM,
    pval = TRUE,
    risk.table = TRUE,
    title = paste("Kaplan-Meier Curve for", gene_label, "(Optimal Cut-off)"),
    xlab = "Time (days)",
    ylab = "Survival Probability"
  )
  
  print(km_plot)
  
  # Return data, model, and cut-off for later use
  return(list(
    data = df_KM,
    fit = fit_gene,
    cutpoint = cut_value
  ))
}


# Kaplan-Meier curve for lncRNA-X
gene_id_lncRNAX <- "gene id"
km_lncRNAX<- make_gene_km(
  gene_id = gene_id_lncRNAX,
  gene_label = "lncRNA-X"
)

df_KM_lncRNAX <- km_lncRNAX$data


# Kaplan-Meier curve for lncRNA-Y
gene_id_lncRNAY <- "gene id"
km_lncRNAY <- make_gene_km(
  gene_id = gene_id_lncRNAY,
  gene_label = "lncRNA-Y"
)

df_KM_lncRNA-Y <- km_RBM$data


################################ 6. Kaplan-Meier curve for two-lncRNA combined signature
# Cox model coefficients based on hazard ratios
# Formula:
# risk score = coefficient1 * lncRNAX expression +
#              coefficient2 * lncRNAY expression

coef_lncRNAX <- log(1.0434240)    
coef_lncRNAY  <- log(0.9904962)    


# Match samples between expression matrix and clinical data
samples_use <- intersect(
  colnames(exp_mat),
  clinical_full$sample
)


# Extract expression values for both lncRNAs
exprlncRNAX <- exp_mat["gene id", samples_use]
exprlncRNAY  <- exp_mat["gene id", samples_use]


# Calculate two-lncRNA risk score
risk_score <- coef_lncRNAX * as.numeric(exprlncRNAX) +
  coef_lncRNAY * as.numeric(exprlncRNAY)


# Create risk score data frame
df_KM_Both <- data.frame(
  sample = samples_use,
  risk_score = risk_score
)


# Merge risk score with clinical survival data
df_KM_Both <- merge(
  df_KM_Both,
  clinical_full,
  by = "sample"
)


# Remove missing values before survival analysis
df_KM_Both <- df_KM_Both %>%
  filter(
    !is.na(risk_score),
    !is.na(OS.time),
    !is.na(OS)
  )


# Determine optimal cut-off for the combined risk score
cutpoint_Both <- surv_cutpoint(
  df_KM_Both,
  time = "OS.time",
  event = "OS",
  variables = "risk_score"
)


# Extract cut-off value
risk_cutoff <- cutpoint_Both$cutpoint[1, "cutpoint"]


# Divide patients into high- and low-risk groups
df_KM_Both$risk_group <- ifelse(
  df_KM_Both$risk_score > risk_cutoff,
  "High",
  "Low"
)


# Fit Kaplan-Meier model for combined risk signature
fit_both <- survfit(
  Surv(OS.time, OS) ~ risk_group,
  data = df_KM_Both
)


# Plot Kaplan-Meier curve for two-lncRNA signature
ggsurvplot(
  fit_both,
  data = df_KM_Both,
  pval = TRUE,
  risk.table = TRUE,
  title = "Kaplan-Meier Curve for Two-lncRNA Combined Signature",
  xlab = "Time (days)",
  ylab = "Survival Probability"
)


################################### 7. Time-dependent ROC curve
# Generate time-dependent ROC curve for 1-, 3-, and 5-year survival prediction
roc_os <- timeROC(
  T = df_KM_Both$OS.time,
  delta = df_KM_Both$OS,
  marker = df_KM_Both$risk_score,
  cause = 1,
  times = c(365, 1095, 1825),
  iid = TRUE
)


# Plot 1-year ROC curve
plot(
  roc_os,
  time = 365,
  col = "#D55E00",
  lwd = 2
)


# Add 3-year ROC curve
plot(
  roc_os,
  time = 1095,
  col = "#0072B2",
  lwd = 2,
  add = TRUE
)


# Add 5-year ROC curve
plot(
  roc_os,
  time = 1825,
  col = "#009E73",
  lwd = 2,
  add = TRUE
)


# Add legend showing AUC values
legend(
  "bottomright",
  legend = c(
    paste0("1-year AUC = ", round(roc_os$AUC[1], 3)),
    paste0("3-year AUC = ", round(roc_os$AUC[2], 3)),
    paste0("5-year AUC = ", round(roc_os$AUC[3], 3))
  ),
  col = c("#D55E00", "#0072B2", "#009E73"),
  lwd = 2,
  bty = "n"
)


################################# 8. Stratified KM analysis: lncRNA-X and lncRNA-Y groups
# Create simplified lncRNA-X grouping data
df_lncRNAX_group <- df_KM_lncRNAX %>%
  select(sample, OS.time, OS, lncRNAX_group = group_opt)


# Create simplified lncRNAY grouping data
df_lncRNAY_group <- df_KM_lncRNAY %>%
  select(sample, lncRNAY_group = group_opt)


# Merge grouping information from both lncRNAs
df_stratify4groups <- df_lncRNAX_group %>%
  inner_join(df_lncRNAY_group, by = "sample")


# Create four expression-combination groups
df_stratify4groups <- df_stratify4groups %>%
  mutate(
    Group4 = case_when(
      lncRNAX_group == "Low"  & lncRNAY_group == "Low"  ~ "lncRNAX Low / lncRNAY Low",
      lncRNAX_group == "Low"  & lncRNAY_group == "High" ~ "lncRNAX Low / lncRNAY High",
      lncRNAX_group == "High" & lncRNAY_group == "Low"  ~ "lncRNAX High / lncRNAY Low",
      lncRNAX_group == "High" & lncRNAY_group == "High" ~ "lncRNAX High / lncRNAY High"
    )
  )


# Fit Kaplan-Meier model using the four combined expression groups
fit_stratify4groups <- survfit(
  Surv(OS.time, OS) ~ Group4,
  data = df_stratify4groups
)


# Plot stratified Kaplan-Meier curve
ggsurvplot(
  fit_stratify4groups,
  data = df_stratify4groups,
  pval = TRUE,
  conf.int = TRUE,
  risk.table = TRUE,
  palette = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3"),
  title = "Stratified Kaplan-Meier Curve: lncRNA-X and lncRNA-Y",
  xlab = "Time (days)",
  ylab = "Survival Probability"
)