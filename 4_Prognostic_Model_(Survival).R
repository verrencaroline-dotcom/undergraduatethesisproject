######### Prognostic Model Construction and Validation 

#Preparing clinical and survival data

##importing clinical data 
library(survminer)
library(survival)
library(dplyr)
library(stringr)

clinical <- read.delim("clinical.tsv", sep="\t", header=TRUE)
clinical_covariates <- clinical_covariates <- clinical[, c(
  "sample",
  "age_at_earliest_diagnosis_in_years.diagnoses.xena_derived",
  "ajcc_pathologic_t.diagnoses",
  "ajcc_pathologic_n.diagnoses",
  "ajcc_pathologic_m.diagnoses",
  "ajcc_pathologic_stage.diagnoses"
   )]
colnames(clinical_covariates) <- c(
  "sample", "age", "T_stage", "N_stage", "M_stage", "Stage")

##merging survival and clinical data
clinical_covariates$sample <- gsub("A$", "", clinical_covariates$sample)
clinical_full <- merge(survival_brac, clinical_covariates, by = "sample")
clinical_full <- clinical_full %>% 
  select(-DSS,-DSS.time,-DFI, -DFI.time, -PFI, -PFI.time)

##formatting and ensuring no missing data for important covariates
clinical_full <- clinical_full %>%
  filter(!is.na(OS) & !is.na(OS.time))
clinical_full <- clinical_full %>%
  filter(!is.na(Stage) & !is.na(age))
clinical_full$Stage <- toupper(clinical_full$Stage)
clinical_full$Stage <- gsub("STAGE ", "", clinical_full$Stage)
clinical_full$Stage <- gsub("[^A-Z0-9]", "", clinical_full$Stage)
clinical_full$Stage <- factor(
  clinical_full$Stage,
  levels = c("I", "II", "IIA", "IIB", "III", "IIIA", "IIIB", "IV"),
  ordered = TRUE
   )
clinical_full$age <- as.numeric(clinical_full$age)
clinical_full$OS.time <- as.numeric(clinical_full$OS.time)
clinical_full$OS <- as.numeric(clinical_full$OS)
clinical_full <- clinical_full[!duplicated(clinical_full$sample), ]

## subsetting data and making full survival data with significant lncRNAs
siglncRNAnames <- unique(sig_cor_annot$gene_id)
exp_lnc <- exp[siglncRNAnames, ]
exp_lnc_t <- as.data.frame(t(exp_lnc))
exp_lnc_t$sample <- rownames(exp_lnc_t)
surv_data <- merge(clinical_full, exp_lnc_t, by = "sample")

## univariate analysis
uni_results <- lapply(siglncRNAnames, function(gene) {
  coxph(
    Surv(OS.time, OS) ~ surv_data[[gene]],
    data = surv_data
  )
})
uni_table <- data.frame(
  gene = siglncRNAnames,
  HR = sapply(uni_results, function(x) exp(coef(x))),
  p.value = sapply(uni_results, function(x)
    summary(x)$coefficients[,"Pr(>|z|)"]
  )
)

# Multiple testing correction (recommended)
uni_table$FDR <- p.adjust(uni_table$p.value, method = "BH")

# Sort by significance
uni_table <- uni_table[order(uni_table$p.value), ]
head(uni_table)

#make table and annotation 
colnames(sig_uni)[colnames(sig_uni) == "gene"] <- "gene_id"
sig_uni_annot <- merge(sig_uni, sig_lnc[, c("gene_id", "gene_name", "gene_type")])
write.csv(write.csv(sig_uni_annot,"./univariate glutathione lncRNAs.csv"))

##Multivariate analysis
sig_uni_genes <- sig_uni_annot$gene_id
multivar_data <- surv_data[, c("OS.time", "OS", sig_uni_genes)]
cox_formula <- as.formula(
       paste("Surv(OS.time, OS) ~", paste(sig_uni_genes, collapse = " + ")))
cox_multi <- coxph(cox_formula, data = multivar_data)
multivar_results <- data.frame(
       gene = rownames(summary(cox_multi)$coefficients),
       HR = summary(cox_multi)$conf.int[,"exp(coef)"],
       Lower95 = summary(cox_multi)$conf.int[,"lower .95"],
       Upper95 = summary(cox_multi)$conf.int[,"upper .95"],
       p.value = summary(cox_multi)$coefficients[,"Pr(>|z|)"])
multivar_results$FDR <- p.adjust(multivar_results$p.value, method = "BH")

colnames(multivar_results)[colnames(multivar_results) == "gene"] <- "gene_id"
multivar_results_annot <- merge(multivar_results, sig_lnc[, c("gene_id", "gene_name", "gene_type")])
multivar_results_annot$independent <- with(multivar_results, ifelse(p.value < 0.05 & FDR < 0.05, TRUE, FALSE))

##Kaplan-Meier Curve for each individual lncRNA
#For ACAP2
gene_id_ACAP <- "ENSG00000229325"
expr_ACAP <- as.numeric(exp_mat[rownames(exp_mat) == gene_id_ACAP, ])
names(expr_ACAP) <- colnames(exp_mat)

df_KM_ACAP <- data.frame(
  sample = names(expr_ACAP),
  expression = expr_ACAP
  )

df_KM_ACAP <- merge(df_KM_ACAP, clinical_full, by = "sample")

cut_ACAP <- surv_cutpoint(
  df_KM_ACAP,
  time = "OS.time",
  event = "OS",
  variables = "expression"
  )
df_KM_ACAP$group_opt <- ifelse(df_KM_ACAP$expression > cut_ACAP$cutpoint$cutpoint, "High", "Low")
fit_ACAP <- survfit(Surv(OS.time, OS) ~ group_opt, data = df_KM_ACAP)

ggsurvplot(
  fit_ACAP, data = df_KM_ACAP, pval = TRUE, risk.table = TRUE,
  title = paste("KM Curve for", gene_id_ACAP, "(Optimal Cutoff)"),
  xlab = "Time",
  ylab = "Survival Probability"
  )

#For RBMA5S1 (just replace everything on top for ACAP with RBM5AS1 info)

##################Kaplan Meier Curve 
coef1 <- log(1.0434240)    # ENSG00000229325
coef2 <- log(0.9904962)    # ENSG00000281691

exprACAP <- as.numeric(exp_mat["ENSG00000229325", ])
exprRBM <- as.numeric(exp_mat["ENSG00000281691", ])

samples_use <- intersect(colnames(exp_mat), clinical_full$sample)
exprACAP <- exprACAP[samples_use]
exprRBM <- exprRBM[samples_use]

risk_score <- coefACAP * exprACAP + coefRBM * exprRBM

df_KM_Both <- data.frame(
  sample = samples_use,
  risk_score = risk_score
)
df_KM_Both <- merge(df_KM_Both, clinical_full, by = "sample")

cutpoint_Both <- surv_cutpoint(
  df_KM_Both,
  time = "OS.time",
  event = "OS",
  variables = "risk_score"
)

df_KM_Both$risk_group <- ifelse(df$risk_score > cutpoint_Both$cutpoint$cutpoint, "High", "Low")

fit_both <- survfit(Surv(OS.time, OS) ~ risk_group, data = df_KM_Both)

ggsurvplot(
  fit_both, data = df_KM_Both, pval = TRUE, risk.table = TRUE,
  title = "KM Curve for 2-lncRNA Combined Signature",
  xlab = "Time",
  ylab = "Survival Probability"
)

####################### ROC Curve
library(timeROC)

roc_os <- timeROC(
  T = df$OS.time,
  delta = df_KM_Both$OS,
  marker = df_KM_Both$risk_score,
  cause = 1,
  times = c(365, 1095, 1825),
  iid = TRUE
  )

plot(roc_os, time = 365, col = "#D55E00", lwd = 2)
plot(roc_os, time = 1095, col = "#0072B2", lwd = 2, add = TRUE)
plot(roc_os, time = 1825, col = "#009E73", lwd = 2, add = TRUE)
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

################## stratified ACAP2-IT1 & RBM5-AS1 groups KM 

df_stratify4groups <- df_KM_ACAP %>%
  inner_join(df_KM_RBM, by="sample")
df_stratify4groups <- df_stratify4groups %>%
  mutate(
    Group4 = case_when(
    group_opt.x=="Low" & group_opt.y=="Low"   ~ "ACAP2 Low / RBM5 Low",
    group_opt.x=="Low" & group_opt.y=="High"  ~ "ACAP2 Low / RBM5 High",
    group_opt.x=="High" & group_opt.y=="Low"  ~ "ACAP2 High / RBM5 Low",
    group_opt.x=="High" & group_opt.y=="High" ~ "ACAP2 High / RBM5 High"))

fit_stratify4groups <- survfit(Surv(OS.time.x, OS.x) ~ Group4, data=df_stratify4groups)

ggsurvplot(
  fit, data=df_stratify4groups,
  pval=TRUE, conf.int=TRUE,
  risk.table=TRUE,
  palette=c("#E41A1C","#377EB8","#4DAF4A","#984EA3"),
  title="Stratified KM: ACAP2-IT1 & RBM5-AS1"
  )


################## plots
## survival status plot
p_survivstat <- ggplot(df_KM_Both, aes(x = seq_along(OS.time), y = OS.time, color = factor(OS))) +
  geom_point(size = 1.5) +
  scale_color_manual(values = c("0" = "blue", "1" = "red"),
                     labels = c("Alive", "Dead")) +
  labs(x = "Patients (low → high risk)", y = "Survival time (days)") +
  theme_minimal()
p_survivstat

##risk score distribution 
