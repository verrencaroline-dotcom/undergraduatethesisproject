# read data
set.seed(1)
library(BiocManager)
library(dplyr)

## reading the TCGA & GTEX labels
iden<- read.table('TCGA_GTEX_category.txt', header= TRUE, sep = '\t') 

## reading the TCGA & GTEX expression data 
exp<- read.table("./exp_new.txt", row.names = 1)
    

#############################################
#take breast data (both normal tissue and cancer tissue)
## filtering for only breast data labels 
brac_id<- iden %>% filter(grepl("breast", TCGA_GTEX_main_category, ignore.case = TRUE))

## importing survival information for GTEX and TCGA from all tissue
survival <- read.table("./survival.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

## filtering only for breast-related survival information 
survival_brac<- survival %>% filter(sample %in% brac_id$sample)

## replaces the column names of the dataframe with the values in the first row because he first row often contains sample IDs or column headers that weren’t properly set when reading the file.
colnames(exp)<- exp[1,]

## This removes the first row (which has just been used as column names).[-1,] means “all rows except the first one”.
exp<- exp[-1,]

## filters exp data for only breast-related samples 
brac_exp <- exp[, colnames(exp) %in% brac_id$sample]

## line simply removes the variable exp from your R environment (to free memory or avoid confusion). After this command, exp no longer exists — only brac_exp remains.
rm(exp)

## converts all columns in brac_exp to numeric values
exp <- as.data.frame(lapply(brac_exp, function(x) as.numeric(as.character(x))))

## restores row names (the gene labels)
rownames(exp)<- rownames(brac_exp)

## restores column names 
colnames(exp)<- colnames(brac_exp)

## removes brac_exp
rm(brac_exp)

## writes the contents of brac_id$sample (a vector of sample IDs) into a text file named exp_id.txt. so text is filled with sample ID.
write(brac_id$sample, "./exp_id.txt")

##back-log transforming data (because the data is in logged form (0-15 range))
exp <- round(((2^exp) - 1), 0)

#############################################
#filtering the tcga samples  

## samples are the column names 
samples <- colnames(exp)

##take samples that only start with TCGA string
tcga_samples <- grep("^TCGA", samples, value = TRUE)

## This extracts the two-character sample type code from the end of each TCGA ID. The two-character sample code identifies type of tumor and if its normal tissue or not. 
tcga_sample_types <- sapply(tcga_samples, function(x) substr(x, nchar(x)-1, nchar(x)))

## displays how many samples are in each group code (normal vs. tumor)
table(tcga_sample_types)
### results = 01 = primary tumor (1092), 06 = metastatic (7)

## finds all non-TCGA columns in your expression matrix =  GTEx samples.
non_tcga_samples <- colnames(exp)[!grepl("^TCGA", colnames(exp))]

## only use primary tumors 01 as metastatic is little and could bias the results 
tcga_tumor_samples <- tcga_samples[tcga_sample_types == "01"]

## changing the exp data to contain filtered data 
samples_to_keep <- c(tcga_tumor_samples, non_tcga_samples)
exp <- exp[, samples_to_keep]

##removing all unnecesssary files as we already have clean exp data. 
rm(tcga_samples, tcga_sample_types, non_tcga_samples,
   tcga_tumor_samples, samples_to_keep)

## writes the contents of column names of exp (a vector of sample IDs) into a text file named exp_used.txt. 
write(colnames(exp), "./exp_used.txt")

#############################################
#grouping of normal and tumor samples for DGE later
library(limma)

## samples are column names of the expression data 
samples <- colnames(exp)

## main data assumes TCGA samples = tumor and all others = normal.
metadata <- data.frame(
  sample = samples,
  condition = ifelse(grepl("TCGA", samples), "Tumor", "Normal")
)
##row names of metadata are the sample
rownames(metadata) <- metadata$sample

##defining factors = tells R which samples are Normal and which are Tumor.
metadata$condition <- factor(metadata$condition,
                             levels = c("Normal", "Tumor"))

##designing model matrix
design <- model.matrix(~ 0 + condition, data = metadata)
##column names of model matrix based on factor 
colnames(design) <- c("Normal", "Tumor")  

## we need to tell limma what comparison we actually want. 
contrast_matrix <- makeContrasts(
  TumorVsNormal = Tumor - Normal,  
  levels = design
  )

#############################################
# cleaning data and normalization of data for differential gene expression
library(edgeR)

## making dge object 
dge <- DGEList(counts = exp, group = metadata$condition)
##filtering/removes gene with very low expression
keep <- filterByExpr(
  dge, 
  design = design,
  min.count = 5,        # Reduce from default 10 → 5 (or lower)
  min.total.count = 10, # Reduce total count cutoff if needed
  large.n = 10,         # Require expression in fewer samples
  min.prop = 0.25       # Keep genes expressed in 25% of smallest group
)
##keeping only filtered genes to be statistically stronger
dge <- dge[keep, , keep.lib.sizes = FALSE]

## normalize counts so that counts are comparable across samples
dge <- calcNormFactors(dge, method = "TMM")  

#############################################
# DGE analysis
## voom transformation. prepares RNA-seq data for limma linear modeling by converts your RNA-seq counts into log2-counts per million (logCPM) while estimating the mean-variance relationship.
v <- voom(dge, design, plot = TRUE)

## fits a linear model for each gene using the voom-transformed data
fit <- lmFit(v, design)
## applies the contrast(s) you defined (Tumor vs Normal)
fit <- contrasts.fit(fit, contrast_matrix)
## performs empirical Bayes shrinkage of the standard errors across all genes
fit <- eBayes(fit)

#############################################
# lnc annotation imports
library(rtracklayer)
library(stringr)

## Import GTF annotation. reads a GTF (Gene Transfer Format) file. contains genome annotations (what is miRNA, lincRNA, etc.)
gtf_data <- import("gencode.v22.annotation.gtf", format = "gtf")
head(gtf_data)

## converts the GRanges object to a regular R data frame, easier to manipulate.
gtf_df <- as.data.frame(gtf_data)
rm(gtf_data)
head(gtf_df)

## create a clean reference table of genes for annotation. selecting only three columns. removes duplicate rows so that each gene appears only once.Creates a clean annotation table where each row corresponds to a unique gene with its ID, name, and type.
annot<- unique(gtf_df[, c("gene_id", "gene_name", "gene_type")])
rm(gtf_df)
head(annot)

## lnc_annot that only a data frame containing only lncRNAs and their subtypes.You can use this to filter your DEGs for lncRNAs or perform downstream analyses specific to lncRNAs. 
lnc_annot<-annot%>%filter(gene_type == "lincRNA" | 
                              gene_type == "sense_intronic" |
                              gene_type == "sense_overlapping" |
                              gene_type == "antisense" |
                              gene_type == "processed_transcript" |
                              gene_type == "3prime_overlapping_ncrna")
                            
#############################################
#clean DEG result

## extracts DGE results from limma fit.
DGEresults <- topTable(fit, coef = "TumorVsNormal", number = Inf, adjust.method = "fdr")
head(DGEresults)

## removing version number from Ensembl IDs
DGEresults$gene_id <- sub("\\..*", "", rownames(DGEresults))

## make sure lnc_annot has matching IDs 
lnc_annot$gene_id <- sub("\\..*", "", lnc_annot$gene_id)
# annot$gene_id_2 <- sub("\\..*", "", annot$gene_id)

## merging DEG table with lnc_annot table
lnc_DEGs <- merge(DGEresults, lnc_annot, by = "gene_id")
#res_anot <- merge(results, annot, by.x = "gene_id", by.y = "gene_id_2", all.x = TRUE)
write.csv(lnc_DEGs,"./DEG_alllncRNAs.csv")
#write.csv(res_anot,"./res/deg/DEG_all.csv")

##finding significantly different lncRNAs between normal and tumor samples. extra filter
sig_genes <- lnc_DEGs[lnc_DEGs$adj.P.Val < 0.05 & abs(lnc_DEGs$logFC) > 1, ]
#sig_genes <- res_anot[res_anot$adj.P.Val < 0.05 & abs(res_anot$logFC) > 1, ]
dim(sig_genes)
head(sig_genes)

pos_gene<-sig_genes[sig_genes$logFC>1,]

write.csv(sig_genes,"./DEG_siglncRNAs.csv")

##################################################
#checking and handling missing values

## counts how many missing values are in each column (logFC, P.Value, gene_name)
sum(is.na(lnc_DEGs$logFC))
sum(is.na(lnc_DEGs$P.Value))
sum(is.na(lnc_DEGs$gene_name))
#sum(is.na(res_anot$logFC))
#sum(is.na(res_anot$P.Value))
#sum(is.na(res_anot$gene_name))

#note: no missing values so don't do next things
#missing<-res_anot[is.na(res_anot$gene_name), ]

#res <- res_anot[!is.na(res_anot$gene_name), ]

#write.csv(missing,"./res/unknown.csv")

##################################################
#Volcano plot for differentially expressed lncRNAs in tumor vs normal tissue
library(EnhancedVolcano)

# Select the top 10 lncRNAs with the smallest adjusted p-values
top<-head(lnc_DEGs[order(lnc_DEGs$adj.P.Val),], 10)

# Replace extremely small raw p-values with 1e-50
lnc_DEGs$P.Value[lnc_DEGs$P.Value < 1e-50] <- 1e-50

# Open a new plotting window with a specified width and height
dev.new(width = 8, height = 6)

# Generate the volcano plot showing differential expression of lncRNAs
EnhancedVolcano(
  lnc_DEGs,
  lab = lnc_DEGs$gene_name,
  x = 'logFC',
  y = 'P.Value',
  title = 'Volcano Plot: Differential Expression of lncRNAs',
  subtitle = 'Adjusted p < 0.001 | |log2FC| > 1',
  pCutoff = 1e-3,
  FCcutoff = 1,
  pointSize = 1.2,
  labSize = 2.5,
  colAlpha = 0.7,
  drawConnectors = TRUE,
  widthConnectors = 0.3,
  max.overlaps = 50,
  selectLab = top$gene_name,  # Automatically label top 10 genes
  ylim = c(0, 60),
  col = c('grey80', 'skyblue3', 'forestgreen', 'firebrick'),
  boxedLabels = FALSE,
  legendPosition = 'right',
  legendLabSize = 9,
  legendIconSize = 3.0,
  gridlines.major = FALSE,
  gridlines.minor = FALSE,
  titleLabSize = 12,
  subtitleLabSize = 9,
  axisLabSize = 10
)

# Save the volcano plot as a high-resolution PNG & PDF file
ggsave("volcano_plot.png", width = 8, height = 6, dpi = 300)
ggsave("volcano_plot.pdf", width = 8, height = 6, dpi = 300)

rm(top)

# Boxplot of expression levels for selected significant lncRNAs
library(reshape2)
library(ggplot2)

# Select the top 10 significant genes based on adjusted p-value
top <- head(sig_genes[order(sig_genes$adj.P.Val), "gene_id"], 10)

# Subset the expression matrix to include only the selected top 10 genes
bp_dat <- exp[top, ]
bp_dat$Probe <- rownames(bp_dat)

# Merge expression data with lncRNA annotation data
bp_dat$Probe <- sub("\\.\\d+$", "", bp_dat$Probe)  
bp_dat <- merge(bp_dat, lnc_annot[, c("gene_id", "gene_name")],
                by.x = "Probe", by.y = "gene_id", all.x = TRUE)
head(bp_dat[, c("Probe", "gene_name")])

# Reshape the expression data from wide format to long format for plotting with ggplot2
bp_long <- melt(bp_dat, id.vars = c("Probe", "gene_name"), 
                variable.name = "Sample", value.name = "Expression")

# Add sample group information
bp_long$Group <- ifelse(grepl("TCGA", bp_long$Sample), "Tumor", "Normal")

# Perform independent t-tests for each gene
p_values <- bp_long %>%
  group_by(Probe, gene_name) %>%
  summarize(
    pvalue = tryCatch(t.test(Expression ~ Group)$p.value, error = function(e) NA),
    .groups = "drop"
  )

# Create facet labels containing the gene name and corresponding t-test p-value 
p_values <- p_values %>%
  mutate(
    label = ifelse(
      is.na(pvalue),
      paste0(gene_name, " (p = NA)"),
      paste0(gene_name, " (p = ", signif(pvalue, 3), ")")
    )
  )

#  Merge the p-value labels back into the long-format expression data
bp_long <- merge(bp_long, p_values[, c("Probe", "label")], by = "Probe")

# Create boxplots comparing expression levels between Tumor and Normal samples
boxplot <- ggplot(bp_long, aes(x = Group, y = Expression, fill = Group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(position = position_jitter(0.2), size = 0.5, alpha = 0.6) +
  facet_wrap(~ label, scales = "free_y") + # Facet titles include p-values
  scale_fill_manual(values = c("Normal" = "skyblue", "Tumor" = "tomato")) +
  labs(
    title = "Tumor vs Normal Expression Levels",
    x = "Group",
    y = "Expression"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    text = element_text(size = 14),
    strip.text = element_text(size = 10) # Adjust facet title size
  )

plot(boxplot)