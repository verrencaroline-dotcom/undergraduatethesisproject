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

# glutathione annotation imports
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

##importing glutathione gene list 
gsh_genes <- read.table("gsh_genes.txt", header = FALSE, stringsAsFactors = FALSE)[,1]

## Map glutathione genes to gene names of annotated data -> to match with gene IDs 
gsh_annot <- annot[annot$gene_name %in% gsh_genes, ]
## Extract corresponding gene_id to that of the gsh gene names 
gsh_ids <- gsh_annot$gene_id
##removing version . so that gene_ids match
gsh_ids <- gsub("\\..*", "", gsh_ids)

# Remove version numbers from the row names of the expression matrix
rownames(exp) <- sub("\\..*", "", rownames(exp))


# Convert the expression data frame into a numeric matrix
exp_mat <- as.matrix(exp)
mode(exp_mat) <- "numeric"


# Remove genes with zero variance across all samples
exp_mat <- exp_mat[rowVars(exp_mat) > 0, ]


# Keep only glutathione-related genes that are present in the expression matrix
gsh_ids <- intersect(gsh_ids, rownames(exp_mat))


# Check how many glutathione-related genes are available for ssGSEA analysis
length(gsh_ids)


# Create a gene set list for ssGSEA
gene_sets <- list(
  Glutathione_Metabolism = gsh_ids
)


# Create the ssGSEA parameter object
ssgsea_param <- ssgseaParam(
  expr = exp_mat,
  geneSets = gene_sets
)


# Run ssGSEA using the GSVA package
ssgsea_scores <- gsva(ssgsea_param)


# Check the dimensions of the ssGSEA score matrix
dim(ssgsea_scores)


# View the first few ssGSEA scores
head(ssgsea_scores)


###################################################################
# Compare glutathione metabolism ssGSEA scores between normal and tumor breast samples

# Import sample category information for TCGA and GTEx samples
iden <- read.table(
  "D:/Thesis/Raw Data/TCGA_GTEX_category.txt",
  header = TRUE,
  sep = "\t"
)


# Select sample IDs corresponding to normal breast tissue from GTEx
normal_breast_samples <- iden$sample[
  iden$TCGA_GTEX_main_category == "GTEX Breast"
]


# Subset ssGSEA scores to include only normal breast samples
ssgsea_normal_breast <- ssgsea_scores[
  , colnames(ssgsea_scores) %in% normal_breast_samples
]


# View the normal breast ssGSEA scores
head(ssgsea_normal_breast)


# Select sample IDs corresponding to breast cancer samples from TCGA
tumor_breast_samples <- iden$sample[
  iden$TCGA_GTEX_main_category == "TCGA Breast Invasive Carcinoma"
]


# Subset ssGSEA scores to include only tumor breast samples
ssgsea_tumor_breast <- ssgsea_scores[
  , colnames(ssgsea_scores) %in% tumor_breast_samples
]


# View the tumor breast ssGSEA scores
head(ssgsea_tumor_breast)


# Create a data frame containing ssGSEA scores and sample group labels
df_ssgsea <- data.frame(
  score = c(ssgsea_normal_breast, ssgsea_tumor_breast),
  group = factor(
    c(
      rep("Normal", length(ssgsea_normal_breast)),
      rep("Tumor",  length(ssgsea_tumor_breast))
    ),
    levels = c("Normal", "Tumor")
  )
)


# Perform a Wilcoxon rank-sum test
wilcox.test(score ~ group, data = df_ssgsea)


# Calculate the median ssGSEA score for each group
tapply(df_ssgsea$score, df_ssgsea$group, median)


# Create a boxplot comparing glutathione metabolism ssGSEA scores between normal and tumor breast samples
ggplot(df_ssgsea, aes(x = group, y = score, fill = group)) +
  
  # Add boxplots to show the median, quartiles, and overall distribution
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  
  # Add jittered points to show individual sample ssGSEA scores
  geom_jitter(width = 0.2, size = 1.2, alpha = 0.6) +
  
  # Add Wilcoxon test p-value to the plot
  stat_compare_means(
    method = "wilcox.test",
    label = "p.format"
  ) +
  
  # Use a clean classic theme
  theme_classic() +
  
  # Add plot title and axis labels
  labs(
    title = "Glutathione Metabolism ssGSEA Scores",
    y = "ssGSEA score",
    x = ""
  ) +
  
  # Remove the legend because the x-axis already shows the groups
  theme(legend.position = "none")


# Create a data frame containing ssGSEA scores for normal breast samples only
df_normal_ssgsea <- data.frame(
  SampleID = colnames(ssgsea_normal_breast),
  Glutathione_Metabolism = as.numeric(ssgsea_normal_breast)
)


# Export the normal breast ssGSEA scores as a CSV file
write.csv(
  df_normal_ssgsea,
  file = "Normal_Breast_ssGSEA.csv",
  row.names = FALSE
)


# Test whether normal breast ssGSEA scores follow a normal distribution
shapiro_test <- shapiro.test(as.numeric(ssgsea_normal_breast))


# Create a QQ plot to visually assess normality of normal breast ssGSEA scores
qqnorm(
  as.numeric(ssgsea_normal_breast),
  main = "QQ Plot: Normal Breast ssGSEA Scores"
)

# Add a reference line to the QQ plot
qqline(
  as.numeric(ssgsea_normal_breast),
  col = "red",
  lwd = 2
)


# Create a histogram to visualize the distribution of normal breast ssGSEA scores
hist(
  as.numeric(ssgsea_normal_breast),
  breaks = 15,
  col = "grey",
  border = "black",
  main = "Histogram: Normal Breast ssGSEA Scores",
  xlab = "Glutathione Metabolism ssGSEA Score",
  ylab = "Frequency"
)
