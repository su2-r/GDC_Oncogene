library(tidyverse)
install.package("BiocManager")
install.packages("TCGAbiolinks")

query <- GDCquery(project = "TCGA-LUAD",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq-FPKM-UQ")
GDCdownload(query)
d_luad0 <- GDCprepare(query)
d_luad = as.data.frame(d_luad0@colData)

# Specify access level (replace with "controlled" if needed), 
# repeat this code for additional data cancer type to compare from GDc portal
# in this case study use TCGA-Lung Adenocarcinoma (LUAD) and TCGA-Lung Squamous Cell Carcinoma (LUSC).
query <- GDCquery(project = "TCGA-LUAD",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "HTSeq-npes-tr",access = "Open")
# Download data
tryCatch({
  GDCdownload(query)
}, error = function(error) {
  message("Error downloading data:", error)
})

# Prepare data
d_luad0 <- GDCprepare(query)

# Convert to data frame
d_luad <- as.data.frame(d_luad0@colData)

# Load the required library
library(SummarizedExperiment)

# Extract gene expression data for LUAD and LUSC
e_luad <- assay(d_luad0)
e_lusc <- assay(d_lusc0)

# Filter for genes on chromosome 1
g_luad <- d_luad0@rowRanges %>% 
  as.data.frame() %>% 
  filter(seqnames == 'chr1') %>% 
  pull(ensembl_gene_id)

g_lusc <- d_lusc0@rowRanges %>% 
  as.data.frame() %>% 
  filter(seqnames == 'chr1') %>% 
  pull(ensembl_gene_id)

# Subset gene expression data to chromosome 1 genes
e_luad <- e_luad[g_luad, ]
e_lusc <- e_lusc[g_lusc, ]

# Save the data
save(d_luad, d_lusc, e_luad, e_lusc, file = 'data.TCGA_LUAD_LUSC.gene_expression.Rdata')

#load data 
head(d_luad)
summary(e_luad)
load('data.TCGA_LUAD_LUSC.gene_expression.Rdata')

#explore
colnames(d_luad)
colnames(d_lusc)

#analize pathologyTP: Primary solid, TumorTR: Recurrent solid, TumorNT: Solid Tissue Normal
#sample code https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes"
d_luad %>% count(shortLetterCode)
d_lusc %>% count(shortLetterCode)

# use code snippet x number of data sets
d_luad %>% 
  count(shortLetterCode) %>% 
  ggplot(., aes(shortLetterCode, n, fill = shortLetterCode)) + 
  geom_bar(stat = "identity", position = position_dodge())

d_luad %>% 
  count(shortLetterCode) %>% 
  ggplot(., aes(shortLetterCode, n, fill = shortLetterCode)) + 
  geom_bar(stat = "identity", position = position_dodge())

bind_rows(
  d_luad %>% 
    mutate(type = 'luad') %>%
    select(type, shortLetterCode, tumor_stage),
  d_lusc %>% 
    mutate(type = 'lusc') %>%
    select(type, shortLetterCode, tumor_stage)
) %>%
  count(shortLetterCode, type) %>% 
  complete(type, shortLetterCode, fill = list(n = 0)) %>%  # Fills in missing combinations with 0 counts
  ggplot(., aes(shortLetterCode, n, fill = shortLetterCode)) + 
  geom_bar(stat = "identity", position = position_dodge()) + 
  facet_wrap(~type, scales = 'free_x', ncol = 5)

# Plot histogram
ggplot(d_luad, aes(age_at_diagnosis)) + geom_histogram(bins=100)

# Change the axis
ggplot(d_luad, aes(age_at_diagnosis/365)) + geom_histogram(bins=100)

ggplot(d_luad, aes(age_at_diagnosis)) + geom_density()

# Histogram and density plot
ggplot(d_luad, aes(age_at_diagnosis / 365)) +
  geom_histogram(bins = 100, color = "black", fill = "white") +
  geom_density(alpha = 0.2, fill = "#FF6666")

# Histogram with density overlay
ggplot(d_luad, aes(age_at_diagnosis / 365)) +
  geom_histogram(aes(y = ..density..), color = "black", fill = "white") +
  geom_density(alpha = 0.2, fill = "#FF6666")

# Histogram of age at diagnosis by gender
ggplot(d_luad, aes(age_at_diagnosis / 365, fill = gender)) + 
  geom_histogram(bins = 100)

# Boxplot of age at diagnosis by gender
ggplot(d_luad, aes(gender, age_at_diagnosis / 365, fill = gender)) + 
  geom_boxplot()

# Violin plot of age at diagnosis by gender
ggplot(d_luad, aes(gender, age_at_diagnosis / 365, fill = gender)) + 
  geom_violin()

# Combined violin and boxplot
ggplot(d_luad, aes(gender, age_at_diagnosis / 365)) + 
  geom_violin(aes(fill = gender)) + 
  geom_boxplot(fill = 'white', width = 0.25)

# Combined violin, boxplot, and dot plot
ggplot(d_luad, aes(gender, age_at_diagnosis / 365)) + 
  geom_violin(aes(fill = gender)) + 
  geom_boxplot(fill = 'white') + 
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 1, alpha = 0.5)

ggplot(d_luad, aes(age_at_diagnosis / 365)) + 
  geom_histogram(bins = 100) +
  labs(
    x = "Age at Diagnosis (years)",
    y = "Number of LUAD patients",
    title = "TCGA: Lung Cancer Adenocarcinoma"
  )








