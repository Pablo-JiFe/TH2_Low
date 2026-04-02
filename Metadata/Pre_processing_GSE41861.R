# Script to perform normalization, annotation and batch correction for PCA and cluster analysis

library(oligo)
library(GEOquery)
library(tidyverse)
library(limma)
library(hgu133plus2.db)
library(gridExtra)



# 1.- Load data and metadata ----------------------------------------------

# 1.1 Get supplementary files
# The use of supp files (.cel) is for future analysis of multiple data bases

#getGEOSuppFiles("GSE41861", baseDir = "D:/th2_low/GSE41861")

# 1.1.2 Untar files

#untar("D:/th2_low/GSE41861/GSE41861_RAW.tar", exdir = "Data/GSE41861/")

# 1.1.3 Listing .cel files

cel_files <- list.celfiles("Data/GSE41861/", full.names = TRUE, listGzipped = TRUE)

# 1.1.4 Reading in cel files

raw_data <- read.celfiles(cel_files)

# 1.2 Metadata

gse <- getGEO("GSE41861", destdir = "C:/R/TH2_LOW/Metadata", getGPL = FALSE)

# 1.2.2 Extract metadata (pData)

metadata <- pData(phenoData(gse[[1]]))

#> pData of raw_data initially only contains the ID for the counts and an index
#> meanwhile metadata contains the full metadata but the IDs differ from the count data IDs

# Object that specifies which GSE is in use

gse_obj <- "GSE41861"

# 3.- Preprocessing metadata --------------------------------------------------

# 3.1 Object with the names of each file

file_name <- sampleNames(raw_data)

df <- as.data.frame(file_name)

#3.2 Clean metadata

# Assuming your data is in a dataframe called 'df' with column 'file_name'
# Use this regex to handle the inconsistent underscores and the "_2" suffix
df_cleaned <- df %>% 
  extract(file_name, 
          into = c("GSM", "Subject", "Middle_Info", "Atopy", "Tissue"),
          # Regex breakdown:
          # ([^_]+) : Match GSM
          # _([^_]+) : Match Subject (AA-XXX)
          # _(.*)_ : Match everything in the middle (Steroids/Severity/NS)
          # (POS|NEG) : Match Atopy
          # _([A-Z]+) : Match Tissue (NASAL/PRLL/PLLL)
          regex = "([^_]+)_([^_]+)_(.*)_(POS|NEG)_([A-Z]+)", 
          remove = FALSE) %>%
  # Now split the middle part (e.g., "OFF_MOD-PER" or just "NS")
  separate(Middle_Info, 
           into = c("Steroids", "Severity"), 
           sep = "_", 
           fill = "right")

df_cleaned$GSM[df_cleaned$Atopy == "NEG" & df_cleaned$Steroids != "NS"]

# View the result
head(df_cleaned)

metadata <- metadata %>%
  mutate(
    age = as.numeric(`age:ch1`),
    age.bin = ifelse(age <= median(age),
                     yes = "Young",
                     no = "Old"),
    age.bin = as.factor(age.bin),
    gender = ifelse(`Sex:ch1` == "M", "Male", "Female"),
    disease = `disease:ch1`,
    tissue = ifelse(`tissue:ch1` == "Upper airway (Nasal)", "Nasal", "Bronchial"),
    `age:ch1` = NULL,
    `disease:ch1` = NULL,
    `tissue:ch1` = NULL,
    `gender:ch1` = NULL,
    characteristics_ch1.1 = NULL,
    characteristics_ch1.2 = NULL,
    characteristics_ch1.3 = NULL,
    characteristics_ch1.4 = NULL,
    characteristics_ch1.1 = NULL
  ) %>% 
  rownames_to_column("GSM") %>% 
  left_join(df_cleaned, by = "GSM") %>% 
  column_to_rownames("GSM")



# 3.2.2 Add the object to the phenotipic data

pData(raw_data) <- pData(raw_data) %>% 
  rownames_to_column("file_name") %>%
  full_join(metadata, by = "file_name", keep = TRUE) %>% 
  column_to_rownames("file_name.y") %>% 
  rename(file_name = file_name.x)
  

metadata_all <- metadata

# 4.- Preprocess data -----------------------------------------------------

# 4.1 Boxplot of raw data

boxplot(raw_data)

# 4.2 Data normalization

norm_data <- rma(raw_data)

# 4.2.2 Boxplot of normalized data

boxplot(norm_data)

# 4.3 Expression matrix

expr_matrix <- exprs(norm_data)


# 5.- Probe ID to symbol --------------------------------------------------

# 5.1 Get mapping to change probe IDs to gene names

probe_gene <- AnnotationDbi::select(
  hgu133plus2.db,
  keys = rownames(expr_matrix),
  columns = "SYMBOL",
  keytype = "PROBEID"
)

# 5.2 Delete NA symbols and duplicate symbols and add to rownames

# 5.2.1 Transform matrix to a tidy data frame and add Symbols

gene_expres_matrix <- 
  expr_matrix %>%
  as.data.frame() %>%
  rownames_to_column("PROBEID") %>%
  inner_join(probe_gene, by = "PROBEID") %>%  # Join with your mapping object
  filter(!is.na(SYMBOL) & SYMBOL != "") %>%   # Remove NAs and empty symbols
  
  # 5.2.2 Calculate variance for each probe across all samples
  
  mutate(variance = apply(dplyr::select(., -PROBEID, -SYMBOL), 1, var)) %>%
  
  # 5.2.3 Keep only the probe with the highest variance per Gene Symbol
  
  group_by(SYMBOL) %>%
  slice_max(order_by = variance, n = 1, with_ties = FALSE) %>% 
  ungroup() %>%
  
  # 5.2.4 Remove columns, add symbol to rownames and reformat to matrix
  
  dplyr::select(-PROBEID, -variance) %>%
  column_to_rownames("SYMBOL") %>%
  as.matrix()




# 6.- Correct Batch effect ------------------------------------------------


library(sva)

# 6.1 Create objects for batch effect

pheno_all <- pData(norm_data) %>%
  mutate(Status = case_when(
    disease == "Control" ~ "Control",
    `severity:ch1` == "Mild" ~ "Asthma_Mild",
    `severity:ch1` == "Moderate" ~ "Asthma_Mod",
    `severity:ch1` == "Severe" ~ "Asthma_Sev",
    TRUE ~ "Other"
  ))

expr_all <- gene_expres_matrix

# 1. Create a single 'Status' variable that captures the hierarchy


# 6.2 Model

mod  <- model.matrix(~ as.factor(disease), data = pheno_all)
mod0 <- model.matrix(~ 1, data = pheno_all)

# 6.3 Find SVs across the whole dataset


expr_all <- as.matrix(expr_all)


# 6.3.2 Sv object

sv_obj <- sva(expr_all, mod, mod0)


########
# 6.3.2.2
# Object for later use with sv and patient IDs
sv_df <- as.data.frame(sv_obj$sv)
rownames(sv_df) <- colnames(gene_expres_matrix)
########

# 6.4 Clean the matrix for visualization/clustering. We use the known 'study' variable and the new SVs

expr_corrected <- removeBatchEffect(expr_all, 
                                    design = mod,
                                    covariates = sv_obj$sv)
  # 6.5 PCA CORRECTED

library(factoextra)


pca <- prcomp(t(expr_corrected))
eigen <- get_eigenvalue(pca)

pca_plot <- function(i, u, e) {
  df <- data.frame(
    PCi = pca$x[, i],
    PCu = pca$x[, u],
    disease = pData(norm_data)$disease,
    tissue = pData(norm_data)$Tissue,
    steroid = pData(norm_data)$Steroids,
    atopy = pData(norm_data)$Atopy
  )
  
  p <- ggplot(df, aes(PCi, PCu, color = .data[[e]])) +
    geom_point(size = 3)  +
    stat_ellipse(geom = "polygon", aes(fill = .data[[e]]), alpha = 0.2, level = 0.9) +
    labs(x = paste0("PC", i, " (", round(eigen[i, 2], 2), "%)"),
         y = paste0("PC", u, " (", round(eigen[u, 2], 2), "%)"),
         fill = as.character(e),
         color = as.character(e))
  print(p)
}

pca_plot(1, 2, e = "atopy")



# 6.5.2 Top contributors to PC1 and PC2

a <- fviz_contrib(pca, choice = "var", axes = 1, top = 30)


b <- fviz_contrib(pca, choice = "var", axes = 2, top = 30)


grid.arrange(a, b, ncol = 2)


# 7.- Creating special objects --------------------------------------------

# 7.0.1 Make sure that rownames of metadata and colnames of express matrix are in same order

all(rownames(pData(norm_data)) == colnames(gene_expres_matrix))

# 7.0.2 Object with asthma patients

asthma_patients <- pData(norm_data)$disease != "Control"

# 7.1 Data with only expression of asthma patients 

# 7.1.2 Object with count of asthma patients not batch corrected

expr_asthma <- gene_expres_matrix[, asthma_patients]

# 7.1.2.2 Object with count of asthma patients batch corrected

expr_asthma.batch <- expr_corrected[, asthma_patients]

# 7.2 Metadata with only asthma patients

metadata_asthma <- pData(norm_data)[asthma_patients, ]


# 8.- Special objects: Nasal samples ------------------------------------------------

nasal_sample.id <- pData(norm_data)$Tissue == "NASAL"


# 8.1 Data with only expression of nasal tissue

# 8.1.2 Object with count of nasal tissue not batch corrected

expr_nasal <- gene_expres_matrix[, nasal_sample.id]

# 8.1.2.2 Object with count of nasal tissue batch corrected

expr_nasal.batch <- expr_corrected[, nasal_sample.id]

# 8.2 Metadata with only nasal tissue

metadata_nasal <- pData(norm_data)[nasal_sample.id, ]



# 9.- Special objects: Non nasal samples ----------------------------------

bronchial_sample.id <- pData(norm_data)$Tissue != "NASAL"


# 9.1 Data with only expression of bronchial tissue

# 9.1.2 Object with count of bronchial tissue not batch corrected

expr_bronchial <- gene_expres_matrix[, bronchial_sample.id]

# 9.1.2.2 Object with count of bronchial tissue batch corrected

expr_bronchial.batch <- expr_corrected[, bronchial_sample.id]

# 9.2 Metadata with only bronchial tissue

metadata_bronchial <- pData(norm_data)[bronchial_sample.id, ]



# // Dictionary // --------------------------------------------------------
#> raw_data <- Data that corresponds to the extracted files, contains the expression matrix as well as 
#>            few other data. Unnormalized
#>            
#> metadata <- Object to save metadata if metadata_all has to be rerun
#> metadata_all <- Metadata with modifications for better lecture and correction of batch effect in this script
#>            and in diff expression
#>            
#> pData(norm_data) <-  Metadata in oligo object
#> 
#> expr_matrix <- Normalized expression matrix with PROBE IDs corresponding to multiple symbols
#> 
#> gene_expres_matrix <- Normalized expression matrix with unique SYMBOL IDs non batch proccessed
#> 
#> expr_corrected <- Expression data corrected by batch with SVA
#> 
#> expr_asthma <- Object with count of only asthma patients not batch corrected
#> 
#> expr_astma.batch <- Object with count of asthma patients batch corrected
#> 
#> metadata_asthma <- Metadata with only asthma patients
#> 
#> gse_obj <- Object that specifies which GSE is in use
