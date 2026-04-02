library(limma)
library(EnhancedVolcano)
library(org.Hs.eg.db)

# Object with patient IDs

metadata_asthma

# 2. Subset both to ONLY include the Asthma patients (excluding healthy)

sv_asthma <- sv_df[rownames(metadata_asthma), , drop = FALSE]

# 3. Rename SV columns safely

colnames(sv_asthma) <- paste0("SV", 1:ncol(sv_asthma))

# 4. Combine and create the design
# Use '0 +' so you get a column for each cluster (easier for contrasts)

col_data <- cbind(metadata_asthma, sv_asthma)

design <- model.matrix(~ 0 + cluster + age + gender + study, data = col_data)

# Clean column names (removes 'cluster' prefix from names)

colnames(design) <- make.names(colnames(design))


count_data <- expr_asthma

all(rownames(metadata_asthma) == colnames(expr_asthma))

# Fit the model

fit <- lmFit(count_data, design)

# Define your comparison between clusters
# Example: Cluster1 vs Cluster2
contrast.matrix <- makeContrasts(
  C1_vs_C2 = cluster1 - cluster2, 
  levels = design
)

fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)

# Get the differentially expressed genes
res <- topTable(fit, coef = "C1_vs_C2", number = Inf, sort.by = "P")

res

res_sig <- res %>%
  filter(adj.P.Val < 0.05 & abs(logFC) > 1)



qqt(fit$t[,1], df = fit$df.prior + fit$df.residual, 
    main = "Final Batch-Adjusted Q-Q Plot")
abline(0, 2, col = "red")



# We use the probe-level fit for the volcano to see all data points

EnhancedVolcano(res,
                lab = rownames(res),
                x = 'logFC',
                y = 'adj.P.Val',
                title = 'Cluster 2 vs Cluster 1',
                pCutoff = 0.05,
                FCcutoff = 1.0)


# // Dictionary // --------------------------------------------------------

#> res <-  object with Differential gene expression of all genes
#> 
#> res_sig <- Similar but only contains those that have an absolut log fold change > 1
#> and a p adjusted value of < 0.05

