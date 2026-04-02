# Script to perform cluster analysis on asthma patients

library(ConsensusClusterPlus)


# 1.- Selecting top genes by variance -------------------------------------

# 1.1 Apply variance to the object of asthma only patients corrected by batch effect

var_genes <- apply(expr_asthma.batch, 1, var)

# 1.2 Sort in decreasing manner based on variance and select top 5000 genes

top_genes <- names(sort(var_genes, decreasing = TRUE))[1:10000]

# 1.3 Asign to object that will be the argument for the clustering

cluster_counts <- expr_asthma.batch[top_genes, ]

# 1.3.2 Convert to matrix

count_matrix <- as.matrix(cluster_counts)


# 2.- Clustering ----------------------------------------------------------

# 2.1 Clustering

results <- ConsensusClusterPlus(
  count_matrix, # Objetoa procesar, debe ser una matriz
  maxK = 5, # Numero maximo de clusters
  reps = 1000, # Cuantas veces se repetira la funcipon
  pItem = 0.8, # Cuantos items (1 = 100%)
  pFeature = 1, # Cuantos features (1 = 100%) 
  title = "Results/Cluster/", # A donde se van a guardar las imagenes
  clusterAlg = "pam", # Cluster jerarquico
  distance = "pearson",
  seed = 1262118388.71279,
  plot = "png"
)


# 2.2 Observe object

results

# 2.3 ntegretad classification likelyhood

icl = calcICL(results, title = "Results/Cluster/", plot = "png")
icl$clusterConsensus

# 2.3.2 Convert to data frame to aggregate

icl_df <- as.data.frame(icl$clusterConsensus)

# 2.3.3 Observe the highest mean for k in icl

aggregate(data = icl_df, clusterConsensus ~ k, FUN = mean)


# 2.4 Convert the results cluster object corresponding to selected k to a column in the metadata

metadata_asthma <- 
  as.data.frame(results[[2]]$consensusClass) %>% # Data frame corresponding to which cluster does each patient belong to
  rename(cluster = "results[[2]]$consensusClass") %>% # Rename column to "cluster"
  rownames_to_column("file_name") %>% # Asign rownames to a column for let join
  left_join(metadata_asthma, by = "file_name") %>% # Join
  column_to_rownames("file_name") %>% # Rename rownames
  mutate(cluster = as.factor(cluster))


# 3.- Observe clusters via PCA --------------------------------------------

# 3.1 PCA of expression

pca_asthma <- prcomp(t(expr_asthma.batch))

# 3.1.2 Eigenvalues of the PC

eigen_as <- get_eigenvalue(pca_asthma)

# 3.2 PCA plot function

pca_plot_as <- function(i, u){
  
  df2 <- data.frame(
    PCi = pca_asthma$x[,i],
    PCu = pca_asthma$x[,u],
    th2 = metadata_asthma$th2_group,
    cluster = factor(metadata_asthma$cluster),
    age = metadata_asthma$age.bin
  )
  
  p <- ggplot(df2, aes(PCi, PCu, color = th2, shape = age)) +
    geom_point(size = 3)  +
    stat_ellipse(geom = "polygon", 
                 aes(fill = cluster, group = cluster), 
                 alpha = 0.2, 
                 level = 0.9,
                 color = NA) +
    scale_fill_brewer(palette = "Set2") +
    labs(x = paste0("PC", i, " (", round(eigen_as[i, 2], 2), "%)"),
         y = paste0("PC", u, " (", round(eigen_as[u, 2], 2), "%)")
    )
  print(p)
}

# 3.2.2 Which PCA to plot

pca_plot_as(1, 2)


# 3.3 Top contributors to PC1 and PC2

a <- fviz_contrib(pca_asthma, choice = "var", axes = 1, top = 30)


b <- fviz_contrib(pca_asthma, choice = "var", axes = 2, top = 30)


grid.arrange(a, b, ncol = 2)


# 4.- Compare with initial publication ------------------------------------

# 4.1 Observe the distribution of the original publications TH2 classification
# within our clusters

metadata_asthma %>% 
  filter(.preserve = c(th2_group, cluster)) %>% 
  count(th2_group, cluster) %>% 
  group_by(cluster) %>% 
  mutate(prop = n / sum(n)) %>% 
  dplyr::select(-n) %>% 
  pivot_wider(names_from = th2_group, values_from = prop)


chisq.test(table(metadata_asthma$cluster, metadata_asthma$th2_group))


# // Dictionary // --------------------------------------------------------

#> metadata_asthma <- Metadata with only asthma patients and contains a column
#> with the cluster that each patient forms part of
