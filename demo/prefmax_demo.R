# Code for running PrefMax on the Muscadine Grape Juices example

library(PrefMaX)

res <- run_prefmax(path = "C:/Pascal Schlich/PRO/Thorsten Guksch",              # Pathname of the directory containing the Excel data file
                   file = "Data_PrefMaX.xlsx",              # Name of the Excel file containing data into 2 sheets
                   sheet_hedo = "Liking",    # Name of the data sheet containing the liking data
                   sheet_qda = "Sensory",    # Name of the data sheet containing the sensory data
                   scale.unit= FALSE,        # FALSE for a PCA on covariance rather than correlations
                   scale_desc = 5,           # Multiply default attribute coordinates by scale_desc
                   ncp=4,                    # Clustering on ncp first components
                   nb.clust = NULL,          # NULL for automatic detection of number of clusters or number of clusters (1 for no clustering)
                   min_r2 = 0.35,
                   label_consumers = NULL,
                   output_dir = "C:/Pascal Schlich/PRO/Thorsten Guksch"      # Store TIFF plots and Excel result file in working directory (default that can be changed)
                   )
