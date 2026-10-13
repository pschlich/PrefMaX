# Code for running PrefMax on the Muscadine Grappe Juices example

library(PrefMaX)

path <- dirname(system.file("extdata", "Muscadine Grape Juices Data for PrefMaX.xlsx", package = "PrefMaX"))
file <- basename(system.file("extdata", "Muscadine Grape Juices Data for PrefMaX.xlsx", package = "PrefMaX"))

res <- run_prefmax(path = path,              # Pathname of the directory containing the Excel data file
                   file = file,              # Name of the Excel file containing data into 2 sheets
                   sheet_hedo = "Liking",    # Name of the data sheet containing the liking data"
                   sheet_qda = "Sensory",    # Name of the data sheet containg the sensory data
                   scale.unit= FALSE,        # FALSE for a PCA on coavariance rather than correlations
                   scale_desc = 2,           # Multiply default attribute coordinates by scale_desc
                   ncp=2,                    # Clustering on ncp first components
                   nb.clust = NULL,          # Number of clusters or NULL for automatic detection or 1 for no clustering
                   min_r2 = NULL,            # Threshold of qality of representation to display or not an attribute
                   label_consumers = NULL,   # No labels if FALSE or NULL and more than 100 consumers
                   output_dir= "C:/Pascal Schlich/PRO/Livre Tec&Doc Analyse Sensorielle/Muscadine")         # NULL for automatic detectiion of number of clusters or number of clusters (1 for no clustering)
