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
                   nb.clust = NULL)          # NULL for automatic detection of number of clusters or number of clusters (1 for no clustering)
                   output_dir = getwd()      # Store TIFF plots and Excel resul file in working directory (default that can be changed)
                   )
