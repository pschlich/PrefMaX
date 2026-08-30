# ==============================================================================
# 1. Calcul de la méthode PrefMaX
# ==============================================================================

#' Calcul de la méthode PrefMaX
#' @param Y Matrice ou data.frame (n x p) des notes hédoniques (n produits en lignes, p consommateurs en colonnes).
#' @param X Matrice ou data.frame (n x q) des moyennes des descripteurs sensoriels (n produits en lignes, q descripteurs en colonnes).
#' @param scale.unit Logique. Si TRUE, l'ACP est normée. Si FALSE (défaut), l'ACP est non normée.
#' @param ncp Nombre de dimensions conservées dans l'ACP (défaut 5).
#'
#' @return Un objet de classe \code{prefmax} contenant une liste avec Y, X, R, XOPT, bounds, weights, et pca.
#' @export
prefmax <- function(Y, X, scale.unit = FALSE, ncp = 5) {
  Y <- as.matrix(Y)
  X <- as.matrix(X)

  if (nrow(X) != nrow(Y)) {
    stop("Y et X doivent avoir le même nombre de lignes (produits).")
  }

  # Alignement automatique des produits par rownames si présents
  if (!is.null(rownames(Y)) && !is.null(rownames(X))) {
    common_prods <- intersect(rownames(X), rownames(Y))
    if (length(common_prods) != nrow(X)) {
      stop("Les noms de produits dans Y et X ne correspondent pas exactement.")
    }
    Y <- Y[common_prods, , drop = FALSE]
    X <- X[common_prods, , drop = FALSE]
  }

  n <- nrow(Y)
  p <- ncol(Y)
  q <- ncol(X)

  if (is.null(rownames(Y))) rownames(Y) <- paste0("P", 1:n)
  if (is.null(colnames(Y))) colnames(Y) <- paste0("Cons_", 1:p)
  if (is.null(colnames(X))) colnames(X) <- paste0("Desc_", 1:q)

  R_mat <- matrix(NA, nrow = p, ncol = q, dimnames = list(colnames(Y), colnames(X)))
  XOPT_mat <- matrix(NA, nrow = p, ncol = q, dimnames = list(colnames(Y), colnames(X)))

  bounds <- data.frame(
    Min = apply(X, 2, min, na.rm = TRUE),
    Max = apply(X, 2, max, na.rm = TRUE),
    row.names = colnames(X)
  )

  for (k in 1:q) {
    x_k <- X[, k]
    min_x <- bounds$Min[k]
    max_x <- bounds$Max[k]

    for (j in 1:p) {
      y_j <- Y[, j]

      fit <- lm(y_j ~ x_k + I(x_k^2))

      r2 <- summary(fit)$r.squared
      R_mat[j, k] <- sqrt(r2)

      coeffs <- coef(fit)
      b <- ifelse(is.na(coeffs["x_k"]), 0, coeffs["x_k"])
      c <- ifelse(is.na(coeffs["I(x_k^2)"]), 0, coeffs["I(x_k^2)"])

      pred_fun <- function(x) predict(fit, newdata = data.frame(x_k = x))
      candidates <- c(min_x, max_x)

      if (c < 0) {
        x_top <- -b / (2 * c)
        if (x_top >= min_x && x_top <= max_x) {
          candidates <- c(candidates, x_top)
        }
      }

      preds <- pred_fun(candidates)
      best_idx <- which.max(preds)
      XOPT_mat[j, k] <- candidates[best_idx]
    }
  }

  weights <- rowMeans(R_mat)

  pca_res <- FactoMineR::PCA(
    XOPT_mat,
    scale.unit = scale.unit,
    ncp = ncp,
    row.w = weights,
    graph = FALSE
  )

  res <- list(
    Y = Y,
    X = X,
    R = R_mat,
    XOPT = XOPT_mat,
    bounds = bounds,
    weights = weights,
    pca = pca_res,
    scale.unit = scale.unit
  )

  class(res) <- "prefmax"
  return(res)
}

# ==============================================================================
# 2. Graphique de l'éboulis des valeurs propres (Scree Plot) pour PrefMaX
# ==============================================================================

#' Graphique de l'éboulis des valeurs propres (Scree Plot) pour PrefMaX
#'
#' @param x Un objet de classe \code{prefmax}.
#' @param title Titre du graphique.
#'
#' @return Un objet \code{ggplot} et affiche le tableau dans la console.
#' @export
screeplot_prefmax <- function(x, title = "PrefMaX - Scree Plot (Eigenvalues)") {
  if (!inherits(x, "prefmax")) {
    stop("L'objet doit être de classe 'prefmax'.")
  }

  eig <- as.data.frame(x$pca$eig)
  colnames(eig) <- c("Eigenvalue", "Variance_Pct", "Cumulative_Pct")

  cat("\n=======================================================\n")
  cat("                EIGENVALUES & VARIANCE                 \n")
  cat("=======================================================\n")
  print(round(eig, 2))
  cat("=======================================================\n\n")

  eig$Axis <- factor(1:nrow(eig))

  # Facteur de conversion pour aligner le % cumulé (0-100%) sur la hauteur des barres (%)
  max_var <- max(eig$Variance_Pct, na.rm = TRUE)
  scale_factor <- max_var / 100

  p <- ggplot2::ggplot(eig, ggplot2::aes(x = Axis)) +
    # Barres bleues pour le % de variance par axe
    ggplot2::geom_bar(ggplot2::aes(y = Variance_Pct), stat = "identity",
                      fill = "steelblue", color = "black", alpha = 0.7, width = 0.6) +
    ggplot2::geom_text(ggplot2::aes(y = Variance_Pct, label = paste0(round(Variance_Pct, 1), "%")),
                       vjust = -0.5, size = 3, fontface = "bold", color = "steelblue") +
    # Courbe rouge pour le % cumulé (rapporté à l'échelle secondaire)
    ggplot2::geom_line(ggplot2::aes(y = Cumulative_Pct * scale_factor, group = 1),
                       color = "firebrick", linewidth = 1) +
    ggplot2::geom_point(ggplot2::aes(y = Cumulative_Pct * scale_factor),
                        color = "firebrick", size = 2.5) +
    ggplot2::geom_text(ggplot2::aes(y = Cumulative_Pct * scale_factor, label = paste0(round(Cumulative_Pct, 1), "%")),
                       vjust = -0.8, size = 2.8, fontface = "bold", color = "firebrick") +
    # Axes principal (gauche) et secondaire (droite)
    ggplot2::scale_y_continuous(
      name = "% of Variance",
      limits = c(0, max_var * 1.18),
      sec.axis = ggplot2::sec_axis(~ . / scale_factor, name = "Cumulative Variance (%)")
    ) +
    ggplot2::labs(
      title = title,
      x = "Dimensions / Axes"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 13),
      axis.title.y.left = ggplot2::element_text(color = "steelblue", face = "bold"),
      axis.title.y.right = ggplot2::element_text(color = "firebrick", face = "bold"),
      axis.text.y.left = ggplot2::element_text(color = "steelblue"),
      axis.text.y.right = ggplot2::element_text(color = "firebrick"),
      panel.border = ggplot2::element_rect(color = "gray80", fill = NA, linewidth = 0.5)
    )

  return(p)
}

# ==============================================================================
# 3. Représentation graphique Biplot de la méthode PrefMaX
# ==============================================================================

#' Représentation graphique Biplot de la méthode PrefMaX
#'
#' @param x Un objet de classe \code{prefmax}.
#' @param axes Un vecteur de deux entiers indiquant les axes de l'ACP à afficher (défaut c(1, 2)).
#' @param scale_desc Facteur de dilatation pour les coordonnées des descripteurs (défaut 1).
#' @param cluster Objet de classe \code{prefmax_cluster}. Si NULL, affichage standard.
#' @param title Le titre du graphique.
#' @param label_size Taille de la police pour les étiquettes.
#' @param max.overlaps Nombre maximal de chevauchements pour ggrepel (défaut Inf).
#' @param step Espacement souhaité entre les lignes de la grille.
#' @param ... Arguments supplémentaires.
#'
#' @return Un objet \code{ggplot} représentant le biplot PrefMaX.
#' @export
plot.prefmax <- function(x, axes = c(1, 2),
                         scale_desc = 1,
                         cluster = NULL,
                         title = "PrefMaX",
                         label_size = 3.5,
                         max.overlaps = Inf,
                         step = NULL, ...) {

  if (!inherits(x, "prefmax")) {
    stop("L'objet doit être de classe 'prefmax'.")
  }

  ax1 <- axes[1]
  ax2 <- axes[2]
  pca_res <- x$pca

  cons_coords <- as.data.frame(pca_res$ind$coord[, c(ax1, ax2)])
  colnames(cons_coords) <- c("Dim1", "Dim2")
  cons_coords$Label <- rownames(cons_coords)
  cons_coords$Fontface <- "plain"
  cons_coords$Size <- label_size

  if (!is.null(cluster)) {
    if (!inherits(cluster, "prefmax_cluster")) {
      stop("L'argument 'cluster' doit être un objet de classe 'prefmax_cluster'.")
    }
    cons_coords$Type <- as.character(cluster$clusters[rownames(cons_coords)])
  } else {
    cons_coords$Type <- "Consumers"
  }

  desc_coords <- as.data.frame(pca_res$var$coord[, c(ax1, ax2)] * scale_desc)
  colnames(desc_coords) <- c("Dim1", "Dim2")
  desc_coords$Label <- rownames(desc_coords)
  desc_coords$Type <- "Attributes"
  desc_coords$Fontface <- "bold"
  desc_coords$Size <- label_size + 0.5

  center_xopt <- pca_res$call$centre
  X_centered <- scale(x$X, center = center_xopt, scale = FALSE)
  V <- pca_res$svd$V[, c(ax1, ax2)]
  prod_coords <- as.data.frame(X_centered %*% V)
  colnames(prod_coords) <- c("Dim1", "Dim2")
  prod_coords$Label <- rownames(prod_coords)
  prod_coords$Type <- "Products"
  prod_coords$Fontface <- "bold"
  prod_coords$Size <- label_size + 0.5

  all_labels <- rbind(cons_coords, desc_coords, prod_coords)

  eig <- pca_res$eig
  var_ax1 <- round(eig[ax1, 2], 1)
  var_ax2 <- round(eig[ax2, 2], 1)

  all_x <- all_labels$Dim1
  all_y <- all_labels$Dim2
  min_val <- min(c(all_x, all_y), na.rm = TRUE)
  max_val <- max(c(all_x, all_y), na.rm = TRUE)

  if (is.null(step)) {
    range_val <- max_val - min_val
    if (range_val <= 4) step <- 0.5 else if (range_val <= 10) step <- 1 else step <- 2
  }
  breaks_seq <- seq(floor(min_val / step) * step, ceiling(max_val / step) * step, by = step)

  if (!is.null(cluster)) {
    k <- cluster$nb.clust
    clust_colors <- hcl(h = seq(15, 375, length = k + 1), l = 45, c = 100)[1:k]
    names(clust_colors) <- paste0("Cluster_", 1:k)

    color_palette <- c(
      "Products"   = "red",
      "Attributes" = "forestgreen",
      clust_colors
    )
  } else {
    color_palette <- c(
      "Products"   = "red",
      "Consumers"  = "blue",
      "Attributes" = "forestgreen"
    )
  }

  p <- ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
    ggplot2::geom_point(data = cons_coords, ggplot2::aes(x = Dim1, y = Dim2, color = Type), alpha = 0.7, size = 2) +
    ggplot2::geom_segment(data = desc_coords, ggplot2::aes(x = 0, y = 0, xend = Dim1, yend = Dim2, color = Type),
                          arrow = ggplot2::arrow(length = ggplot2::unit(0.2, "cm")), linewidth = 0.7) +
    ggplot2::geom_point(data = prod_coords, ggplot2::aes(x = Dim1, y = Dim2, color = Type), size = 3.5, shape = 17)

  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p + ggrepel::geom_text_repel(
      data = all_labels,
      ggplot2::aes(x = Dim1, y = Dim2, label = Label, color = Type, fontface = Fontface, size = Size),
      max.overlaps = max.overlaps,
      box.padding = 0.35,
      point.padding = 0.3,
      show.legend = FALSE
    ) + ggplot2::scale_size_identity()
  } else {
    p <- p + ggplot2::geom_text(
      data = all_labels,
      ggplot2::aes(x = Dim1, y = Dim2, label = Label, color = Type, fontface = Fontface),
      vjust = -0.7, size = label_size, show.legend = FALSE
    )
  }

  p <- p +
    ggplot2::scale_color_manual(values = color_palette, name = NULL) +
    ggplot2::scale_x_continuous(breaks = breaks_seq) +
    ggplot2::scale_y_continuous(breaks = breaks_seq) +
    ggplot2::coord_fixed(ratio = 1) +
    ggplot2::labs(
      title = title,
      x = paste0("Axis ", ax1, " (", var_ax1, " %)"),
      y = paste0("Axis ", ax2, " (", var_ax2, " %)")
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(title = NULL),
      shape = ggplot2::guide_legend(title = NULL)
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14),
      legend.position = "bottom",
      panel.border = ggplot2::element_rect(color = "gray80", fill = NA, linewidth = 0.5),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA)
    )

  return(p)
}

# ==============================================================================
# 4. Classification des consommateurs dans l'espace PrefMaX
# ==============================================================================

#' Classification des consommateurs dans l'espace PrefMaX
#'
#' @param x Un objet de classe \code{prefmax}.
#' @param ncp Nombre de dimensions de l'ACP à prendre en compte (défaut 2).
#' @param nb.clust Nombre de clusters souhaité. Si NULL, déterminé par silhouette.
#' @param max.clust Nombre maximal de clusters à tester (défaut 6).
#'
#' @return Un objet de classe \code{prefmax_cluster}.
#' @export
prefmax_cluster <- function(x, ncp = 2, nb.clust = NULL, max.clust = 6) {
  if (!inherits(x, "prefmax")) {
    stop("L'objet doit être de classe 'prefmax'.")
  }

  max_available_ncp <- ncol(x$pca$ind$coord)
  ncp_used <- min(ncp, max_available_ncp)
  coords <- x$pca$ind$coord[, 1:ncp_used, drop = FALSE]

  p <- nrow(coords)
  q <- ncol(x$XOPT)

  sil_df <- NULL

  if (is.null(nb.clust)) {
    if (!requireNamespace("cluster", quietly = TRUE)) {
      stop("Le package 'cluster' est requis. Veuillez l'installer ou spécifier 'nb.clust'.")
    }

    sil_scores <- numeric(max.clust - 1)
    for (k in 2:max.clust) {
      km_test <- kmeans(coords, centers = k, nstart = 25)
      sil <- cluster::silhouette(km_test$cluster, dist(coords))
      sil_scores[k - 1] <- mean(sil[, 3])
    }
    nb.clust <- which.max(sil_scores) + 1
    message(paste0("Nombre de clusters déterminé automatiquement par la silhouette (sur ", ncp_used, " axes) : K = ", nb.clust))

    sil_df <- data.frame(
      Clusters = 2:max.clust,
      Average_Silhouette = round(sil_scores, 4),
      Selected = ifelse(2:max.clust == nb.clust, "YES", "")
    )
  }

  set.seed(123)
  km_res <- kmeans(coords, centers = nb.clust, nstart = 50)

  clusters <- factor(km_res$cluster, labels = paste0("Cluster_", 1:nb.clust))
  names(clusters) <- rownames(coords)

  XOPT_mat <- x$XOPT
  R_mat <- x$R

  optima_by_cluster <- matrix(NA, nrow = nb.clust, ncol = q,
                              dimnames = list(paste0("Cluster_", 1:nb.clust), colnames(XOPT_mat)))

  for (c_idx in 1:nb.clust) {
    clust_name <- paste0("Cluster_", c_idx)
    cons_in_clust <- names(clusters)[clusters == clust_name]

    for (k in 1:q) {
      xopt_vals <- XOPT_mat[cons_in_clust, k]
      r_weights <- R_mat[cons_in_clust, k]

      if (sum(r_weights, na.rm = TRUE) > 0) {
        optima_by_cluster[c_idx, k] <- sum(xopt_vals * r_weights, na.rm = TRUE) / sum(r_weights, na.rm = TRUE)
      } else {
        optima_by_cluster[c_idx, k] <- mean(xopt_vals, na.rm = TRUE)
      }
    }
  }

  res <- list(
    clusters = clusters,
    nb.clust = nb.clust,
    ncp = ncp_used,
    centers = km_res$centers,
    optima_means = as.data.frame(optima_by_cluster),
    silhouette_table = sil_df,
    prefmax_obj = x
  )

  class(res) <- "prefmax_cluster"
  return(res)
}

# ==============================================================================
# 5. Caractérisation des clusters par ANOVA pondérée
# ==============================================================================

#' Caractérisation des clusters de consommateurs par ANOVA pondérée
#'
#' @param object Un objet de classe \code{prefmax_cluster}.
#' @param ... Arguments supplémentaires.
#'
#' @return Data.frame avec ANOVA et optima moyens par cluster.
#' @export
summary.prefmax_cluster <- function(object, ...) {
  if (!inherits(object, "prefmax_cluster")) {
    stop("L'objet doit être de classe 'prefmax_cluster'.")
  }

  xopt_mat <- object$prefmax_obj$XOPT
  r_mat    <- object$prefmax_obj$R
  clusters <- object$clusters
  nb_clust <- object$nb.clust
  q        <- ncol(xopt_mat)
  desc_names <- colnames(xopt_mat)

  cluster_cols <- paste0("Cluster_", 1:nb_clust)
  res_df <- data.frame(
    Attribute = desc_names,
    F_stat = numeric(q),
    probF = numeric(q),
    matrix(NA, nrow = q, ncol = nb_clust, dimnames = list(NULL, cluster_cols)),
    stringsAsFactors = FALSE
  )

  for (k in 1:q) {
    var_name <- desc_names[k]

    df_k <- data.frame(
      XOPT = xopt_mat[, k],
      Weight = r_mat[, k],
      Cluster = clusters
    )

    df_k <- df_k[!is.na(df_k$XOPT) & !is.na(df_k$Weight) & df_k$Weight > 0, ]

    if (nrow(df_k) > nb_clust && length(unique(df_k$Cluster)) > 1) {
      fit <- lm(XOPT ~ Cluster, data = df_k, weights = Weight)
      anova_res <- anova(fit)

      res_df$F_stat[k] <- anova_res["Cluster", "F value"]
      res_df$probF[k]  <- anova_res["Cluster", "Pr(>F)"]
    } else {
      res_df$F_stat[k] <- NA
      res_df$probF[k]  <- NA
    }

    for (c_idx in 1:nb_clust) {
      c_name <- paste0("Cluster_", c_idx)
      res_df[k, c_name] <- object$optima_means[c_name, var_name]
    }
  }

  res_df <- res_df[order(res_df$F_stat, decreasing = TRUE, na.last = TRUE), ]

  res_df$F_stat <- round(res_df$F_stat, 2)
  res_df$probF  <- format.pval(res_df$probF, digits = 3, eps = 0.001)
  for (col in cluster_cols) {
    res_df[[col]] <- round(res_df[[col]], 2)
  }

  rownames(res_df) <- NULL
  return(res_df)
}

# ==============================================================================
# 6. Pipeline globale run_prefmax
# ==============================================================================

#' Key-in-hand pipeline for the PrefMaX method
#'
#' @param path Directory containing the input Excel data file.
#' @param file Excel file name.
#' @param sheet_hedo Sheet name or index for liking data (default 1).
#' @param sheet_qda Sheet name or index for QDA data (default 2).
#' @param scale.unit Logical. If TRUE, standardized PCA. If FALSE (default), unstandardized PCA.
#' @param ncp Number of dimensions to retain for PCA and clustering (default 2).
#' @param nb.clust Number of clusters (if NULL, auto silhouette; if 1, no clustering).
#' @param scale_desc Scaling factor for attribute arrows (default 1).
#' @param title_standard Base title for standard biplot.
#' @param title_cluster Base title for clustered biplot.
#' @param res_dpi Resolution of saved TIFF images (default 300 DPI).
#' @param output_dir Directory where output files (Excel and TIFF) will be saved. Default is current working directory \code{getwd()}. Set to \code{NULL} to disable file saving.
#'
#' @return An object of class \code{prefmax_pipeline}.
#' @export
#'
#' @examples
#' library(PrefMaX)
#'
#' path <- dirname(system.file("extdata", "Muscadine Grape Juices Data for PrefMaX.xlsx", package = "PrefMaX"))
#' file <- basename(system.file("extdata", "Muscadine Grape Juices Data for PrefMaX.xlsx", package = "PrefMaX"))
#'
#' res <- run_prefmax(
#'   path = path,              # Directory containing the example Excel file
#'   file = file,              # Name of the example Excel file
#'   sheet_hedo = "Liking",    # Name of liking sheet
#'   sheet_qda = "Sensory",    # Name of sensory sheet
#'   scale.unit = FALSE,       # FALSE for PCA on covariance
#'   scale_desc = 2,           # Multiply default attribute coordinates by 2
#'   ncp = 2,                  # Clustering based on the two first PCs
#'   nb.clust = NULL,          # NULL for automatic cluster detection
#'   output_dir = getwd()      # Save generated files in current working directory
#' )

run_prefmax <- function(path,
                        file,
                        sheet_hedo = 1,
                        sheet_qda = 2,
                        scale.unit = FALSE,
                        ncp = 2,
                        nb.clust = NULL,
                        scale_desc = 1,
                        title_standard = "PrefMaX Standard Biplot",
                        title_cluster = "PrefMaX Clustered Biplot",
                        res_dpi = 300,
                        output_dir = getwd()) {

  if (!dir.exists(path)) {
    stop("The specified input directory does not exist: ", path)
  }

  if (!is.null(output_dir) && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  if (!requireNamespace("readxl", quietly = TRUE) || !requireNamespace("openxlsx", quietly = TRUE)) {
    stop("Packages 'readxl' and 'openxlsx' are required.")
  }

  full_file_path <- file.path(path, file)
  if (!file.exists(full_file_path)) {
    stop("Excel file '", file, "' not found in directory: ", path)
  }

  message("--> Loading Excel data...")
  Y_df <- readxl::read_excel(full_file_path, sheet = sheet_hedo)
  X_df <- readxl::read_excel(full_file_path, sheet = sheet_qda)

  Y <- as.matrix(Y_df[, -1])
  rownames(Y) <- as.character(as.data.frame(Y_df)[, 1])

  X <- as.matrix(X_df[, -1])
  rownames(X) <- as.character(as.data.frame(X_df)[, 1])

  n_products   <- nrow(X)
  n_consumers  <- ncol(Y)
  n_attributes <- ncol(X)
  prod_names   <- rownames(X)
  attr_names   <- colnames(X)

  cat("\n=======================================================\n")
  cat("               DATA VERIFICATION SUMMARY               \n")
  cat("=======================================================\n")
  cat("• Products     (", n_products, ") : ", paste(prod_names, collapse = ", "), "\n", sep = "")
  cat("• Consumers    (", n_consumers, ")\n", sep = "")
  cat("• Attributes   (", n_attributes, ") : ", paste(attr_names, collapse = ", "), "\n", sep = "")
  cat("• PCA Type     : ", ifelse(scale.unit, "Standardized (scale.unit = TRUE)", "Unstandardized (scale.unit = FALSE)"), "\n", sep = "")
  cat("• Dimensions   : ncp = ", ncp, "\n", sep = "")
  cat("=======================================================\n\n")

  message("--> Running PrefMaX analysis...")
  res_pref <- prefmax(Y, X, scale.unit = scale.unit, ncp = ncp)

  p_scree <- screeplot_prefmax(res_pref)
  print(p_scree) # Affichage immédiat du screeplot à l'écran

  if (!is.null(output_dir)) {
    ggplot2::ggsave(filename = file.path(output_dir, "prefmax_screeplot.tiff"), plot = p_scree, device = "tiff",
                    width = 8, height = 5, units = "in", dpi = res_dpi, compression = "lzw", bg = "white")
  }

  wb <- openxlsx::createWorkbook()

  # --------------------------------------------------------------------------
  # Préparation des DataFrames pour les 4 feuilles de base
  # --------------------------------------------------------------------------
  df_data_verif <- data.frame(
    Parametre = c("Repertoire (PATH)", "Fichier source", "Nombre de Produits (N)", "Nombre de Consommateurs (P)", "Nombre de Descripteurs (Q)", "Type ACP", "Dimensions (ncp)"),
    Valeur = c(path, file, n_products, n_consumers, n_attributes,
               ifelse(scale.unit, "Standardisee (scale.unit = TRUE)", "Non standardisee (scale.unit = FALSE)"),
               ncp),
    stringsAsFactors = FALSE
  )

  df_eigenvalues <- as.data.frame(res_pref$pca$eig)
  colnames(df_eigenvalues) <- c("Eigenvalue", "Variance_Pct", "Cumulative_Pct")
  df_eigenvalues <- round(df_eigenvalues, 4)

  df_xopt <- as.data.frame(res_pref$XOPT)

  df_r_sqrt <- as.data.frame(res_pref$R)
  df_r_sqrt$`Subject Weight` <- res_pref$weights

  # Écriture des 4 feuilles de base
  openxlsx::addWorksheet(wb, "Data_Verification")
  openxlsx::writeData(wb, "Data_Verification", df_data_verif, startRow = 1, startCol = 1)
  openxlsx::writeData(wb, "Data_Verification", "Domaine de variation des descripteurs (Min / Max) :", startRow = 10, startCol = 1)
  openxlsx::writeData(wb, "Data_Verification", res_pref$bounds, startRow = 11, startCol = 1, rowNames = TRUE)

  openxlsx::addWorksheet(wb, "Eigenvalues")
  openxlsx::writeData(wb, "Eigenvalues", df_eigenvalues, rowNames = TRUE)

  openxlsx::addWorksheet(wb, "XOPT")
  openxlsx::writeData(wb, "XOPT", df_xopt, rowNames = TRUE)

  openxlsx::addWorksheet(wb, "R_sqrt")
  openxlsx::writeData(wb, "R_sqrt", df_r_sqrt, rowNames = TRUE)

  plots_standard_list <- list()
  plots_cluster_list  <- list()

  # --------------------------------------------------------------------------
  # CAS 1 : nb.clust = 1 (Pas de clustering)
  # --------------------------------------------------------------------------
  if (!is.null(nb.clust) && nb.clust == 1) {
    message("--> Option nb.clust = 1 detected: Global analysis without clustering.")

    res_clust <- NULL
    table_res <- NULL

    if (!is.null(output_dir)) {
      excel_out_name <- paste0("prefmax_global_results_", format(Sys.Date(), "%Y%m%d"), ".xlsx")
      openxlsx::saveWorkbook(wb, file = file.path(output_dir, excel_out_name), overwrite = TRUE)
    }

    # Génération uniquement des Biplots Standard
    for (i in 1:(ncp - 1)) {
      ax_pair <- c(i, i + 1)
      pair_name <- paste0("axes_", i, "_", i + 1)
      t_title <- paste0(title_standard, " (Axes ", i, "-", i + 1, ")")

      p_std <- plot(res_pref, axes = ax_pair, scale_desc = scale_desc, title = t_title)
      plots_standard_list[[pair_name]] <- p_std
      print(p_std)

      if (!is.null(output_dir)) {
        file_std <- paste0("prefmax_biplot_standard_axes_", i, "_", i + 1, ".tiff")
        ggplot2::ggsave(filename = file.path(output_dir, file_std), plot = p_std, device = "tiff",
                        width = 10, height = 7, units = "in", dpi = res_dpi, compression = "lzw", bg = "white")
      }
    }

    # --------------------------------------------------------------------------
    # CAS 2 : Clustering des consommateurs (nb.clust NULL ou > 1)
    # --------------------------------------------------------------------------
  } else {
    message("--> Performing consumer clustering on ", ncp, " dimensions...")
    res_clust <- prefmax_cluster(res_pref, ncp = ncp, nb.clust = nb.clust)

    message("--> Characterizing clusters via weighted ANOVA...")
    table_res <- summary(res_clust)

    df_silhouette <- res_clust$silhouette_table
    if (is.null(df_silhouette)) {
      df_silhouette <- data.frame(
        Note = paste0("Clustering performed with user-defined K = ", res_clust$nb.clust, " on ", ncp, " dimensions")
      )
    }

    df_consumers <- data.frame(
      Consumer = colnames(Y),
      Cluster = as.character(res_clust$clusters),
      stringsAsFactors = FALSE
    )

    # Ajout des 3 feuilles de clustering
    openxlsx::addWorksheet(wb, "Silhouette_Scores")
    openxlsx::writeData(wb, "Silhouette_Scores", df_silhouette)

    openxlsx::addWorksheet(wb, "Consumer_Clusters")
    openxlsx::writeData(wb, "Consumer_Clusters", df_consumers)

    openxlsx::addWorksheet(wb, "Cluster_Characterization")
    openxlsx::writeData(wb, "Cluster_Characterization", table_res)

    if (!is.null(output_dir)) {
      excel_out_name <- paste0("prefmax_cluster_results_", format(Sys.Date(), "%Y%m%d"), ".xlsx")
      openxlsx::saveWorkbook(wb, file = file.path(output_dir, excel_out_name), overwrite = TRUE)
    }

    # Génération des Biplots Standard ET des Biplots avec Clusters
    for (i in 1:(ncp - 1)) {
      ax_pair <- c(i, i + 1)
      pair_name <- paste0("axes_", i, "_", i + 1)

      t_title_std <- paste0(title_standard, " (Axes ", i, "-", i + 1, ")")
      p_std <- plot(res_pref, axes = ax_pair, scale_desc = scale_desc, title = t_title_std)
      plots_standard_list[[pair_name]] <- p_std
      print(p_std)

      if (!is.null(output_dir)) {
        file_std <- paste0("prefmax_biplot_standard_axes_", i, "_", i + 1, ".tiff")
        ggplot2::ggsave(filename = file.path(output_dir, file_std), plot = p_std, device = "tiff",
                        width = 10, height = 7, units = "in", dpi = res_dpi, compression = "lzw", bg = "white")
      }

      t_title_clust <- paste0(title_cluster, " (Axes ", i, "-", i + 1, ")")
      p_clust <- plot(res_pref, axes = ax_pair, cluster = res_clust, scale_desc = scale_desc, title = t_title_clust)
      plots_cluster_list[[pair_name]] <- p_clust
      print(p_clust)

      if (!is.null(output_dir)) {
        file_clust <- paste0("prefmax_biplot_clusters_axes_", i, "_", i + 1, ".tiff")
        ggplot2::ggsave(filename = file.path(output_dir, file_clust), plot = p_clust, device = "tiff",
                        width = 10, height = 7, units = "in", dpi = res_dpi, compression = "lzw", bg = "white")
      }
    }
  }

  if (!is.null(output_dir)) {
    message("\n==> Execution completed successfully! Excel workbook and TIFF plots saved in: ", output_dir)
  } else {
    message("\n==> Execution completed successfully! (No files saved on disk)")
  }

  out <- list(
    prefmax = res_pref,
    cluster = res_clust,
    summary_table = table_res,
    plots_standard = plots_standard_list,
    plots_cluster = plots_cluster_list
  )
  class(out) <- "prefmax_pipeline"
  return(invisible(out))
}
