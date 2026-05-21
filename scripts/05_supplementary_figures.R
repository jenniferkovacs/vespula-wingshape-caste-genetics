# =============================================================================
# 05_supplementary_figures.R
# Vespula maculifrons Wing Shape Study
# Murray, Orr, Goodisman, Gosse, Hill, Kovacs
#
# PURPOSE: Generate supplementary display items:
#
#   Supplementary Figure 2: Deformation grids showing mean shape differences
#     between colonies within each caste (20 grids: 10 colonies x 2 pairwise
#     directions, repeated for gynes, workers, and males).
#
#   Supplementary Figure 3: Deformation grids showing mean shape differences
#     between patrilines within colonies within each female caste (up to 60
#     grids across all colony x patriline combinations for gynes and workers).
#
#   Also exports landmark displacement tables (CSV) for all pairwise
#   comparisons at the caste, colony, and patriline levels.
#
# INPUT:   output/01_gpa_workspace.RData (Y.gpa, gdf, PCA, variance_props)
#
# OUTPUT:  Deformation grid PDFs
#          CSV displacement tables:
#            disp_array_caste_wide.csv
#            worker_colony_disp.csv
#            gyne_colony_disp.csv
#            male_colony_disp.csv
#            all_worker_patriline_colonies_disp.csv
#            all_gyne_patriline_colonies_disp.csv
#
# DEPENDENCIES: 00_setup.R, 01_data_import_and_GPA.R
# =============================================================================

source("scripts/00_setup.R")
load("output/01_gpa_workspace.RData")

# Helper: compute per-landmark displacement magnitude between two shapes
compute_disp <- function(shapeA, shapeB) sqrt(rowSums((shapeA - shapeB)^2))

# Landmark links for deformation grids
ref   <- mshape(Y.gpa$coords)
links <- define.links(ref)

# =============================================================================
# SUPPLEMENTARY TABLE: Caste-level landmark displacement
# =============================================================================

# Compute mean shapes per caste
means_by_caste <- list(
  gyne   = mshape(Y.gpa$coords[,, Y.gpa$caste == "gyne"]),
  worker = mshape(Y.gpa$coords[,, Y.gpa$caste == "worker"]),
  male   = mshape(Y.gpa$coords[,, Y.gpa$caste == "male"])
)

castes  <- names(means_by_caste)
num_lms <- nrow(means_by_caste[[1]])

# 3D array: landmarks x FromCaste x ToCaste
disp_array_caste <- array(
  NA,
  dim      = c(num_lms, length(castes), length(castes)),
  dimnames = list(Landmark  = paste0("LM", 1:num_lms),
                  FromCaste = castes,
                  ToCaste   = castes)
)

for (i in seq_along(castes)) {
  for (j in seq_along(castes)) {
    if (i != j) {
      disp_array_caste[, i, j] <- compute_disp(means_by_caste[[castes[i]]],
                                                means_by_caste[[castes[j]]])
    }
  }
}

# Reshape to long then wide format for export
disp_df_caste <- melt(disp_array_caste,
                      varnames = c("Landmark", "FromCaste", "ToCaste"))
names(disp_df_caste)[4] <- "Displacement"
disp_df_caste <- disp_df_caste[!is.na(disp_df_caste$Displacement), ]

disp_wide_caste <- disp_df_caste %>%
  unite(Comparison, FromCaste, ToCaste, sep = "_vs_") %>%
  pivot_wider(names_from = Landmark, values_from = Displacement)

write.csv(disp_wide_caste, "output/disp_array_caste_wide.csv", row.names = FALSE)
cat("Caste displacement table saved.\n")

# =============================================================================
# SUPPLEMENTARY FIGURE 2: Colony-level deformation grids within each caste
# =============================================================================

colony_disp_wide_list <- list()

for (caste_name in c("gyne", "worker", "male")) {

  caste_filter       <- gdf$Caste == caste_name
  coords_sub         <- gdf$coords[,, caste_filter]
  colony_sub         <- droplevels(gdf$Colony[caste_filter])
  coords_by_colony   <- coords.subset(coords_sub, group = colony_sub)
  colony_names       <- names(coords_by_colony)
  means_by_colony    <- lapply(coords_by_colony, mshape)
  num_cols           <- length(colony_names)

  # Displacement array for this caste
  disp_arr <- array(
    NA,
    dim      = c(num_lms, num_cols, num_cols),
    dimnames = list(Landmark   = paste0("LM", 1:num_lms),
                    FromColony = colony_names,
                    ToColony   = colony_names)
  )

  for (i in seq_len(num_cols)) {
    for (j in seq_len(num_cols)) {
      if (i != j) {
        disp_arr[, i, j] <- compute_disp(means_by_colony[[i]],
                                         means_by_colony[[j]])
      }
    }
  }

  # Export displacement table
  df_long <- melt(disp_arr, varnames = c("Landmark", "FromColony", "ToColony"))
  names(df_long)[4] <- "Displacement"
  df_long <- df_long[!is.na(df_long$Displacement), ]
  df_wide <- df_long %>%
    unite(Comparison, FromColony, ToColony, sep = "_vs_") %>%
    pivot_wider(names_from = Landmark, values_from = Displacement)

  colony_disp_wide_list[[caste_name]] <- df_wide
  write.csv(df_wide,
            paste0("output/", caste_name, "_colony_disp.csv"),
            row.names = FALSE)

  # Generate deformation grids (PDF)
  pdf(paste0("output/SuppFig2_colony_deformations_", caste_name, ".pdf"),
      width = 8, height = 8)

  # Plot all unique pairwise colony comparisons
  pairs_done <- character(0)
  for (i in seq_len(num_cols)) {
    for (j in seq_len(num_cols)) {
      if (i >= j) next
      pair_key <- paste(sort(c(colony_names[i], colony_names[j])), collapse = "_vs_")
      if (pair_key %in% pairs_done) next
      pairs_done <- c(pairs_done, pair_key)

      plotRefToTarget(
        means_by_colony[[i]],
        means_by_colony[[j]],
        method = "TPS",
        links  = links,
        mag    = 2,
        main   = paste0(caste_name, ": colony ", colony_names[i],
                        " → ", colony_names[j])
      )
    }
  }
  dev.off()
  cat(sprintf("Supp Fig 2 (%s): colony deformation grids saved.\n", caste_name))
}

# =============================================================================
# SUPPLEMENTARY FIGURE 3: Patriline-level deformation grids (gynes and workers)
# =============================================================================
# Males excluded: haploid, no patriline assignment possible.

patriline_disp_wide_list <- list()

for (caste_name in c("gyne", "worker")) {

  caste_filter       <- gdf$Caste == caste_name
  coords_sub         <- gdf$coords[,, caste_filter]
  colony_sub         <- droplevels(gdf$Colony[caste_filter])
  patriline_sub      <- droplevels(gdf$Patriline[caste_filter])
  coords_by_colony   <- coords.subset(coords_sub, group = colony_sub)
  colony_names       <- names(coords_by_colony)

  caste_colony_disp_list <- list()

  pdf(paste0("output/SuppFig3_patriline_deformations_", caste_name, ".pdf"),
      width = 8, height = 8)

  for (col_name in colony_names) {

    col_filter         <- (colony_sub == col_name)
    coords_in_colony   <- coords_by_colony[[col_name]]
    patriline_in_col   <- droplevels(patriline_sub[col_filter])
    coords_by_pat      <- coords.subset(coords_in_colony, group = patriline_in_col)
    patriline_names    <- names(coords_by_pat)

    # Need at least 2 patrilines for pairwise comparison
    if (length(patriline_names) < 2) next

    means_by_pat <- lapply(coords_by_pat, mshape)
    num_pats     <- length(patriline_names)

    disp_arr_pat <- array(
      NA,
      dim      = c(num_lms, num_pats, num_pats),
      dimnames = list(Landmark      = paste0("LM", 1:num_lms),
                      FromPatriline = patriline_names,
                      ToPatriline   = patriline_names)
    )

    for (i in seq_len(num_pats)) {
      for (j in seq_len(num_pats)) {
        if (i != j) {
          disp_arr_pat[, i, j] <- compute_disp(means_by_pat[[i]],
                                               means_by_pat[[j]])
        }
      }
    }

    # Export this colony's patriline displacement table
    df_long_pat <- melt(disp_arr_pat,
                        varnames = c("Landmark", "FromPatriline", "ToPatriline"))
    names(df_long_pat)[4] <- "Displacement"
    df_long_pat <- df_long_pat[!is.na(df_long_pat$Displacement), ]
    df_wide_pat <- df_long_pat %>%
      unite(Comparison, FromPatriline, ToPatriline, sep = "_vs_") %>%
      pivot_wider(names_from = Landmark, values_from = Displacement)

    caste_colony_disp_list[[col_name]] <- df_wide_pat

    # Deformation grids for this colony's patriline pairs
    for (i in seq_len(num_pats)) {
      for (j in seq_len(num_pats)) {
        if (i >= j) next
        plotRefToTarget(
          means_by_pat[[i]],
          means_by_pat[[j]],
          method = "TPS",
          links  = links,
          mag    = 2,
          main   = paste0(caste_name, " colony ", col_name, ": patriline ",
                          patriline_names[i], " → ", patriline_names[j])
        )
      }
    }
  }

  dev.off()

  # Combine all colonies into one data frame per caste
  caste_all_df <- bind_rows(
    lapply(names(caste_colony_disp_list), function(col) {
      df <- caste_colony_disp_list[[col]]
      if (is.null(df)) return(NULL)
      df$Colony <- col
      df
    })
  ) %>% select(Colony, everything())

  patriline_disp_wide_list[[caste_name]] <- caste_all_df
  write.csv(caste_all_df,
            paste0("output/all_", caste_name, "_patriline_colonies_disp.csv"),
            row.names = FALSE)

  cat(sprintf("Supp Fig 3 (%s): patriline deformation grids and displacement table saved.\n",
              caste_name))
}

cat("\nSupplementary figure generation complete.\n")
