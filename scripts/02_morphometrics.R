# =============================================================================
# 02_morphometric_analyses.R
# Vespula maculifrons Wing Shape Study
# Murray, Orr, Goodisman, Gosse, Hill, Kovacs
#
# PURPOSE: Perform all geometric morphometric analyses:
#   - Procrustes ANOVA to test for caste, colony, and patriline effects
#     on wing shape (Table 1 in manuscript)
#   - PCA-based visualization of caste-level shape variation (Figure 2A)
#   - Deformation grid visualizations of mean shape differences among castes
#     (Figure 2 B-D)
#   - Ten-panel PCA plots by colony and patriline (Figure 4)
#
# INPUT:   output/01_gpa_workspace.RData (Y.gpa, gdf, PCA, variance_props)
#
# OUTPUT:  Table 1 results (printed to console and saved as CSV)
#          Figure 2 (caste PCA + deformation grids)
#          Figure 4 (10-panel patriline PCA)
#
# DEPENDENCIES: 00_setup.R, 01_data_import_and_GPA.R
# =============================================================================

source("scripts/00_setup.R")
load("output/01_gpa_workspace.RData")

# ---- 1. Procrustes ANOVA: Table 1A - All castes ----------------------------
# Tests whether caste and colony explain significant variation in wing shape.
# Colony is nested within caste to account for the hierarchical structure.
# 1000 permutations used for significance testing.

fit_all <- procD.lm(
  coords ~ Caste + Colony/Caste,
  data  = gdf,
  iter  = 999,
  print.progress = FALSE
)
summary(fit_all)
# Table 1A: 

# Save Table 1A results
table1A <- as.data.frame(summary(fit_all)$table)
write.csv(table1A, "output/Table1A_allcastes_procD.csv")

# ---- 2. Procrustes ANOVA: Table 1B - Gynes only ----------------------------
# Tests whether colony and patriline explain wing shape variation within gynes.
# Patriline is nested within colony.

gyne_idx <- which(gdf$Caste == "gyne")


gdf_gyne <- geomorph.data.frame(
  coords    = gdf$coords[,,gyne_idx],
  Csize     = gdf$Csize[gyne_idx],
  Colony    = droplevels(gdf$Colony[gyne_idx]),
  Patriline = droplevels(gdf$Patriline[gyne_idx])
)

fit_gyne <- procD.lm(
  coords ~ Colony + Patriline,
  data  = gdf_gyne,
  iter  = 999,
  print.progress = FALSE
)

summary(fit_gyne)
# Table 1B: 

table1B <- as.data.frame(summary(fit_gyne)$table)
write.csv(table1B, "output/Table1B_gynes_procD.csv")

# ---- 3. Procrustes ANOVA: Table 1C - Workers only --------------------------

worker_idx <- which(gdf$Caste == "worker")

gdf_worker <- geomorph.data.frame(
  coords    = gdf$coords[,,worker_idx],
  Csize     = gdf$Csize[worker_idx],
  Colony    = droplevels(gdf$Colony[worker_idx]),
  Patriline = droplevels(gdf$Patriline[worker_idx])
)

fit_worker <- procD.lm(
  coords ~ Colony + Patriline,
  data  = gdf_worker,
  iter  = 999,
  print.progress = FALSE
)
summary(fit_worker)
# Table 1C: Colony 

table1C <- as.data.frame(summary(fit_worker)$table)
write.csv(table1C, "output/Table1C_workers_procD.csv")

# ---- 4. Figure 2A: PCA of wing shape by caste ------------------------------
# Plot PC1 vs PC2 for all 965 individuals, colored by caste.
# Symbol size and transparency help with overplotting.

pc_df <- data.frame(
  PC1   = PCA$x[, 1],
  PC2   = PCA$x[, 2],
  Caste = Y.gpa$caste
)

fig2A <- ggplot(pc_df, aes(x = PC1, y = PC2, colour = Caste)) +
  geom_point(alpha = 0.65, size = 1.4) +
  scale_colour_manual(values = caste_cols) +
  labs(
    x = sprintf("PC1 (%.2f%%)", variance_props[1]),
    y = sprintf("PC2 (%.2f%%)", variance_props[2]),
    colour = "Caste"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right")

print(fig2A)

# ---- 5. Deformation grids: Figure 2 B-D ------------------------------------
# Visualize mean shape differences between pairs of castes using TPS grids
# with heatmap coloring by landmark displacement magnitude.

# Compute mean shapes for each caste
mean_gyne   <- mshape(Y.gpa$coords[,, Y.gpa$caste == "gyne"])
mean_worker <- mshape(Y.gpa$coords[,, Y.gpa$caste == "worker"])
mean_male   <- mshape(Y.gpa$coords[,, Y.gpa$caste == "male"])

# Define landmark links for plotting (connect adjacent vein intersection points)
ref   <- mshape(Y.gpa$coords)
links <- define.links(ref)

# Helper function: TPS deformation grid with heatmap overlay
# ref_shape  = reference (starting) mean configuration
# tar_shape  = target mean configuration
# title_text = panel label string
plot_deformation <- function(ref_shape, tar_shape, title_text, mag = 2) {
  disp    <- sqrt(rowSums((tar_shape - ref_shape)^2))
  col_vec <- heat.colors(length(disp))[rank(disp)]
  
  plot(ref_shape[, 1], ref_shape[, 2],
       type = "n", asp = 1,
       xlab = "X", ylab = "Y",
       main = title_text)
  
  plotRefToTarget(ref_shape, tar_shape,
                  method = "TPS",
                  links  = links,
                  mag    = mag)
  
  points(ref_shape[, 1], ref_shape[, 2],
         pch = 21, bg = col_vec, cex = 2, lwd = 0.5)
}

# Figure 2B: Gyne → Male
par(mfrow = c(1, 3))
plot_deformation(mean_gyne,   mean_male,   "Gyne → Male")
# Figure 2C: Worker → Male
plot_deformation(mean_worker, mean_male,   "Worker → Male")
# Figure 2D: Gyne → Worker
plot_deformation(mean_gyne,   mean_worker, "Gyne → Worker")
par(mfrow = c(1, 1))

# ---- 6. Figure 4: Ten-panel PCA by colony and patriline --------------------
# Each panel shows PC1 vs PC2 for one colony.
# Symbol shape = caste (triangle: gyne, circle: worker, square: male)
# Color = patriline identity within that colony

# Build plotting data frame
fig4_df <- data.frame(
  PC1       = PCA$x[, 1],
  PC2       = PCA$x[, 2],
  Caste     = Y.gpa$caste,
  Colony    = Y.gpa$colony,
  Patriline = Y.gpa$patriline
)

# Caste symbol shapes
caste_shapes <- c("gyne" = 24, "worker" = 21, "male" = 22)

# Patriline colors: show up to 6 patrilines per colony
# Males are shown in a neutral grey since they have no patriline assignment
patriline_cols <- c(
  "1" = "#1b9e77", "2" = "#d95f02", "3" = "#7570b3",
  "4" = "#e7298a", "5" = "#66a61e", "6" = "#e6ab02",
  "male" = "grey60"
)

# Extract within-colony patriline number (the part after the underscore)
fig4_df <- fig4_df %>%
  mutate(
    pat_num = case_when(
      grepl("_", as.character(Patriline)) ~
        sub("^[0-9]+_", "", as.character(Patriline)),
      TRUE ~ "male"  # males coded as colony number only, no underscore
    ),
    pat_num = factor(pat_num,
                     levels = c("1","2","3","4","5","6","male"))
  )

# Generate 10-panel plot (2 rows x 5 columns)
colony_levels <- levels(Y.gpa$colony)

fig4 <- ggplot(fig4_df,
               aes(x = PC1, y = PC2,
                   shape = Caste,
                   fill  = pat_num)) +
  geom_point(size = 1.8, alpha = 0.8, colour = "black", stroke = 0.3) +
  scale_shape_manual(values = caste_shapes) +
  scale_fill_manual(values  = patriline_cols,
                    na.value = "grey80",
                    name = "Patriline") +
  facet_wrap(~ Colony, nrow = 2, ncol = 5) +
  labs(
    x = sprintf("PC1 (%.1f%%)", variance_props[1]),
    y = sprintf("PC2 (%.1f%%)", variance_props[2]),
    shape = "Caste"
  ) +
  theme_bw(base_size = 9) +
  theme(
    legend.position = "right",
    strip.background = element_rect(fill = "grey90"),
    panel.grid.minor = element_blank()
  )

print(fig4)

ggsave("output/Figure4_patriline_PCA.pdf", fig4,
       width = 12, height = 6, units = "in")

cat("\nMorphometric analyses complete.\n")
cat("Table 1 results saved to output/Table1A_allcastes_procD.csv, etc.\n")
cat("Figure 4 saved to output/Figure4_patriline_PCA.pdf\n")