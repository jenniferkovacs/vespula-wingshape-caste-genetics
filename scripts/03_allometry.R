# =============================================================================
# 03_allometry.R
# Vespula maculifrons Wing Shape Study
# Murray, Orr, Goodisman, Gosse, Hill, Kovacs
#
# PURPOSE: Investigate whether wing shape changes systematically with body size
#          (centroid size) within and among castes. Produces Figure 3.
#
# APPROACH:
#   1. Fit a pooled multivariate Procrustes allometry model across all castes
#      with a caste x log(Csize) interaction term to test whether allometric
#      trajectories differ among castes.
#   2. Extract 1D regression scores (RegScore) that project each individual's
#      wing shape onto the pooled allometric axis.
#   3. Regress RegScore on log(Csize) separately within each caste to estimate
#      caste-specific allometric slopes and their significance.
#   4. Produce Figure 3: scatter of RegScore vs log(Csize) with fitted lines
#      per caste.
#
# NOTE ON INTERPRETATION:
#   RegScore values are centered on zero by construction. Negative values
#   indicate shapes closer to the "small-size" end of the allometric trajectory;
#   positive values indicate shapes closer to the "large-size" end. The x-axis
#   (log centroid size) is also mean-centered within the pooled model, so
#   negative values simply reflect specimens smaller than the grand mean size.
#   This centering does not affect slope estimates or significance tests.
#
# NOTE ON GYNE PATTERN:
#   Gyne allometry shows a strong positive slope, but the result is partly
#   driven by one colony at the large-size extreme of the gyne distribution.
#   Results should be interpreted with appropriate caution (see manuscript
#   Discussion). A sensitivity analysis excluding that colony is recommended.
#
# INPUT:   output/01_gpa_workspace.RData (Y.gpa, gdf, PCA, variance_props)
#
# OUTPUT:  Figure 3 (PDF and PNG)
#          Allometry slope table (CSV)
#
# DEPENDENCIES: 00_setup.R, 01_data_import_and_GPA.R
# =============================================================================

source("scripts/00_setup.R")
load("output/01_gpa_workspace.RData")

# ---- 1. Test for caste differences in allometric trajectories ---------------
# Multivariate Procrustes linear model with caste x log(Csize) interaction.
# A significant interaction term confirms that allometric slopes differ
# among castes (i.e., trajectories are not parallel in shape space)

# Remove outlier gynes before fitting the allometry model
# 258_G1 (index 849) and 258_G10 (index 850): unusually small gynes from
# colony 258 that distort the pooled allometric axis
outlier_idx <- c(849, 850)
cat("Removing:", as.character(Y.gpa$ID[outlier_idx]), "\n")

gdf <- geomorph.data.frame(
  coords    = gdf$coords[,, -outlier_idx],
  Csize     = gdf$Csize[-outlier_idx],
  Caste     = droplevels(gdf$Caste[-outlier_idx]),
  Colony    = droplevels(gdf$Colony[-outlier_idx]),
  Patriline = droplevels(gdf$Patriline[-outlier_idx])
)

# Verify removal
cat(sprintf("gdf_clean: %d individuals (was 965)\n", dim(gdf_clean$coords)[3]))


groupAllometry <- procD.lm(
  coords ~ log(Csize) * Caste,
  data  = gdf,
  iter  = 999,
  RRPP  = TRUE,
  print.progress = FALSE
)
summary(groupAllometry)
# The significant Caste x log(Csize) interaction confirms divergent allometries.

# ---- 2. Extract 1D regression scores (RegScore) ----------------------------
# plotAllometry() projects each individual's shape onto the single
# multivariate allometric axis (common across all castes).
# The resulting RegScore provides a 1D summary of size-related shape change.

Pallom <- plotAllometry(
  groupAllometry,
  size   = gdf$Csize,
  logsz  = TRUE,
  method = "RegScore",
  pch    = 19,
  col    = as.numeric(gdf$Caste)
)

RegScore <- Pallom$RegScore  # one score per individual (n = 965)

# ---- 3. Build data frame for caste-specific regressions --------------------

df_allom <- data.frame(
  AllomScore = RegScore,
  lnCsize    = log(gdf$Csize),
  Caste      = gdf$Caste,
  Colony     = gdf$Colony,
  ID         = Y.gpa$ID        # add ID so we can filter by name
)


df_allom <- df_allom[!df_allom$ID %in% outlier_IDs, ]
cat(sprintf("df_allom after outlier removal: %d individuals\n", nrow(df_allom)))

# ---- 4. Fit caste-specific linear models -----------------------------------
# Regress RegScore on log(Csize) separately for each caste.
# This estimates the slope (and its SE, t, p) for each caste individually.

lm_gyne   <- lm(AllomScore ~ lnCsize,
                data = subset(df_allom, Caste == "gyne"))
lm_worker <- lm(AllomScore ~ lnCsize,
                data = subset(df_allom, Caste == "worker"))
lm_male   <- lm(AllomScore ~ lnCsize,
                data = subset(df_allom, Caste == "male"))

# Print slope summaries
cat("\n=== Gyne allometry (slope = size-shape relationship) ===\n")
print(summary(lm_gyne)$coefficients)
cat("\n95% CI for gyne slope:\n")
print(confint(lm_gyne)["lnCsize", ])

cat("\n=== Worker allometry ===\n")
print(summary(lm_worker)$coefficients)
cat("\n95% CI for worker slope:\n")
print(confint(lm_worker)["lnCsize", ])

cat("\n=== Male allometry ===\n")
print(summary(lm_male)$coefficients)
cat("\n95% CI for male slope:\n")
print(confint(lm_male)["lnCsize", ])

# Summary table of slopes
slope_tab <- data.frame(
  Caste = c("gyne", "worker", "male"),
  Slope = c(coef(lm_gyne)["lnCsize"],
            coef(lm_worker)["lnCsize"],
            coef(lm_male)["lnCsize"]),
  SE    = c(summary(lm_gyne)$coef["lnCsize","Std. Error"],
            summary(lm_worker)$coef["lnCsize","Std. Error"],
            summary(lm_male)$coef["lnCsize","Std. Error"]),
  t     = c(summary(lm_gyne)$coef["lnCsize","t value"],
            summary(lm_worker)$coef["lnCsize","t value"],
            summary(lm_male)$coef["lnCsize","t value"]),
  p     = c(summary(lm_gyne)$coef["lnCsize","Pr(>|t|)"],
            summary(lm_worker)$coef["lnCsize","Pr(>|t|)"],
            summary(lm_male)$coef["lnCsize","Pr(>|t|)"]),
  CI_lo = c(confint(lm_gyne)["lnCsize",1],
            confint(lm_worker)["lnCsize",1],
            confint(lm_male)["lnCsize",1]),
  CI_hi = c(confint(lm_gyne)["lnCsize",2],
            confint(lm_worker)["lnCsize",2],
            confint(lm_male)["lnCsize",2])
)

print(slope_tab)
write.csv(slope_tab, "output/allometry_slopes.csv", row.names = FALSE)

# ---- 5. Generate predicted regression lines per caste ----------------------

newdat <- expand.grid(
  lnCsize = seq(min(df_allom$lnCsize), max(df_allom$lnCsize), length.out = 200),
  Caste   = levels(df_allom$Caste)
)

# Predict from the pooled interaction model for line drawing
lm_allom <- lm(AllomScore ~ lnCsize * Caste, data = df_allom)
newdat$fit <- predict(lm_allom, newdata = newdat)

# ---- 6. Figure 3: Caste-specific allometric trajectories -------------------
# Points = individual RegScores colored by caste
# Lines = fitted caste-specific regression lines

pooled <- ggplot(df_allom, aes(x = lnCsize, y = AllomScore, colour = Caste)) +
  geom_point(alpha = 0.5, size = 1.2) +
  geom_line(data = newdat,
            aes(y = fit, colour = Caste),
            linewidth = 1.2) +
  scale_colour_manual(
    values = caste_cols,
    labels = c("gyne" = "Gyne", "worker" = "Worker", "male" = "Male")
  ) +
  labs(
    x      = "log(Centroid Size)",
    y      = "Allometric score (RegScore)",
    colour = "Caste",
    title  = "Caste-specific allometric trajectories for wing shape"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right")

print(fig3)

ggsave("output/Figure3_allometry.pdf", fig3,
       width = 7, height = 5, units = "in")
ggsave("output/Figure3_allometry.png", fig3,
       width = 7, height = 5, units = "in", dpi = 300)

# ---- 7. Sensitivity analysis: gyne allometry without outlier colony --------
# The manuscript notes that one colony drives the extreme end of the gyne
# allometric pattern. This block reruns the gyne model excluding that colony
# to assess robustness. Identify the outlier colony visually from Figure 3
# before setting outlier_colony below.

outlier_colony <- "235"

df_gyne_sens <- subset(df_allom, Caste == "gyne" & Colony != outlier_colony)

if (nrow(df_gyne_sens) > 0) {
  lm_gyne_sens <- lm(AllomScore ~ lnCsize, data = df_gyne_sens)
  cat("\n=== Gyne allometry: sensitivity analysis without colony", outlier_colony, "===\n")
  print(summary(lm_gyne_sens)$coefficients)
  cat("95% CI:", confint(lm_gyne_sens)["lnCsize", ], "\n")
  
  # Compare slope with and without the outlier colony
  cat("\n--- Comparison ---\n")
  cat(sprintf("Full gyne slope:    %.4f (SE = %.4f, p = %.4g)\n",
              coef(lm_gyne)["lnCsize"],
              summary(lm_gyne)$coef["lnCsize","Std. Error"],
              summary(lm_gyne)$coef["lnCsize","Pr(>|t|)"]))
  cat(sprintf("Without colony %s: %.4f (SE = %.4f, p = %.4g)\n",
              outlier_colony,
              coef(lm_gyne_sens)["lnCsize"],
              summary(lm_gyne_sens)$coef["lnCsize","Std. Error"],
              summary(lm_gyne_sens)$coef["lnCsize","Pr(>|t|)"]))
}
cat("\nAllometry analyses complete.\n")


# Within-colony gyne slopes
df_gyne_only <- subset(df_allom, Caste == "gyne")

within_colony_slopes <- df_gyne_only %>%
  group_by(Colony) %>%
  summarise(
    n     = n(),
    slope = coef(lm(AllomScore ~ lnCsize))[["lnCsize"]],
    p     = summary(lm(AllomScore ~ lnCsize))$coef["lnCsize","Pr(>|t|)"]
  ) %>%
  arrange(slope)

print(within_colony_slopes)



# Within-colony slopes for all three castes
get_within_slopes <- function(caste_name) {
  df_sub <- subset(df_allom, Caste == caste_name)
  df_sub %>%
    group_by(Colony) %>%
    summarise(
      n     = n(),
      slope = coef(lm(AllomScore ~ lnCsize))[["lnCsize"]],
      se    = summary(lm(AllomScore ~ lnCsize))$coef["lnCsize", "Std. Error"],
      p     = summary(lm(AllomScore ~ lnCsize))$coef["lnCsize", "Pr(>|t|)"],
      .groups = "drop"
    ) %>%
    mutate(
      Caste = caste_name,
      sig   = case_when(p < 0.05 ~ "*", p < 0.10 ~ "†", TRUE ~ "")
    )
}

slopes_all <- bind_rows(
  get_within_slopes("gyne"),
  get_within_slopes("worker"),
  get_within_slopes("male")
) %>%
  mutate(Caste = factor(Caste, levels = c("gyne", "worker", "male")))

# Plot: one panel per caste, colonies on y-axis, slope on x-axis
fig3 <- ggplot(slopes_all, aes(x = slope, y = reorder(Colony, slope), colour = Caste)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_errorbarh(aes(xmin = slope - 1.96 * se,
                     xmax = slope + 1.96 * se),
                 height = 0.3) +
  geom_point(size = 2.5) +
  geom_text(aes(label = sig), hjust = -0.5, size = 4, colour = "black") +
  scale_colour_manual(values = caste_cols) +
  facet_wrap(~ Caste, nrow = 1,
             labeller = labeller(Caste = c(gyne = "Gyne",
                                           worker = "Worker",
                                           male = "Male"))) +
  labs(x = "Within-colony allometric slope",
       y = "Colony") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none",
        strip.background = element_rect(fill = "grey90"))



print(fig3)

ggsave("output/Figure3_allometry.pdf", fig3,
       width = 7, height = 5, units = "in")
ggsave("output/Figure3_allometry.png", fig3,
       width = 7, height = 5, units = "in", dpi = 300)
