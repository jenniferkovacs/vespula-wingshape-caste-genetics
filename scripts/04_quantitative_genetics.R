# =============================================================================
# 04_quantitative_genetics.R
# Vespula maculifrons Wing Shape Study
# Murray, Orr, Goodisman, Gosse, Hill, Kovacs
#
# PURPOSE: Quantify the heritable and cross-caste genetic structure of wing
#          shape using two complementary approaches:
#
#   PART A – DISTANCE-BASED ANALYSIS (Figure 5A–C)
#     Worker–gyne Euclidean distances in PC1–PC4 space are modeled with
#     linear mixed-effects models (lme4). Full-sister worker–gyne pairs
#     (same patriline) are compared to half-sister pairs (different patrilines
#     within the same colony). A significant fixed effect of relationship
#     indicates a cross-caste "family signature" in wing shape.
#
#   PART B – G-MATRIX ANALYSIS (Figure 5A–C)
#     A 4-trait Bayesian mixed-effects model (MCMCglmm) estimates the
#     patriline-level genetic variance-covariance matrix (G-matrix) for:
#       - Worker PC1, Worker PC2, Gyne PC1, Gyne PC2
#     Cross-caste genetic correlations (rG) between worker and gyne PCs are
#     derived from the posterior distribution of the G-matrix.
#
#   PART C – HERITABILITY (Figure 5D)
#     Narrow-sense heritability (h²) of PC1 is estimated separately for
#     workers and gynes using the posterior variance components from the
#     4-trait model.
#
# NOTE ON POWER:
#   With approximately 2–4 patrilines per colony across 10 colonies, the
#   design has limited power to precisely estimate genetic correlations.
#   The near-zero rG estimates with wide HPD intervals are consistent with
#   weak genetic constraint but could also reflect insufficient power to
#   detect correlations of small-to-moderate magnitude. See manuscript
#   Discussion and Limitations section for full treatment.
#
# NOTE ON MALE EXCLUSION:
#   Males in Hymenoptera are haploid and develop from unfertilized eggs.
#   Patriline assignment and heritability estimation are therefore not
#   applicable to males. All analyses in this script use workers and gynes only.
#
# RUNTIME WARNING:
#   The MCMCglmm model (mod_4_final) takes approximately 30–60 minutes
#   to run with nitt = 800000. A saved RData file is provided so you can
#   load the fitted model directly without re-running:
#     load("output/04_mod4final.RData")
#
# INPUT:   output/01_gpa_workspace.RData (Y.gpa, gdf, PCA, variance_props)
#
# OUTPUT:  mod_4_final (MCMCglmm object, saved to output/04_mod4final.RData)
#          Figure 5 (PDF and PNG)
#          CSV tables of fixed effects, variance components, rG estimates
#
# DEPENDENCIES: 00_setup.R, 01_data_import_and_GPA.R
# =============================================================================

source("scripts/00_setup.R")
load("output/01_gpa_workspace.RData")

# =============================================================================
# PART A: Build individual-level data frame with PC scores
# =============================================================================

# Use PC1–PC4 (covering > 70% of total shape variation)
k <- 4
PC_scores <- as.data.frame(PCA$x[, 1:k])
colnames(PC_scores) <- paste0("PC", 1:k)

# Individual-level metadata + PC scores
ind_df <- data.frame(
  ID        = as.factor(Y.gpa$ID),
  caste     = droplevels(as.factor(Y.gpa$caste)),
  patriline = droplevels(as.factor(Y.gpa$patriline)),
  colony    = droplevels(as.factor(Y.gpa$colony)),
  PC_scores
)

# Restrict to workers and gynes (males excluded; see NOTE above)
ind_wg <- ind_df %>%
  filter(caste %in% c("worker", "gyne")) %>%
  droplevels()

cat(sprintf("Worker-gyne dataset: %d individuals (%d workers, %d gynes)\n",
            nrow(ind_wg),
            sum(ind_wg$caste == "worker"),
            sum(ind_wg$caste == "gyne")))

# =============================================================================
# PART B: Worker–gyne Procrustes distance model
# =============================================================================

# ---- B1. Compute all worker–gyne Procrustes distances within colonies -------
# For each colony, form all worker × gyne pairs and compute the full
# Procrustes distance between their aligned landmark configurations.
# Label each pair as "full" (same patriline) or "half" (different patrilines).

compute_pairs_in_colony <- function(df_colony) {
  w_ids <- which(Y.gpa$ID %in% df_colony$ID[df_colony$caste == "worker"])
  g_ids <- which(Y.gpa$ID %in% df_colony$ID[df_colony$caste == "gyne"])

  if (length(w_ids) == 0 | length(g_ids) == 0) return(NULL)

  grid <- expand.grid(w = w_ids, g = g_ids)

  # Full Procrustes distance between already-aligned shapes
  dists <- vapply(seq_len(nrow(grid)), function(i) {
    W <- Y.gpa$coords[,, grid$w[i]]
    G <- Y.gpa$coords[,, grid$g[i]]
    sqrt(sum((W - G)^2))
  }, numeric(1))

  w_meta <- ind_wg[match(Y.gpa$ID[grid$w], ind_wg$ID), ]
  g_meta <- ind_wg[match(Y.gpa$ID[grid$g], ind_wg$ID), ]

  rel <- ifelse(w_meta$patriline == g_meta$patriline, "full", "half")

  data.frame(
    colony       = w_meta$colony,
    pat_w        = w_meta$patriline,
    pat_g        = g_meta$patriline,
    ID_w         = w_meta$ID,
    ID_g         = g_meta$ID,
    distance     = dists,
    relationship = factor(rel, levels = c("full", "half"))
  )
}

by_colony <- split(ind_wg, ind_wg$colony)
pair_list <- lapply(by_colony, compute_pairs_in_colony)
pair_list <- pair_list[!vapply(pair_list, is.null, logical(1))]
pair_df   <- dplyr::bind_rows(pair_list)

# Create patriline-pair factor for random effects model
pair_df <- pair_df %>%
  mutate(
    pat_pair = interaction(
      pmin(as.character(pat_w), as.character(pat_g)),
      pmax(as.character(pat_w), as.character(pat_g)),
      drop = TRUE
    )
  )

cat(sprintf("\nTotal worker-gyne pairs: %d\n", nrow(pair_df)))
cat("Relationship breakdown:\n")
print(table(pair_df$relationship))

# ---- B2. Mixed models on distances -----------------------------------------

# Model 1: colony as random effect only
mod_dist <- lmer(
  distance ~ relationship + (1 | colony),
  data = pair_df
)
summary(mod_dist)
# "relationshiphalf" estimate = extra distance for half-sister pairs

# Model 2: colony + patriline-pair random effects
mod_dist2 <- lmer(
  distance ~ relationship + (1 | colony) + (1 | pat_pair),
  data = pair_df
)
summary(mod_dist2)

# Variance component proportions
vc <- VarCorr(mod_dist2)
vc_df <- as.data.frame(vc)
vc_df$PropVar <- vc_df$vcov / sum(vc_df$vcov)
print(vc_df)

# Then extract variance components by group name
rand_df <- tibble::tibble(
  component = c("Patriline pair", "Colony", "Residual"),
  var = c(
    vc_df$vcov[vc_df$grp == "pat_pair"],
    vc_df$vcov[vc_df$grp == "colony"],
    vc_df$vcov[vc_df$grp == "Residual"]
  )
) %>% mutate(prop = var / sum(var))

# Save fixed effects table
fix_tab <- broom.mixed::tidy(mod_dist2, effects = "fixed") %>%
  mutate(term = recode(term,
                       "(Intercept)"      = "Intercept (full sisters)",
                       "relationshiphalf" = "relationship = half sisters"))
write.csv(fix_tab,  "output/Table_FixedEffects_DistanceModel.csv", row.names = FALSE)
write.csv(vc_df,    "output/Table_VarComp_DistanceModel.csv",      row.names = FALSE)

# ---- B3. Build distances_df for plotting -----------------------------------
distances_df <- pair_df %>%
  mutate(
    relationship = recode(relationship,
                          full = "Full sisters",
                          half = "Half sisters"),
    relationship = factor(relationship,
                          levels = c("Full sisters", "Half sisters"))
  )

# Model-predicted means for overlay on violin plot
fix_df <- tibble::tibble(
  relationship = c("Full sisters", "Half sisters"),
  estimate = c(
    fixef(mod_dist2)["(Intercept)"],
    fixef(mod_dist2)["(Intercept)"] + fixef(mod_dist2)["relationshiphalf"]
  ),
  se = c(
    summary(mod_dist2)$coefficients["(Intercept)", "Std. Error"],
    summary(mod_dist2)$coefficients["relationshiphalf", "Std. Error"]
  )
)
# Then extract variance components by group name
rand_df <- tibble::tibble(
  component = c("Patriline pair", "Colony", "Residual"),
  var = c(
    vc_df$vcov[vc_df$grp == "pat_pair"],
    vc_df$vcov[vc_df$grp == "colony"],
    vc_df$vcov[vc_df$grp == "Residual"]
  )
) %>% mutate(prop = var / sum(var))

# =============================================================================
# PART C: 4-trait MCMCglmm G-matrix model
# =============================================================================

# ---- C1. Build 4-trait stacked dataset (worker_PC1, worker_PC2, gyne_PC1, gyne_PC2)

w  <- ind_wg %>% filter(caste == "worker")
g  <- ind_wg %>% filter(caste == "gyne")

trait_levels <- c("worker_PC1", "worker_PC2", "gyne_PC1", "gyne_PC2")

w1 <- w %>% mutate(Trait = factor("worker_PC1", levels = trait_levels), y = PC1)
w2 <- w %>% mutate(Trait = factor("worker_PC2", levels = trait_levels), y = PC2)
g1 <- g %>% mutate(Trait = factor("gyne_PC1",   levels = trait_levels), y = PC1)
g2 <- g %>% mutate(Trait = factor("gyne_PC2",   levels = trait_levels), y = PC2)

mcmc_data_4 <- bind_rows(w1, w2, g1, g2) %>%
  mutate(
    Trait     = droplevels(Trait),
    patriline = droplevels(patriline),
    colony    = droplevels(colony)
  )

cat(sprintf("\n4-trait stacked dataset: %d rows\n", nrow(mcmc_data_4)))
print(table(mcmc_data_4$Trait))

# ---- C2. Patriline relatedness matrix K ------------------------------------
# Workers and gynes within a colony share a single mother (the founding queen)
# but have different fathers (polyandry). Half-sibling relatedness among
# individuals from different patrilines within a colony = 0.25
# (1/2 of the maternal contribution for diploid females).
# Self-relatedness on the diagonal = 1.

pat_levels <- levels(mcmc_data_4$patriline)
n_pat      <- length(pat_levels)

Kmat <- matrix(0, nrow = n_pat, ncol = n_pat,
               dimnames = list(pat_levels, pat_levels))
diag(Kmat) <- 1

# Map patriline to colony (first occurrence)
pat_to_colony <- tapply(mcmc_data_4$colony,
                        mcmc_data_4$patriline,
                        function(x) x[1])

# Set off-diagonals for patrilines within the same colony
for (i in seq_len(n_pat)) {
  for (j in seq_len(n_pat)) {
    if (i == j) next
    if (pat_to_colony[pat_levels[i]] == pat_to_colony[pat_levels[j]]) {
      Kmat[i, j] <- 0.25  # maternal half-sib relatedness
    }
  }
}

# Ensure positive definiteness and invert for MCMCglmm ginverse argument
Knear <- as.matrix(nearPD(Kmat)$mat)
Kinv  <- solve(Knear)
Kinv  <- as(Kinv, "dgCMatrix")
names(dimnames(Kinv)) <- c("patriline", "patriline")

# ---- C3. Prior specification for 4-trait G-matrix model --------------------
# G1: patriline-level 4x4 variance-covariance matrix (the G-matrix of interest)
# G2: colony-level 4x4 variance-covariance matrix (accounts for shared environment)
# R:  residual 4x4 variance-covariance matrix
# Weakly informative priors: diagonal V scaled to expected small variances,
# nu = 5 degrees of freedom (relatively uninformative for this dataset size).

prior_4 <- list(
  G = list(
    G1 = list(V = diag(4) * 0.02, nu = 5),  # patriline
    G2 = list(V = diag(4) * 0.02, nu = 5)   # colony
  ),
  R = list(V = diag(4) * 0.98, nu = 5)       # residual
)

# ---- C4. Fit MCMCglmm model ------------------------------------------------
# WARNING: This takes ~30–60 minutes. Set RUN_MODEL = FALSE to load saved
# results instead.

RUN_MODEL <- TRUE   # set to FALSE to load saved model

if (RUN_MODEL) {

  cat("\nFitting 4-trait MCMCglmm model (this may take 30-60 minutes)...\n")

  set.seed(42)  # for reproducibility of MCMC sampling

  mod_4_final <- MCMCglmm(
    fixed    = y ~ Trait - 1,
    random   = ~ us(Trait):patriline + us(Trait):colony,
    rcov     = ~ us(Trait):units,
    family   = "gaussian",
    data     = mcmc_data_4,
    ginverse = list(patriline = Kinv),
    prior    = prior_4,
    nitt     = 800000,  # total iterations
    burnin   = 200000,  # burn-in iterations discarded
    thin     = 300,     # thinning interval (≈ 2000 posterior samples)
    verbose  = TRUE
  )

  save(mod_4_final, file = "output/04_mod4final.RData")
  cat("Model saved to output/04_mod4final.RData\n")

} else {
  load("output/04_mod4final.RData")
  cat("Loaded saved model from output/04_mod4final.RData\n")
}

# ---- C5. MCMC diagnostics --------------------------------------------------
# Inspect mixing (trace plots) and effective sample sizes.
# Effective size should be >> 200 for reliable posterior estimates.

cat("\nEffective sample sizes (should be > 200):\n")
print(effectiveSize(mod_4_final$VCV))

cat("\nAutocorrelation at lag 1 (should be < 0.1 after thinning):\n")
print(autocorr.diag(mod_4_final$VCV)[2, ])  # lag-1 row

# Visual trace plots (uncomment to view)
# plot(mod_4_final$VCV[, c("Traitworker_PC1:Traitworker_PC1.patriline",
#                           "Traitgyne_PC1:Traitgyne_PC1.patriline")])

# ---- C6. Extract G-matrix: variances and genetic correlations ---------------

vcv_names <- colnames(mod_4_final$VCV)
pat_cols  <- grep("patriline$", vcv_names, value = TRUE)

traits  <- c("worker_PC1", "worker_PC2", "gyne_PC1", "gyne_PC2")
n_trait <- length(traits)

# Helper: rebuild 4x4 covariance matrix from one named VCV vector
rebuild_mat <- function(vec, traits) {
  m <- matrix(NA, nrow = length(traits), ncol = length(traits),
              dimnames = list(traits, traits))
  for (i in seq_along(traits)) {
    for (j in seq_along(traits)) {
      lab_ij <- paste0("Trait", traits[i], ":Trait", traits[j], ".patriline")
      lab_ji <- paste0("Trait", traits[j], ":Trait", traits[i], ".patriline")
      if (lab_ij %in% names(vec)) {
        m[i, j] <- vec[lab_ij]
      } else if (lab_ji %in% names(vec)) {
        m[i, j] <- vec[lab_ji]
      }
    }
  }
  m
}

# Helper: 95% credible interval
CI_fun <- function(x) as.numeric(quantile(x, probs = c(0.025, 0.975), na.rm = TRUE))

# Compute posterior covariance and correlation matrices for each MCMC sample
pat_CIs <- apply(mod_4_final$VCV[, pat_cols], 1, function(v) {
  cov_mat <- rebuild_mat(v, traits)
  cor_mat <- cov2cor(cov_mat)
  list(cov_mat = cov_mat, cor_mat = cor_mat)
})

n_iter    <- length(pat_CIs)
cov_array <- array(NA, dim = c(n_trait, n_trait, n_iter),
                   dimnames = list(traits, traits, NULL))
cor_array <- array(NA, dim = c(n_trait, n_trait, n_iter),
                   dimnames = list(traits, traits, NULL))

for (k in seq_len(n_iter)) {
  cov_array[,, k] <- pat_CIs[[k]]$cov_mat
  cor_array[,, k] <- pat_CIs[[k]]$cor_mat
}

# Genetic correlation table (upper triangle, with 95% HPD)
cor_list <- list()
idx      <- 1
for (i in seq_along(traits)) {
  for (j in seq_along(traits)) {
    if (j > i) {
      cor_samples <- cor_array[i, j, ]
      ci          <- CI_fun(cor_samples)
      cor_list[[idx]] <- data.frame(
        Trait1 = traits[i],
        Trait2 = traits[j],
        rG     = mean(cor_samples, na.rm = TRUE),
        rG_l95 = ci[1],
        rG_u95 = ci[2]
      )
      idx <- idx + 1
    }
  }
}
G_cor_tab <- do.call(rbind, cor_list)

cat("\n=== Patriline-level genetic correlations (rG) ===\n")
print(G_cor_tab %>% mutate(across(where(is.numeric), ~ round(.x, 4))))
write.csv(G_cor_tab, "output/Table_Gmatrix_Correlations.csv", row.names = FALSE)
# Key cross-caste correlations reported in manuscript:
cat(sprintf(
  "\nWorker PC1 vs Gyne PC1: rG = %.4f (95%% HPD: %.2f to %.2f)\n",
  G_cor_tab$rG[G_cor_tab$Trait1 == "worker_PC1" & G_cor_tab$Trait2 == "gyne_PC1"],
  G_cor_tab$rG_l95[G_cor_tab$Trait1 == "worker_PC1" & G_cor_tab$Trait2 == "gyne_PC1"],
  G_cor_tab$rG_u95[G_cor_tab$Trait1 == "worker_PC1" & G_cor_tab$Trait2 == "gyne_PC1"]
))
cat(sprintf(
  "Worker PC2 vs Gyne PC2: rG = %.4f (95%% HPD: %.2f to %.2f)\n",
  G_cor_tab$rG[G_cor_tab$Trait1 == "worker_PC2" & G_cor_tab$Trait2 == "gyne_PC2"],
  G_cor_tab$rG_l95[G_cor_tab$Trait1 == "worker_PC2" & G_cor_tab$Trait2 == "gyne_PC2"],
  G_cor_tab$rG_u95[G_cor_tab$Trait1 == "worker_PC2" & G_cor_tab$Trait2 == "gyne_PC2"]
))

# Variance table (diagonal of G-matrix)
var_list <- lapply(seq_along(traits), function(i) {
  var_samples <- cov_array[i, i, ]
  ci          <- CI_fun(var_samples)
  data.frame(Trait = traits[i],
             Var_G = mean(var_samples, na.rm = TRUE),
             Var_l95 = ci[1], Var_u95 = ci[2])
})
G_var_tab <- do.call(rbind, var_list)
write.csv(G_var_tab, "output/Table_Gmatrix_Variances.csv", row.names = FALSE)

# ---- C7. Heritability of PC1 within each female caste ----------------------
# h² = VA / (VA + VC + VR)
# VA = patriline variance, VC = colony variance, VR = residual variance

get_h2 <- function(trait_name) {
  vp <- mod_4_final$VCV[, paste0("Trait", trait_name, ":Trait", trait_name, ".patriline")]
  vc <- mod_4_final$VCV[, paste0("Trait", trait_name, ":Trait", trait_name, ".colony")]
  vr <- mod_4_final$VCV[, paste0("Trait", trait_name, ":Trait", trait_name, ".units")]
  h2_post <- vp / (vp + vc + vr)
  list(
    mean = mean(h2_post),
    HPD  = HPDinterval(as.mcmc(h2_post), prob = 0.95)
  )
}

h2_worker <- get_h2("worker_PC1")
h2_gyne   <- get_h2("gyne_PC1")

cat(sprintf("\nWorker PC1 heritability: h² = %.3f (95%% HPD: %.3f - %.3f)\n",
            h2_worker$mean, h2_worker$HPD[1], h2_worker$HPD[2]))
cat(sprintf("Gyne PC1 heritability:   h² = %.3f (95%% HPD: %.3f - %.3f)\n",
            h2_gyne$mean,   h2_gyne$HPD[1],   h2_gyne$HPD[2]))

herit_df <- tibble::tibble(
  trait = c("Worker PC1", "Gyne PC1"),
  h2    = c(h2_worker$mean,   h2_gyne$mean),
  lower = c(h2_worker$HPD[1], h2_gyne$HPD[1]),
  upper = c(h2_worker$HPD[2], h2_gyne$HPD[2])
)

# =============================================================================
# PART D: Figure 5 assembly
# =============================================================================

# -- Panel A: Variance components from distance model (bar chart)
pA1 <- ggplot(rand_df,
              aes(x = prop, y = factor(component, levels = rev(component)))) +
  geom_col(fill = "grey40") +
  scale_x_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = "Proportion of variance", y = NULL,
       title = "A) Random effects (distance model)") +
  theme_bw(base_size = 11) +
  theme(panel.grid.major.y = element_blank())

# -- Panel B: Fixed effects - full vs half sister distance
pA2 <- ggplot(fix_df, aes(x = relationship, y = estimate)) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = estimate - 1.96 * se,
                    ymax = estimate + 1.96 * se),
                width = 0) +
  labs(x = NULL, y = "Model-predicted distance",
       title = "B) Fixed effects") +
  theme_bw(base_size = 11)

# -- Panel C: Violin + boxplot of observed distances by relatedness
pC <- ggplot(distances_df,
             aes(x = relationship, y = distance, fill = relationship)) +
  geom_violin(trim = TRUE, alpha = 0.5, color = NA) +
  geom_boxplot(width = 0.15, outlier.size = 0.4,
               alpha = 0.9, color = "black") +
  geom_point(data = fix_df,
             aes(x = relationship, y = estimate),
             inherit.aes = FALSE, size = 2.5, color = "black") +
  scale_fill_manual(values = c("Full sisters"  = "#1b9e77",
                                "Half sisters" = "#d95f02")) +
  labs(x = NULL, y = "Worker–gyne Procrustes distance",
       title = "C) Cross-caste distances by relatedness") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")

# -- Panel D: Heritability estimates for worker and gyne PC1
pD <- ggplot(herit_df, aes(x = trait, y = h2)) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = NULL, y = expression(h^2),
       title = "D) Patriline-level heritabilities (PC1)") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

# Assemble Figure 5.1
fig5.1 <- (pA1 / pA2) | pC | pD
fig5.1 <- fig5 + plot_annotation(tag_levels = "A")
print(fig5.1)

ggsave("output/Figure5-1_quantitative_genetics.pdf", fig5.1,
       width = 12, height = 6, units = "in")
ggsave("output/Figure5-1_quantitative_genetics.png", fig5.1,
       width = 12, height = 6, units = "in", dpi = 300)

cat("\nQuantitative genetics analyses complete.\n")
cat("Figures saved to output/Figure5_quantitative_genetics.pdf/.png\n")



# ============================================================
# Figure 5: Quantitative genetic structure of wing shape
# Panel A: Full vs half sister Procrustes distances (violin)
# Panel B: Cross-caste genetic correlations from G-matrix
# Panel C: Patriline-level heritability for PC1
# ============================================================

library(ggplot2)
library(patchwork)
library(dplyr)

# ---- Panel A: Violin plot with model means overlaid ------------------------

pA <- ggplot(distances_df,
             aes(x = relationship, y = distance, fill = relationship)) +
  geom_violin(trim = TRUE, alpha = 0.5, color = NA) +
  geom_boxplot(width = 0.12, outlier.size = 0.3,
               alpha = 0.9, color = "black") +
  geom_point(data = fix_df,
             aes(x = relationship, y = estimate),
             inherit.aes = FALSE,
             size = 3, color = "black") +
  geom_errorbar(data = fix_df,
                aes(x = relationship,
                    ymin = estimate - 1.96 * se,
                    ymax = estimate + 1.96 * se),
                inherit.aes = FALSE,
                width = 0.06, linewidth = 0.8, color = "black") +
  annotate("segment",
           x = 1, xend = 2,
           y = max(distances_df$distance) * 0.95,
           yend = max(distances_df$distance) * 0.95,
           linewidth = 0.5) +
  annotate("segment",
           x = 1, xend = 1,
           y = max(distances_df$distance) * 0.95,
           yend = max(distances_df$distance) * 0.92,
           linewidth = 0.5) +
  annotate("segment",
           x = 2, xend = 2,
           y = max(distances_df$distance) * 0.95,
           yend = max(distances_df$distance) * 0.92,
           linewidth = 0.5) +
  annotate("text",
           x = 1.5, y = max(distances_df$distance) * 0.97,
           label = "***",
           size = 5, hjust = 0.5) +
  scale_fill_manual(values = c("Full sisters"  = "#1b9e77",
                               "Half sisters" = "#d95f02")) +
  labs(x = NULL,
       y = "Worker–gyne Procrustes distance",
       title = "A) Cross-caste shape similarity by relatedness") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")

# ---- Panel B: Genetic correlations from G-matrix ---------------------------
# Build a clean data frame of all cross-caste and within-caste rG estimates

G_cor_plot <- G_cor_tab %>%
  filter(
    !(Trait1 == "worker_PC1" & Trait2 == "gyne_PC2"),
    !(Trait1 == "worker_PC2" & Trait2 == "gyne_PC1")
  ) %>%
  mutate(
    pair_label = case_when(
      Trait1 == "worker_PC1" & Trait2 == "gyne_PC1"   ~ "wPC1 – gPC1",
      Trait1 == "worker_PC2" & Trait2 == "gyne_PC2"   ~ "wPC2 – gPC2",
      Trait1 == "worker_PC1" & Trait2 == "worker_PC2" ~ "wPC1 – wPC2",
      Trait1 == "gyne_PC1"   & Trait2 == "gyne_PC2"   ~ "gPC1 – gPC2"
    ),
    comparison_type = case_when(
      Trait1 == "worker_PC1" & Trait2 == "gyne_PC1" ~ "Cross-caste",
      Trait1 == "worker_PC2" & Trait2 == "gyne_PC2" ~ "Cross-caste",
      TRUE ~ "Within-caste"
    ),
    pair_label = factor(pair_label,
                        levels = c("wPC1 – gPC1", "wPC2 – gPC2",
                                   "wPC1 – wPC2", "gPC1 – gPC2"))
  ) %>%
  filter(!is.na(pair_label))

pB <- ggplot(G_cor_plot,
             aes(x = rG, y = pair_label, colour = comparison_type)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_errorbarh(aes(xmin = rG_l95, xmax = rG_u95),
                 height = 0.3, linewidth = 0.7) +
  geom_point(size = 3) +
  scale_colour_manual(values = c("Cross-caste"  = "#7570b3",
                                 "Within-caste" = "#999999"),
                      name = NULL) +
  scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
  labs(x = expression("Genetic correlation (" * r[G] * ")"),
       y = NULL,
       title = "B) Patriline-level genetic correlations") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")
# ---- Panel C: Heritability estimates ---------------------------------------

pC <- ggplot(herit_df,
             aes(x = trait, y = h2, colour = trait)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                width = 0.12, linewidth = 0.8) +
  geom_text(aes(y = upper + 0.05, label = "***"),
            colour = "black", size = 4) +
  scale_colour_manual(values = c("Worker PC1" = "#C46A2A",
                                 "Gyne PC1"   = "#2F6F5E")) +
  scale_y_continuous(limits = c(0, 1),
                     breaks = seq(0, 1, 0.25)) +
  labs(x = NULL,
       y = expression("Heritability (" * h^2 * ")"),
       title = "C) Patriline-level heritability (PC1)") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 20, hjust = 1))

# ---- Assemble Figure 5 -----------------------------------------------------

fig5 <- pA | pB | pC

ggsave("output/Figure5_quantitative_genetics.pdf", fig5,
       width = 12, height = 5, units = "in")
ggsave("output/Figure5_quantitative_genetics.png", fig5,
       width = 12, height = 5, units = "in", dpi = 300)


