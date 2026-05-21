# =============================================================================
# 00_setup.R
# Vespula maculifrons Wing Shape Study
# Murray, Orr, Goodisman, Gosse, Hill, Kovacs
#
# PURPOSE: Install and load all packages required for this analysis.
#          Run this script once before running any other script.
#
# PACKAGES:
#   geomorph    - geometric morphometrics (GPA, PCA, Procrustes ANOVA)
#   StereoMorph - landmark digitizing and shape file reading
#   MCMCglmm    - Bayesian mixed-effects models for heritability and G-matrix
#   Matrix      - sparse matrix operations (relatedness matrix inversion)
#   lme4        - linear mixed-effects models (distance analyses)
#   ggplot2     - plotting
#   patchwork   - multi-panel figure assembly
#   dplyr       - data manipulation
#   tidyr       - data reshaping
#   coda        - MCMC diagnostics and HPD intervals
#   reshape2    - array reshaping for supplementary displacement tables
#   broom.mixed - tidy model output from lme4
#   scales      - axis formatting in ggplot2
# =============================================================================

# ---- Install packages (run once, then comment out) --------------------------

# install.packages("geomorph")
# install.packages("StereoMorph")
# install.packages("MCMCglmm")
# install.packages("Matrix")
# install.packages("lme4")
# install.packages("ggplot2")
# install.packages("patchwork")
# install.packages("dplyr")
# install.packages("tidyr")
# install.packages("coda")
# install.packages("reshape2")
# install.packages("broom.mixed")
# install.packages("scales")

# For geomorph from GitHub (stable branch):
# install.packages("devtools")
# devtools::install_github("geomorphR/geomorph", ref = "Stable", build_vignettes = TRUE)

# ---- Load packages ----------------------------------------------------------

library(geomorph)    # v4.0.9
library(StereoMorph) # v1.6.7
library(MCMCglmm)    # v2.36
library(Matrix)
library(lme4)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(coda)
library(reshape2)
library(broom.mixed)
library(scales)

# ---- Consistent color palette used across all figures ----------------------
# These colors are used in Figures 2, 3, 4, and 5

caste_cols <- c(
  "gyne"   = "#2F6F5E",  # Forest Green
  "worker" = "#C46A2A",  # Burnt Sienna
  "male"   = "#A5C9F0"   # Light Sky Blue
)

# Trait colors for worker vs gyne PCs (used in Figure 5)
trait_cols <- c(
  "worker_PC1" = "#C46A2A",
  "worker_PC2" = "#C46A2A",
  "gyne_PC1"   = "#2F6F5E",
  "gyne_PC2"   = "#2F6F5E"
)
