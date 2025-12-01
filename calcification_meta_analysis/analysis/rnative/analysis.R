# ----------------------------
# Import packages and functions
# ----------------------------

load_packages <- function(packages) {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg, repos = "https://cloud.r-project.org")
    }
    library(pkg, character.only = TRUE)
  }
}

# load all required packages
load_packages(c("metafor", "rprojroot", "ggplot2", "dplyr", "ggstance"))

paper_conf_root <- find_root(
  criterion = has_dir("calcification_meta_analysis"),
  path = getwd()
)
functions_path <- file.path(paper_conf_root, "calcification_meta_analysis", "analysis", "rnative", "functions.R")
source(functions_path)
# make the 'figures' directory if it doesn't already exist
dir.create(file.path(paper_conf_root, "figures"), showWarnings = FALSE)

# ----------------------------
# Prepare data
# ----------------------------

dat <- read.csv(file.path(paper_conf_root, "data", "analysis_ready_calcification_data.csv"))

effect_type <- "st_relative_calcification"
# rename effect sizes and their variance for convenience
colnames(dat)[which(colnames(dat) == effect_type)] <- "yi"
colnames(dat)[which(colnames(dat) == paste0(effect_type, "_var"))] <- "vi"
# rename key moderators for convenience
dat$dt <- dat$delta_t
dat$dph <- dat$delta_ph

# remove extreme climatology corresponding to 10th and 90th percentiles at 2100 under SSP5-8.5
max_dt <- 4.71518004416422
min_dt <- -0.8509536893003847
max_ph <- 0.017361047261669333
min_ph <- -0.4439020009266399
clim_dat <- dat[(dat$dt > min_dt) & (dat$dt < max_dt) &
  (dat$dph > min_ph) & (dat$dph < max_ph), ]

# ----------------------------
# Define model
# ----------------------------

# nested random structure accounting for study and sample level variance
cg_random_structure <- ~ 1 | doi / ID
# linear additive model
mods_formula <- ~ dt + dph - 1

# ----------------------------
# Run a few models to compare outputs with the Python implementation
# ----------------------------
truncated_data <- dat[1:200, ]

# raw model (no moderators)
raw_model <- rma.mv(yi = truncated_data$yi, V = truncated_data$vi, random = ~ 1 | doi, data = truncated_data, method = "REML")
print("Raw model:")
print(raw_model)

# simple linear additive model with moderators
mods_model <- rma.mv(yi = truncated_data$yi, V = truncated_data$vi, random = ~ 1 | doi, data = truncated_data, mods = mods_formula, method = "REML")
print("Model with moderators:")
print(mods_model)

# modifying random effect structure
mods_random_model <- rma.mv(yi = truncated_data$yi, V = truncated_data$vi, random = cg_random_structure, data = truncated_data, mods = mods_formula, method = "REML")
print("Model with modified random effect structure:")
print(mods_random_model)

# -------------
# FIGURE S2 – Sensitivity analysis: climate filtering effect
# -------------

clim_out <- clim_filtering(dat, mods_formula, cg_random_structure, grouping_var = "core_grouping", dt_bounds = c(min_dt, max_dt), dph_bounds = c(min_ph, max_ph), plot = TRUE)
png(file.path(paper_conf_root, "figures", "FIG_S2.jpg"), width = 8, height = 8, units = "in", res = 300)
print(plot_coefficient_comparison(clim_out$results))
dev.off()

# now exclude Foraminifera since too affected by climate filtering
final_dat <- subset(clim_dat, tolower(core_grouping) != "foraminifera")

# -------------
# FIGURE S6 – Sensitivity analysis: excluding extreme effect sizes/variances
# -------------
extreme_filter_results <- extreme_filtering(final_dat, mods_formula, cg_random_structure, plot = FALSE, extreme_limits = c(-1000, 1000))
png(file.path(paper_conf_root, "figures", "FIG_S6.jpg"), width = 8, height = 8, units = "in", res = 300)
plot_extreme_filtering(extreme_filter_results$results)
dev.off()

# -------------
# FIGURE S7 – Cook's distance analysis: excluding extreme effect sizes/variances
# -------------
cooks_filter_results <- influence_filtering(final_dat, mods_formula, cg_random_structure, plot = TRUE)
png(file.path(paper_conf_root, "figures", "FIG_S7.jpg"), width = 8, height = 8, units = "in", res = 300)
print(plot_coefficient_comparison(cooks_filter_results$results))
dev.off()

# -------------
# FIGURE S8 – Exploring publication bias
# -------------

# Rhosenthal's fail-safe N test
fsn_result <- fsn(final_dat$yi, final_dat$vi)
print("Fail-safe N:")
print(fsn_result)

# Egger's regression test for funnel plot asymmetry
# raw model (no moderators)
egger_test_no_mods_result <- regtest(rma(yi = final_dat$yi, vi = final_dat$vi), model = "rma")
print("Egger's regression test with no moderators:")
print(egger_test_no_mods_result)

# Egger's regression test with moderators
egger_test_mods_result <- regtest(rma(yi, vi, mods = mods_formula, data = final_dat), model = "rma")
print("Egger's regression test with moderators:")
print(egger_test_mods_result)

# funnel plot
png(file.path(paper_conf_root, "figures", "FIG_S8.jpg"), width = 10, height = 10, units = "in", res = 300)
funnel(final_dat$yi, final_dat$vi,
  shade = c("white", "gray55", "gray75"),
  yaxs = "i", xaxs = "i",
  legend = TRUE, back = "gray90", hlines = NULL,
  xlab = "Percentage change in calcification rate",
  ylab = "Standard Error",
  level = c(.1, .05, .01),
  las = 1, digits = list(1L, 0),
  ylim = c(0, 88),
  xlim = c(-250, 250),
  mgp = c(3, 1, 0),
)
print(funnel_plot)
dev.off()
