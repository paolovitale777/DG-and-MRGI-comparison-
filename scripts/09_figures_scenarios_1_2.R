library(dplyr)
library(tidyr)
library(ggplot2)
dir.create("results", showWarnings = FALSE)

# DG vs MRGI, fixed selection sizes: F6 (Elite Yield Trial) mean genetic value and genetic
# variance over time, plus the table of selection intensity (% selected) per step.
# Run from the repository root after the simulation scripts.

nReps        <- 5                 # must match nReps in simulation scripts
nYearsBurnIn <- 5                 # must match the burn-in scripts
yearsTotal   <- nYearsBurnIn + 25 # must match the burn-in scripts

# Folder with the .rds outputs and selection logs of the four scenarios (relative to the repository root)
dataDir <- "output"
# Folder with the selection logs used for the intensity table (relative to the repository root)
logDir  <- "output"

# Row label = required-gain scenario; method label = file prefix (scenarioName)
gainScenarios <- list(
  "Equal required gains (1, 1, 1, 1)" =
    c("DGI" = "Multitrait_DG_index_same_weights",
      "MRGI"          = "Multitrait_MRGI_index_same_weights"),
  "Uneven required gains (1, 0.1, 0.8, 0.1)" =
    c("DGI" = "Multitrait_DG_index_different_weights",
      "MRGI"          = "Multitrait_MRGI_index_uneven_weights"))
methodLevels <- c("DGI", "MRGI")


## Read in output data (scenarios, methods, reps, traits) ----
rawData <- list()
for (SCEN in names(gainScenarios)) {
  for (METHOD in names(gainScenarios[[SCEN]])) {
    for (REP in 1:nReps) {
      tmp  <- readRDS(file.path(dataDir, paste0(gainScenarios[[SCEN]][[METHOD]], "_", REP, ".rds")))
      keep <- grep("^(mean|var)_F6_EYT_t[1-5]$", names(tmp), value = TRUE)
      tmp  <- tmp[tmp$year >= nYearsBurnIn, c("year", keep)]
      tmp$scenario <- SCEN
      tmp$method   <- METHOD
      tmp$rep      <- REP
      rawData[[length(rawData) + 1]] <- tmp
    }
  }
}
rawData <- bind_rows(rawData)


## Long format ----
long <- rawData %>%
  pivot_longer(matches("^(mean|var)_F6_EYT_t[1-5]$"),
               names_to      = c("stat", "trait"),
               names_pattern = "^(mean|var)_F6_EYT_t([1-5])$") %>%
  pivot_wider(names_from = stat, values_from = value) %>%
  mutate(year     = year - nYearsBurnIn,
         trait    = factor(paste0("Trait ", trait)),
         method   = factor(method,   levels = methodLevels),
         scenario = factor(scenario, levels = names(gainScenarios)))


## Mean and SE across replicates ----
summ <- long %>%
  group_by(scenario, method, trait, year) %>%
  summarise(mean_mean = mean(mean),
            mean_se   = sd(mean) / sqrt(n()),
            var_mean  = mean(var),
            var_se    = sd(var) / sqrt(n()),
            .groups = "drop")


## Shared look ----
themeFig <- theme_bw() +
  theme(legend.position  = "bottom",
        plot.title       = element_text(face = "bold", size = 14, hjust = 0.5),
        strip.text.y     = element_text(face = "bold"),
        strip.text.x     = element_text(face = "bold"))
methodColours <- c("DGI" = "#D55E00", "MRGI" = "#0072B2")


## Figure 1: F6 mean genetic value, rows = scenario, columns = trait ----
p_mean <- ggplot(summ, aes(x = year, y = mean_mean, colour = method, fill = method)) +
  geom_ribbon(aes(ymin = mean_mean - mean_se, ymax = mean_mean + mean_se),
              alpha = 0.2, linetype = 0) +
  geom_line(linewidth = 1) +
  facet_grid(scenario ~ trait, scales = "free_y",
             labeller = labeller(scenario = label_wrap_gen(22))) +
  scale_colour_manual(values = methodColours) +
  scale_fill_manual(values = methodColours) +
  scale_x_continuous("Years after burn-in", limits = c(0, 25)) +
  scale_y_continuous("Mean genetic value (F6 - EYT)") +
  labs(title  = "Mean genetic value of the F6 population: DGI vs MRGI",
       colour = "Method", fill = "Method") +
  themeFig
p_mean
ggsave("results/Fig_F6_mean_DGI_vs_MRGI.png", p_mean, width = 11, height = 7, dpi = 300)


## Figure 2: F6 genetic variance, same layout ----
p_var <- ggplot(summ, aes(x = year, y = var_mean, colour = method, fill = method)) +
  geom_ribbon(aes(ymin = var_mean - var_se, ymax = var_mean + var_se),
              alpha = 0.2, linetype = 0) +
  geom_line(linewidth = 1) +
  facet_grid(scenario ~ trait, scales = "free_y",
             labeller = labeller(scenario = label_wrap_gen(22))) +
  scale_colour_manual(values = methodColours) +
  scale_fill_manual(values = methodColours) +
  scale_x_continuous("Years after burn-in", limits = c(0, 25)) +
  scale_y_continuous("Genetic variance (F6 - EYT)", limits = c(0, NA)) +
  labs(title  = "Genetic variance of the F6 population: DGI vs MRGI",
       colour = "Method", fill = "Method") +
  themeFig
p_var
ggsave("results/Fig_F6_variance_DGI_vs_MRGI.png", p_var, width = 11, height = 7, dpi = 300)


## Table: realized vs required selection intensity per step ----
# One row per scenario, method and selection step, from that method's own runs.
# Realized  = intensity actually applied: nInd / N (fixed by design, same for both methods).
# Required  = intensity the method's index needs to deliver the target of the step
#             (1 base genetic SD x shape; DG: exact target, MRGI: minimum target).
# k: mean over the 25 years within each rep, then mean +- SD across reps; the required %
#    is the % to select at that mean k (normal index values, truncation selection).
# Target_delivered_pct: expected gain at the realized intensity / target, worst trait.
# Index_accuracy: correlation between the index and its genetic value, sqrt(b'Gb / b'Pb).
source("scripts/MRGI_functions.R")   # propFromIntensity()
stageLabels <- c(parents = "Parents (F5_YT + F6_EYT)", F7 = "F7 (from F6_EYT)",
                 F6 = "F6 (from F5_YT)", F5 = "F5 (from headrows)")
logs <- list()
for (SCEN in names(gainScenarios)) {
  for (METHOD in names(gainScenarios[[SCEN]])) {
    for (REP in 1:nReps) {
      x <- readRDS(file.path(logDir, paste0(gainScenarios[[SCEN]][[METHOD]],
                                            "_selection_log_", REP, ".rds")))
      logs[[length(logs) + 1]] <- data.frame(
        scenario = SCEN, method = METHOD, rep = REP, stage = x$stage,
        pct_realized = x$pct_applied, k_realized = x$k_applied,
        k_required   = if (METHOD == "MRGI") x$k_required_MRGI else x$k_required_DG,
        delivered    = x$min_ratio_E_over_M,
        accuracy     = x$index_accuracy)
    }
  }
}
logs <- bind_rows(logs)

fmtPct <- function(p) ifelse(p >= 0.01, formatC(p, format = "f", digits = 2),
                      ifelse(p > 0, formatC(p, format = "e", digits = 1), "< 1e-300"))

intensityTable <- logs %>%
  group_by(scenario, method, stage, rep) %>%                  # mean over years per rep
  summarise(across(c(pct_realized, k_realized, k_required, delivered, accuracy), mean),
            .groups = "drop") %>%
  group_by(scenario, method, stage) %>%                       # across reps
  summarise(Realized_pct = mean(pct_realized),
            Realized_k   = mean(k_realized),
            k_req_mean   = mean(k_required),
            k_req_sd     = sd(k_required),
            Target_delivered_pct = round(100 * mean(delivered)),
            Index_accuracy       = round(mean(accuracy), 2),
            .groups = "drop") %>%
  mutate(Required_k   = sprintf("%.2f +- %.2f", k_req_mean, k_req_sd),
         Required_pct = fmtPct(100 * propFromIntensity(k_req_mean)),
         Realized_pct = round(Realized_pct, 2),
         Realized_k   = round(Realized_k, 2),
         stage    = factor(stageLabels[stage], levels = stageLabels),
         scenario = factor(scenario, levels = names(gainScenarios)),
         method   = factor(method,   levels = methodLevels)) %>%
  arrange(scenario, stage, method) %>%
  select(Scenario = scenario, Method = method, Stage = stage,
         Realized_pct, Realized_k, Required_pct, Required_k,
         Target_delivered_pct, Index_accuracy)

print(as.data.frame(intensityTable), row.names = FALSE)
write.csv(intensityTable, "results/Table_realized_vs_required_intensity.csv", row.names = FALSE)
