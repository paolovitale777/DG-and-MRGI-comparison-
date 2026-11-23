## Table S2: genetic gain after 25 years of index selection and selection intensity per step
# Run from the repository root after scripts 03 to 08. Writes results/Table_S2_results.csv.
# Gain = change in the mean genetic value of stage 2 (F6, elite yield trial) lines from the
# end of burn-in (year 5) to year 30, in genetic SD units of the stage 2 lines at the end of
# burn-in (and in trait units). Required k = intensity each index needs to deliver its target
# in one selection step (mean over 25 years within a replicate, then mean and SD across reps).

source("scripts/MRGI_functions.R")   # propFromIntensity()
nReps   <- 5
outDir  <- "output"
yBurnIn <- 5
yEnd    <- 30

runs <- list(
  list(sc = "1", m = "DGI",  f = "Multitrait_DG_index_same_weights"),
  list(sc = "1", m = "MRGI", f = "Multitrait_MRGI_index_same_weights"),
  list(sc = "2", m = "DGI",  f = "Multitrait_DG_index_different_weights"),
  list(sc = "2", m = "MRGI", f = "Multitrait_MRGI_index_uneven_weights"),
  list(sc = "3", m = "DGI",  f = "Multitrait_DG_index_realistic"),
  list(sc = "3", m = "MRGI", f = "Multitrait_MRGI_index_realistic"))

out <- list()
for (r in runs) {
  res <- lapply(1:nReps, function(i) readRDS(file.path(outDir, paste0(r$f, "_", i, ".rds"))))
  g  <- sapply(res, function(o) sapply(1:4, function(t)
          (o[[paste0("mean_F6_EYT_t", t)]][yEnd] - o[[paste0("mean_F6_EYT_t", t)]][yBurnIn]) /
           sqrt(o[[paste0("var_F6_EYT_t", t)]][yBurnIn])))
  gu <- sapply(res, function(o) sapply(1:4, function(t)
          o[[paste0("mean_F6_EYT_t", t)]][yEnd] - o[[paste0("mean_F6_EYT_t", t)]][yBurnIn]))
  L <- do.call(rbind, lapply(1:nReps, function(i)
         cbind(rep = i, readRDS(file.path(outDir, paste0(r$f, "_selection_log_", i, ".rds"))))))
  own <- if (r$m == "MRGI") L$k_required_MRGI else L$k_required_DG
  kr  <- tapply(own, list(L$stage, L$rep), mean)
  ka  <- tapply(L$k_applied, L$stage, mean)
  pa  <- tapply(L$pct_applied, L$stage, mean)
  for (st in rownames(kr)) {
    km <- mean(kr[st, ])
    out[[length(out) + 1]] <- data.frame(
      scenario = r$sc, index = r$m, stage = st,
      pct_applied = pa[st], k_applied = ka[st],
      k_required = km, k_required_sd = sd(kr[st, ]), pct_required = 100 * propFromIntensity(km),
      gain_SD_t1 = mean(g[1, ]), gain_SD_t1_sd = sd(g[1, ]),
      gain_SD_t2 = mean(g[2, ]), gain_SD_t2_sd = sd(g[2, ]),
      gain_SD_t3 = mean(g[3, ]), gain_SD_t3_sd = sd(g[3, ]),
      gain_SD_t4 = mean(g[4, ]), gain_SD_t4_sd = sd(g[4, ]),
      gain_units_t1 = mean(gu[1, ]), gain_units_t2 = mean(gu[2, ]),
      gain_units_t3 = mean(gu[3, ]), gain_units_t4 = mean(gu[4, ]))
  }
}
tab <- do.call(rbind, out); rownames(tab) <- NULL
dir.create("results", showWarnings = FALSE)
write.csv(tab, "results/Table_S2_results.csv", row.names = FALSE)
print(tab, digits = 3)
