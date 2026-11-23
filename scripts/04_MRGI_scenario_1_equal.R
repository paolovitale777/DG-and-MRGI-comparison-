
### Minimum Required Gain Index (MRGI), fixed selection intensity ###
# Same structure as the desired gain (DG) scenario; only the index changes:
#   DG  : b = G^-1 d                         -> E proportional to d (exact pattern)
#   MRGI: b = argmin c'Pc  s.t.  G c >= M    -> E >= M, lowest intensity (Itoh & Yamada 1988)
# Here the number selected per stage is fixed, so only the SHAPE of minGains matters
# (scaling M scales b, and the ranking of lines is unchanged).

library(AlphaSimR)
dir.create("output", showWarnings = FALSE)   # simulation outputs go here

# Give your scenario a name
scenarioName <- "Multitrait_MRGI_index_same_weights"
source("scripts/sim_param_scenarios_1_2.R")
source("scripts/MRGI_functions.R")   # mrgiWeights()

# Minimum required gains (relative, all traits to be increased; equal gains, 1,1,1,1)
minGains   <- rep(1, nTraits)
useSDunits <- TRUE           # TRUE  = M is in genetic SD units of the base pop (end of burn-in)
                             # FALSE = M is in trait units

# MRGI phenotypic weights for one stage (G, P of the candidates)
mrgiStageWeights <- function(G, P, M) {
  mrgiWeights(as.matrix(G), as.matrix(P), M, compareDGI = FALSE)$b
}

# Loop over replicates
for(REP in 1:nReps){

  ## Load data ----
  load(paste0("output/Multitrait_burnIn",REP,".RData"))

  # Change scenario
  output$scenario <- scenarioName

  # Same seed for all four scenarios at a given rep: reproducible runs, and DG and MRGI
  # start each rep from the same burn-in AND the same random stream (paired comparison)
  set.seed(REP)

  # Breeding objective, fixed for the whole run and the same at every stage
  G0 <- varG(mergePops(list(F5_YT, F6_EYT)))
  M  <- if(useSDunits) minGains * sqrt(diag(G0)) else minGains

  # Headrows were never phenotyped in the burn-in: residual variance from H2_F5s
  var_e_F4 <- var_e(varG(F4), H2_F5s)
  F4 <- setPheno(F4, varE = var_e_F4, reps = 1)

  selLog <- list()   # one row per selection step and year

  # F6 (Elite Yield Trial) statistics: columns + value at the end of burn-in
  recordF6 <- function(output, pop, yr) {
    for(t in 1:nTraits){
      output[[paste0("mean_F6_EYT_t",  t)]][yr] <- meanG(pop)[t]
      output[[paste0("var_F6_EYT_t",   t)]][yr] <- varG(pop)[t, t]
      output[[paste0("genic_F6_EYT_t", t)]][yr] <- genicVarG(pop)[t]
    }
    output
  }
  for(t in 1:nTraits) for(st in c("mean", "var", "genic"))
    output[[paste0(st, "_F6_EYT_t", t)]] <- NA_real_
  output <- recordF6(output, F6_EYT, nYearsBurnIn)

  ## Run the breeding program for (yearsTotal - nYearsBurnIn) = 25 years ----
  for(year in (nYearsBurnIn + 1):yearsTotal){

    print(paste0(scenarioName, " year ", year))
    # select new parents in current year

    # select new parents in current year
    # nParents <- sample(c(100:120), 1) #
    parental_lines <- mergePops(list(F5_YT, F6_EYT))
    # candidates from two stages with different error variances: observed P
    b <- mrgiStageWeights(varG(parental_lines), var(parental_lines@pheno), M)
    selLog[[length(selLog) + 1]] <- cbind(rep = REP, year = year,
      stepInfo("MRGI", "parents", N = parental_lines@nInd, nInd = nParentsTotal, b = b,
               G = varG(parental_lines), P = var(parental_lines@pheno), M = M))

    parent_selected <- selectInd(pop   = parental_lines,
                                 nInd  = nParentsTotal,
                                 trait = selIndex,
                                 b     = b)


    # MRGI weights with known (simulated) G and P = G + E of the stage
    b <- mrgiStageWeights(varG(F6_EYT), varG(F6_EYT) + var_e_F6_EYT, M)
    selLog[[length(selLog) + 1]] <- cbind(rep = REP, year = year,
      stepInfo("MRGI", "F7", N = F6_EYT@nInd, nInd = n_Ind_F7_EYT, b = b,
               G = varG(F6_EYT), P = varG(F6_EYT) + var_e_F6_EYT, M = M))

    F7_AdvEYT <- selectInd(pop = F6_EYT,
                           nInd=n_Ind_F7_EYT,
                           trait = selIndex,
                           b=b)

    # F7: Advanced Elite Yield Trial
    F7_AdvEYT <- self(F7_AdvEYT)
    F7_AdvEYT <- setPheno(F7_AdvEYT,
                          varE = var_e_F7_EYT,
                          reps = 1)  #6*2,




    # MRGI weights with known (simulated) G and P = G + E of the stage
    b <- mrgiStageWeights(varG(F5_YT), varG(F5_YT) + var_e_F5_YT, M)
    selLog[[length(selLog) + 1]] <- cbind(rep = REP, year = year,
      stepInfo("MRGI", "F6", N = F5_YT@nInd, nInd = n_Ind_F6_EYT, b = b,
               G = varG(F5_YT), P = varG(F5_YT) + var_e_F5_YT, M = M))


    F6_EYT <- selectInd(pop = F5_YT,
                        nInd=n_Ind_F6_EYT,
                        trait = selIndex,
                        b=b)


    # F6: Elite Yield Trials
    F6_EYT <- self(F6_EYT)
    F6_EYT <- setPheno(F6_EYT,
                       varE = var_e_F6_EYT,
                       reps = 1)  # 6*2, # 4 locs, 2 reps, 1x plot size


    # F5 - Obregon - season 6 (yr3)
    b <- mrgiStageWeights(varG(F4), varG(F4) + var_e_F4, M)
    selLog[[length(selLog) + 1]] <- cbind(rep = REP, year = year,
      stepInfo("MRGI", "F5", N = F4@nInd, nInd = n_Ind_F5_YT, b = b,
               G = varG(F4), P = varG(F4) + var_e_F4, M = M))
    F5_YT <- selectInd(pop = F4,
                       nInd = n_Ind_F5_YT,
                       trait = selIndex,
                       b = b)
    F5_YT <- setPheno(F5_YT,
                      varE = var_e_F5_YT,
                      reps = 1)  #4*1, # 1 plot

    
    #Headrows
    F4 <- selectInd(DH, n_Ind_F4, use="rand")
    F4 <- setPheno(F4, varE = var_e_F4, reps = 1)
    
    
    DH=makeDH(F1,nDH)# only 100 needed under current simulation framework.
    
    # Crossing block (yr1)
    F1 <- randCross(parent_selected,
                    nCrosses = nCrossesBrnIn,
                    nProgeny = 1)

    meanG(parent_selected)
    # Report performance of new F1s
    output$mean_parents_t1[year] <- meanG(parent_selected)[1]
    output$mean_parents_t2[year] <- meanG(parent_selected)[2]
    output$mean_parents_t3[year] <- meanG(parent_selected)[3]
    output$mean_parents_t4[year] <- meanG(parent_selected)[4]
    output$var_parents_t1[year] <- varG(parent_selected)[1]
    output$var_parents_t2[year] <- varG(parent_selected)[2,2]
    output$var_parents_t3[year] <- varG(parent_selected)[3,3]
    output$var_parents_t4[year] <- varG(parent_selected)[4,4]

    # Report performance of new F8s
    output <- recordF6(output, F6_EYT, year)   # F6 lines of this year

    output$mean_F5_YT_t1[year] <- meanG(F5_YT)[1]
    output$mean_F5_YT_t2[year] <- meanG(F5_YT)[2]
    output$mean_F5_YT_t3[year] <- meanG(F5_YT)[3]
    output$mean_F5_YT_t4[year] <- meanG(F5_YT)[4]
    output$var_F5_YT_t1[year] <- varG(F5_YT)[1]
    output$var_F5_YT_t2[year] <- varG(F5_YT)[2,2]
    output$var_F5_YT_t3[year] <- varG(F5_YT)[3,3]
    output$var_F5_YT_t4[year] <- varG(F5_YT)[4,4]
    output$genic_F5_YT_t1[year] <- genicVarG(F5_YT)[1]
    output$genic_F5_YT_t2[year] <- genicVarG(F5_YT)[2]
    output$genic_F5_YT_t3[year] <- genicVarG(F5_YT)[3]
    output$genic_F5_YT_t4[year] <- genicVarG(F5_YT)[4]


    output$mean_F7_AdvEYT_t1[year] <- meanG(F7_AdvEYT)[1]
    output$mean_F7_AdvEYT_t2[year] <- meanG(F7_AdvEYT)[2]
    output$mean_F7_AdvEYT_t3[year] <- meanG(F7_AdvEYT)[3]
    output$mean_F7_AdvEYT_t4[year] <- meanG(F7_AdvEYT)[4]
    output$var_F7_AdvEYT_t1[year] <- varG(F7_AdvEYT)[1]
    output$var_F7_AdvEYT_t2[year] <- varG(F7_AdvEYT)[2,2]
    output$var_F7_AdvEYT_t3[year] <- varG(F7_AdvEYT)[3,3]
    output$var_F7_AdvEYT_t4[year] <- varG(F7_AdvEYT)[4,4]

    output$genic_F7_AdvEYT_t1[year] <- genicVarG(F7_AdvEYT)[1]
    output$genic_F7_AdvEYT_t2[year] <- genicVarG(F7_AdvEYT)[2]
    output$genic_F7_AdvEYT_t3[year] <- genicVarG(F7_AdvEYT)[3]
    output$genic_F7_AdvEYT_t4[year] <- genicVarG(F7_AdvEYT)[4]


  }

  # Save output
  cat(paste0("saving REP ", REP, "\n"))
  saveRDS(output,paste0("output/", scenarioName,"_",REP,".rds"))
  saveRDS(do.call(rbind, selLog), paste0("output/", scenarioName, "_selection_log_", REP, ".rds"))
}


# Record software versions used for these runs (AlphaSimR version changes the results)
writeLines(capture.output(sessionInfo()), paste0("output/", scenarioName, "_sessionInfo.txt"))

head(output)

rm(list=ls())

ls()
