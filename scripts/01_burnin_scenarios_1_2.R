
rm(list=ls())   # clear environment

####### BurnIn #######

### LINE BREEDING PROGRAM EXAMPLE in AlphaSimR (c.werner@cgiar.org)

library(AlphaSimR)
dir.create("output", showWarnings = FALSE)   # simulation outputs go here

## Simulation parameters ----
# Models a single trait representing dry yield

source("scripts/sim_param_scenarios_1_2.R")


# REP <- 1

for(REP in 1:nReps){
  output <- data.frame(year = 1:yearsTotal,
                       rep = rep(REP, yearsTotal),
                       scenario = rep("BURNIN", yearsTotal),
                       mean_parents_t1 = numeric(yearsTotal),
                       mean_parents_t2 = numeric(yearsTotal),
                       mean_parents_t3 = numeric(yearsTotal),
                       mean_parents_t4 = numeric(yearsTotal),
                       var_parents_t1 = numeric(yearsTotal),
                       var_parents_t2 = numeric(yearsTotal),
                       var_parents_t3 = numeric(yearsTotal),
                       var_parents_t4 = numeric(yearsTotal),
                       genic_parents_t1 = numeric(yearsTotal),
                       genic_parents_t2 = numeric(yearsTotal),
                       genic_parents_t3 = numeric(yearsTotal),
                       genic_parents_t4 = numeric(yearsTotal),
                       mean_F5_YT_t1 = numeric(yearsTotal),
                       mean_F5_YT_t2 = numeric(yearsTotal),
                       mean_F5_YT_t3 = numeric(yearsTotal),
                       mean_F5_YT_t4 = numeric(yearsTotal),
                       var_F5_YT_t1 = numeric(yearsTotal),
                       var_F5_YT_t2 = numeric(yearsTotal),
                       var_F5_YT_t3 = numeric(yearsTotal),
                       var_F5_YT_t4 = numeric(yearsTotal),
                       genic_F5_YT_t1 = numeric(yearsTotal),
                       genic_F5_YT_t2 = numeric(yearsTotal),
                       genic_F5_YT_t3 = numeric(yearsTotal),
                       genic_F5_YT_t4 = numeric(yearsTotal),
                       mean_F7_AdvEYT_t1 = numeric(yearsTotal),
                       mean_F7_AdvEYT_t2 = numeric(yearsTotal),
                       mean_F7_AdvEYT_t3 = numeric(yearsTotal),
                       mean_F7_AdvEYT_t4 = numeric(yearsTotal),
                       var_F7_AdvEYT_t1 = numeric(yearsTotal),
                       var_F7_AdvEYT_t2 = numeric(yearsTotal),
                       var_F7_AdvEYT_t3 = numeric(yearsTotal),
                       var_F7_AdvEYT_t4 = numeric(yearsTotal),
                       genic_F7_AdvEYT_t1 = numeric(yearsTotal),
                       genic_F7_AdvEYT_t2 = numeric(yearsTotal),
                       genic_F7_AdvEYT_t3 = numeric(yearsTotal),
                       genic_F7_AdvEYT_t4 = numeric(yearsTotal),
                       stringsAsFactors=FALSE)


  ## Alternative approach: randomly sampled alleles to create an unrelated founder population.
  founderPop <- quickHaplo(nInd = nParentsBrnIn,
                           nChr = nChr,
                           segSites = nSegSites)

  SP <- SimParam$new(founderPop)

  SP$addTraitA(nQtlPerChr = nQtl,   # Additive trait + simple GxE interaction effect
                mean = genMeans,
                var = varGens,
               corA= cors)

  # add SNP chip
  SP$addSnpChip(nSnpPerChr = nSnp)         # if GS is simulated.

  # turn "founderPop" with all its characteristics
  # into an actual population object that can be used for simulation
  Parents <- newPop(founderPop)

  ####################################################################################
  ## 1. FILL THE BREEDING PROGRAM PIPELINE ---- ######################################
  ####################################################################################

  for(cycle in 1:6){

    # Crossing block (yr1)
    F1 <- randCross(Parents,
                    nCrosses = nCrossesBrnIn,
                    nProgeny = 1)


    # F2 - (yr1)
    if(cycle < 6){

      #Double Haploids
      DH=makeDH(F1,nDH)# only 100 needed under current simulation framework.

    }

    
    if(cycle < 5){
      
      #Headrows
      F4 <- selectInd(DH, n_Ind_F4, use="rand")
    }
    
    if(cycle < 4){
      F5_YT <- selectInd(F4, n_Ind_F5_YT, use="rand")
      F5_YT <- self(F5_YT)

      if(cycle == 1){          # setting var_e in first cycle
        var_g_F5_YT <- varG(F5_YT)
        var_e_F5_YT <- var_e(var_g_F5_YT, H2_F5_YTs)
      }
      F5_YT <- setPheno(F5_YT,
                        varE = var_e_F5_YT,
                        reps = 1)

    }
   # F6: Preliminary Yield Trials - (yr4)
    if(cycle < 3){
      F6_EYT <- selectInd(F5_YT, n_Ind_F6_EYT, use="rand")
      F6_EYT <- self(F6_EYT)

      if(cycle == 1){
        var_g_F6_EYT <- varG(F6_EYT)
        var_e_F6_EYT <- var_e(var_g_F6_EYT, H2_F6_EYTs)
      }
      diag(var_g_F6_EYT) / (diag(var_g_F6_EYT) + diag(var_e_F6_EYT))

      F6_EYT <- setPheno(F6_EYT,
                         varE = var_e_F6_EYT,
                         reps = 1)  # 6*2, # 4 locs, 2 reps, 4x plot size
      }

    # F7: Advanced Yield Trials -(yr5)
    if(cycle < 2){

      F7_AdvEYT <- selectInd(F6_EYT, n_Ind_F7_EYT, use="rand")
      F7_AdvEYT <- self(F7_AdvEYT)

      if(cycle == 1){
        var_g_F7_EYT <- varG(F7_AdvEYT)
        var_e_F7_EYT <- var_e(var_g_F7_EYT, H2_F7_EYTs)
      }
      diag(var_g_F7_EYT) / (diag(var_g_F7_EYT) + diag(var_e_F7_EYT))

      F7_AdvEYT <- setPheno(F7_AdvEYT,
                            varE = var_e_F7_EYT,
                            reps = 1)  # 6*2, # 3 locs, 2 reps, 12x plot size
      }

    print(paste0("REP ", REP, " | Fill Pipeline cycle ", cycle))
    #print(paste0("REP ", REP, " | Burn-in year ", year, " of ", nYearsBurnIn))  ??? There is no year, this must have given you an error.
  }


  all_stages <- vector("list", yearsTotal)

  P <- runif(yearsTotal)

  ## Run 10 years of historic (burn-in) breeding program ----
  for(year in 1:nYearsBurnIn){

    print(paste0("Burn-in year ", year))

    parental_lines<-mergePops(list(F5_YT, F6_EYT))
    Parents <- selectInd(pop = parental_lines,
                         nInd=nParentsTotal,
                         use="rand")


    # F7: Advanced Elite Yield Trial
    F7_AdvEYT <- selectInd(F6_EYT, n_Ind_F7_EYT, use="rand")
    F7_AdvEYT <- self(F7_AdvEYT)
    F7_AdvEYT <- setPheno(F7_AdvEYT,
                          varE = var_e_F7_EYT,
                          reps = 1)  #6*2,



    # F6: Elite Yield Trials
    F6_EYT <- selectInd(F5_YT, n_Ind_F6_EYT, use="rand")
    F6_EYT <- self(F6_EYT)
    F6_EYT <- setPheno(F6_EYT,
                       varE = var_e_F6_EYT,
                       reps = 1)  # 6*2, # 4 locs, 2 reps, 1x plot size


    # F5 - Obregon - season 6 (yr3)
    F5_YT <- selectInd(F4, n_Ind_F5_YT, use="rand")
    F5_YT <- setPheno(F5_YT,
                      varE = var_e_F5_YT,
                      reps = 1)  #4*1, # 1 plot

    
    #Headrows
    F4 <- selectInd(DH, n_Ind_F4, use="rand")
    
    
    DH=makeDH(F1,nDH)# only 100 needed under current simulation framework.
    # Crossing block (yr1)
    F1 <- randCross(Parents,
                    nCrosses = nCrossesBrnIn,
                    nProgeny = 1)


    # TRAINING POPULATION FOR GS: USE ALL AVAILABLE PHENOTYPIC DATA
    # ON GENOTYPED INDIVIDUALS IN THE TRAINING POPULATION
     all_stages[[year]] <- mergePops(list(F5_YT, F6_EYT, F7_AdvEYT))


     # Report performance of new F1s
     output$mean_parents_t1[year] <- meanG(Parents)[1]
     output$mean_parents_t2[year] <- meanG(Parents)[2]
     output$mean_parents_t3[year] <- meanG(Parents)[3]
     output$mean_parents_t4[year] <- meanG(Parents)[4]
     output$var_parents_t1[year] <- varG(Parents)[1]
     output$var_parents_t2[year] <- varG(Parents)[2,2]
     output$var_parents_t3[year] <- varG(Parents)[3,3]
     output$var_parents_t4[year] <- varG(Parents)[4,4]

     # Report performance of new F8s
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
  save.image(paste0("output/Multitrait_burnIn", REP, ".RData"))

}

head(output)

rm(list=ls())




