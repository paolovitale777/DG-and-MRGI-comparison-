

### Functions ###

### Error variance required to obtain target Heritability
#var_e <- function(var_g, H2) {
#  (var_g / H2) - var_g
#}
var_e <- function(var_g, H2) {
  var_g_diag <- diag(var_g)
  var_e_diag <- (var_g_diag / H2) - var_g_diag
  diag(var_e_diag)
}


## Simulation parameters ----
# Models a single trait representing dry yield

# Loop over replicates
# Each replicate models 20 years of breeding program and saves the output

nYearsBurnIn <- 5
yearsTotal <- nYearsBurnIn + 25   # -> 10 future breeding years

### Those values are closer to our recommendations about optimal sampling
### of parents and crosses, while increasing selection intensity AND halving
### the amount of initial plot resources.
nParentsBrnIn <- 60
nCrossesBrnIn <- 100

nParentsTotal <- 60  # NOT IN THE PHENOTYPIC SELECTION SCRIPT - The number of parents is the only parameter changed.
nCrosses <- 100

nDH<- 100
n_Ind_F4<-3000
n_Ind_F5_YT <- 1000
n_Ind_F6_EYT <- 200
n_Ind_F7_EYT <- 50

# Number of stochastic simulation replications
nReps <- 5

# Genetic mean
genMeans <- c(4,7,5,6)
nTraits  <- length(genMeans)   # number of traits, used by the selection scripts

# Total genetic variance
varGens <- c(0.5, 0.8, 0.6,0.2)


# Stage-based H2 from experimental MET data, set in the founders during pipeline filling.
# When these are set, reps in setPheno is set to 1.
# CURRENTLY ONLY EXEMPLARY PLACEHOLDER VALUES.

H2_F5s     <- c(0.05,0.05,0.05, 0.05)
H2_F5_YTs  <- c(0.5,0.6, 0.55, 0.2)
H2_F6_EYTs <- c(0.7, 0.8, 0.75,0.3)
H2_F7_EYTs <- c(0.7, 0.8, 0.75, 0.3)

# trait genetic correlations


cors <- matrix(c( 1.0, 0.8, 0.5, 0.7,
                  0.8, 1.0, 0.6, 0.5,
                  0.5, 0.6, 1.0, 0.7,
                  0.7, 0.5, 0.7, 1.0), nrow = 4, byrow = TRUE)
#

# Residual error for reference plot (e.g., smallest plot size unit)
#plotError <- 2 # H2 = 0.2 (= 1 / (1+4))
# Residual error for plants
#plantError <- 99 # H2 = 0.01, low so that family selection has reasonable accuracy

# Number of chromosomes
nChr <- 21

# Number of QTL, SNP, and segregating sites per chromosome
nQtl <- 300 # 20 chr x 100 QTL = 2000 alleles with an effect on your quantitative trait!
nSnp <- 300 # SNP chip with 6000 markers      # if GS is simulated
nSegSites <- nQtl + nSnp


