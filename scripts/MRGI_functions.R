
### Minimum Required Gain Index (MRGI) ###
# Martini, Vitale, Ceron-Rojas, Montesinos-Lopez, Crossa, Ortiz
# "Minimum Required Gains: A Practical Reframing of Desired-Gain Selection Indices"
#
# Model: y = mu + g + e, g ~ N(0, G), e ~ N(0, E), P = G + E
# Expected response of an LPSI with weights b at intensity k (Eq. 3):
#   E_{b,k} = k G b / sqrt(b' P b)
#
# DGI  (Pesek & Baker 1969, Proposition 1): b_d = G^-1 E_d,  k_d = sqrt(b_d' P b_d)
# MRGI (Itoh & Yamada 1988):                b_opt = argmin_c c' P c  s.t.  G c >= M
#                                           k_opt = sqrt(b_opt' P b_opt), E_opt = G b_opt
#
# All gains are lower bounds in the "more is better" direction. For a trait where a
# decrease is desired, flip the sign of that trait (use `direction` below).
#
# Source this file; it only defines functions. Requires: quadprog (AlphaSimR only for
# selectMRGI()).

library(quadprog)


## Selection intensity <-> selected proportion (normal index, truncation selection) ----

# Selection intensity k for a selected proportion p
intensityFromProp <- function(p) {
  stopifnot(all(p > 0 & p < 1))
  dnorm(qnorm(1 - p)) / p
}

# Selected proportion p that gives intensity k (inverse of intensityFromProp)
propFromIntensity <- function(k) {
  # k(x) = dnorm(x) / (1 - pnorm(x)) for truncation point x; solved on the log scale so
  # that very high intensities (tiny p, e.g. k = 16) still work
  kOfX <- function(x) exp(dnorm(x, log = TRUE) - pnorm(x, lower.tail = FALSE, log.p = TRUE))
  sapply(k, function(kk) {
    if (kk <= 0) return(1)
    x <- uniroot(function(x) kOfX(x) - kk, lower = -40, upper = kk + 1, tol = 1e-12)$root
    pnorm(x, lower.tail = FALSE)
  })
}


## Internal helpers ----

.checkInputs <- function(G, P, M) {
  G <- as.matrix(G); P <- as.matrix(P); M <- as.numeric(M)
  t <- length(M)
  if (!all(dim(G) == t) || !all(dim(P) == t))
    stop("G and P must be t x t, with t = length(M)")
  if (!isSymmetric(unname(G), tol = 1e-8) || !isSymmetric(unname(P), tol = 1e-8))
    stop("G and P must be symmetric")
  if (min(eigen(P, symmetric = TRUE, only.values = TRUE)$values) <= 0)
    stop("P must be positive definite")
  if (all(M == 0)) stop("M must not be the zero vector")
  list(G = G, P = P, M = M)
}

# Expected response for weights b at intensity k (Eq. 3)
expectedGain <- function(b, G, P, k = 1) {
  as.vector(k * G %*% b / sqrt(drop(t(b) %*% P %*% b)))
}


## Desired gain index (exact target) ----

dgiWeights <- function(G, P, Ed, direction = NULL) {
  x <- .checkInputs(G, P, Ed)
  s <- if (is.null(direction)) rep(1, length(x$M)) else sign(direction)
  S <- diag(s, length(s))
  Gs <- S %*% x$G %*% S; Ps <- S %*% x$P %*% S

  b <- as.vector(solve(Gs, x$M))
  k <- sqrt(drop(t(b) %*% Ps %*% b))
  b <- s * b                                   # back to original trait scale
  list(method = "DGI", b = b, k = k, p = propFromIntensity(k),
       E = expectedGain(b, x$G, x$P, k))
}


## Minimum required gain index ----

#' @param G         t x t genetic covariance matrix
#' @param P         t x t phenotypic covariance matrix (G + E)
#' @param M         length-t vector of minimum required gains (>= 0, in the desired direction)
#' @param direction optional length-t vector of +1 / -1: -1 means a decrease is desired for that
#'                  trait (M is then given as the absolute minimum decrease)
#' @param compareDGI also return the exact-target DGI with E_d = M
#' @return list with b (phenotypic weights, original trait scale), k (required selection
#'         intensity), p (selected proportion under truncation), E (expected gain at k),
#'         binding (traits whose constraint is active, i.e. E == M), and optionally dgi
# maintain : optional logical vector; TRUE = keep that trait unchanged (expected gain = 0,
#            an equality constraint as in a restricted index). Its M must be 0 and its
#            direction is irrelevant. FALSE = minimum gain M (inequality, E >= M).
mrgiWeights <- function(G, P, M, direction = NULL, compareDGI = TRUE, maintain = NULL) {
  x <- .checkInputs(G, P, M)
  t <- length(x$M)
  eq <- if (is.null(maintain)) rep(FALSE, t) else as.logical(maintain)
  if (length(eq) != t) stop("maintain must be a logical vector of length t")
  if (any(x$M[eq] != 0)) stop("traits to maintain must have M = 0")
  s <- if (is.null(direction)) rep(1, t) else sign(direction)
  if (length(s) != t || any(s == 0)) stop("direction must be a vector of +1/-1, length t")
  S <- diag(s, t)
  Gs <- S %*% x$G %*% S; Ps <- S %*% x$P %*% S

  # minimize c'Pc  s.t.  (G c)_i = 0 for maintained traits,  (G c)_i >= M_i otherwise
  # solve.QP: min 1/2 c'Dc - d'c  s.t.  A'c >= b0 (first meq rows as equalities)
  #           -> D = 2P, d = 0, A = t(G), equality constraints placed first
  ord <- c(which(eq), which(!eq))
  qp <- solve.QP(Dmat = 2 * Ps, dvec = rep(0, t), Amat = t(Gs)[, ord, drop = FALSE],
                 bvec = x$M[ord], meq = sum(eq))

  b <- qp$solution
  k <- sqrt(drop(t(b) %*% Ps %*% b))           # since G b_opt >= M, this k delivers E >= M
  Eopt_s <- as.vector(Gs %*% b)
  b <- s * b

  out <- list(method = "MRGI", b = b, k = k, p = propFromIntensity(k),
              E = s * Eopt_s,
              binding = seq_len(t) %in% ord[qp$iact],
              maintained = eq,
              M = s * x$M)
  if (compareDGI) out$dgi <- dgiWeights(x$G, x$P, x$M, direction = s)
  class(out) <- c("mrgi", "list")
  out
}

print.mrgi <- function(x, digits = 3, ...) {
  tab <- data.frame(M = x$M, b_MRGI = x$b, E_MRGI = x$E)
  if (!is.null(x$dgi)) { tab$b_DGI <- x$dgi$b; tab$E_DGI <- x$dgi$E }
  tab <- round(tab, digits)
  tab$binding <- x$binding
  print(tab)
  cat(sprintf("\nMRGI: k = %.3f  (select best %.3g %%)\n", x$k, 100 * x$p))
  if (!is.null(x$dgi))
    cat(sprintf("DGI : k = %.3f  (select best %.3g %%)\n", x$dgi$k, 100 * x$dgi$p))
  invisible(x)
}

# Candidates needed to keep m selected individuals at intensity k (N = m / p)
candidatesNeeded <- function(m, k) ceiling(m / propFromIntensity(k))


## AlphaSimR wrapper ----

# Selects nInd individuals from pop on the MRGI. The weights only set the direction of
# the index; the intensity actually applied is the one implied by nInd / pop@nInd.
# The function reports whether that intensity reaches the required k_opt, i.e. whether
# the minimum gains are expected to be met in this selection step.
#   G    : genetic covariance (default: true varG(pop), as in the DG scenario)
#   P    : phenotypic covariance (default: G + varE if varE given, else var(pop@pheno))
#   use  : "pheno" (LPSI, as in the manuscript) or "ebv"
#   returnInfo : TRUE returns list(pop, info) instead of the selected pop only
selectMRGI <- function(pop, nInd, M, G = NULL, P = NULL, varE = NULL,
                       direction = NULL, use = "pheno", verbose = TRUE,
                       returnInfo = FALSE, simParam = NULL) {
  if (is.null(simParam)) simParam <- get("SP", envir = .GlobalEnv)
  if (is.null(G)) G <- varG(pop)
  if (is.null(P)) P <- if (!is.null(varE)) G + varE else var(pop@pheno)

  w <- mrgiWeights(G, P, M, direction = direction, compareDGI = FALSE)
  kApplied <- intensityFromProp(nInd / pop@nInd)

  if (verbose) {
    cat(sprintf("MRGI: required k = %.2f, applied k = %.2f (%d of %d) -> %s\n",
                w$k, kApplied, nInd, pop@nInd,
                ifelse(kApplied >= w$k, "minimum gains expected to be met",
                       "minimum gains NOT reachable at this intensity")))
  }

  sel <- selectInd(pop, nInd = nInd, trait = selIndex, b = w$b, use = use,
                   simParam = simParam)
  if (!returnInfo) return(sel)
  list(pop = sel,
       info = c(w, list(kApplied = kApplied,
                        Eapplied = expectedGain(w$b, G, P, kApplied))))
}


## Index-dependent intensity (as in the manuscript) ----

# The index sets the required intensity k for the absolute gains M of this selection step,
# and the number kept follows from it: nInd = N * p(k), bounded by what the program can
# test (nMax) and a minimum population size (nMin).
#   method : "MRGI" (E >= M, minimum k) or "DGI" (E = M exactly, Pesek & Baker)
#   M      : absolute minimum (MRGI) or exact (DGI) gains for this selection step
#   G, P   : genetic and phenotypic covariance of the candidates
# Returns list(pop = selected, info = one-row data.frame) for reporting.
selectIndexDependent <- function(pop, M, method = c("MRGI", "DGI"), G = NULL, P = NULL,
                                 varE = NULL, nMax = pop@nInd, nMin = 2,
                                 direction = NULL, use = "pheno", simParam = NULL) {
  method <- match.arg(method)
  if (is.null(simParam)) simParam <- get("SP", envir = .GlobalEnv)
  if (is.null(G)) G <- varG(pop)
  if (is.null(P)) P <- if (!is.null(varE)) G + varE else var(pop@pheno)

  w <- if (method == "MRGI") mrgiWeights(G, P, M, direction = direction, compareDGI = FALSE)
       else dgiWeights(G, P, M, direction = direction)

  N      <- pop@nInd
  nWant  <- N * w$p
  nMax   <- min(nMax, N)
  nInd   <- max(nMin, min(nMax, round(nWant)))
  status <- if (round(nWant) > nMax) "capped_nMax"         # index wants more than resources
            else if (round(nWant) < nMin) "floored_nMin"   # index wants fewer than allowed
            else "index"
  kApplied <- if (nInd < N) intensityFromProp(nInd / N) else 0

  sel <- selectInd(pop, nInd = nInd, trait = selIndex, b = w$b, use = use,
                   simParam = simParam)

  s <- if (is.null(direction)) rep(1, length(M)) else sign(direction)
  Eexp <- if (kApplied > 0) expectedGain(w$b, G, P, kApplied) else rep(0, length(M))
  info <- data.frame(method = method, N = N, k_required = w$k, p_required = w$p,
                     nInd = nInd, k_applied = kApplied, status = status,
                     min_ratio_E_over_M = min(s * Eexp / M))
  list(pop = sel, info = info)
}


## Selection-step log for fixed-size scenarios ----

# One row describing a selection step: how many candidates (N) and selected (nInd),
# the intensity applied (from nInd / N), and the intensity the index would NEED to
# deliver the target M in this step, plus the expected gain at the applied intensity.
#   b : index weights used;  G, P : genetic / phenotypic covariance of the candidates
#   M : target of the step (DG: desired gains d; MRGI: minimum gains M)
# k_required = smallest k with k * G b / sd(I) >= M  (for the DG this is k_d of Prop. 1,
#              for the MRGI it is k_opt); it scales linearly with M, so the DG/MRGI
#              ratio of k_required does not depend on the units chosen for M.
stepInfo <- function(method, stage, N, nInd, b, G, P, M, direction = NULL, maintain = NULL) {
  # M         : target magnitudes (>= 0) in the desired direction of each trait
  # direction : +1 = increase wanted, -1 = decrease wanted (default all +1)
  #             A trait with M = 0 must not move in the wrong direction.
  G <- as.matrix(G); P <- as.matrix(P)
  s  <- if (is.null(direction)) rep(1, length(M)) else sign(direction)
  M  <- abs(M)
  r  <- as.vector(G %*% b) / sqrt(drop(t(b) %*% P %*% b))  # gain per unit of intensity
  rs <- s * r                                               # in the desired direction
  pos <- M > 0; zero <- M == 0
  tol <- 1e-8 * max(abs(rs))
  feasible <- all(rs[pos] > 0) && all(rs[zero] >= -tol)
  kReq <- if (feasible) max(M[pos] / rs[pos]) else Inf
  kApp <- if (nInd < N) intensityFromProp(nInd / N) else 0
  E <- kApp * r
  out <- data.frame(method = method, stage = stage, N = N, nInd = nInd,
                    p_applied = nInd / N, k_applied = kApp,
                    k_required = kReq,
                    p_required = if (is.finite(kReq)) propFromIntensity(kReq) else 0,
                    index_accuracy = sqrt(drop(t(b) %*% G %*% b) / drop(t(b) %*% P %*% b)),
                    min_ratio_E_over_M = min(kApp * rs[pos] / M[pos]))
  for (t in seq_along(E)) out[[paste0("E_t", t)]] <- E[t]

  # Both indices on the SAME candidates and target: intensity and % selected each would
  # need to deliver M in this step (DG: E = M exactly; MRGI: E >= M)
  if (any(M > 0)) {
    kDG   <- dgiWeights(G, P, M, direction = s)$k
    kMRGI <- mrgiWeights(G, P, M, direction = s, compareDGI = FALSE, maintain = maintain)$k
    out$k_required_DG     <- kDG
    out$k_required_MRGI   <- kMRGI
    out$pct_required_DG   <- 100 * propFromIntensity(kDG)
    out$pct_required_MRGI <- 100 * propFromIntensity(kMRGI)
  }
  out$pct_applied <- 100 * nInd / N
  out
}
