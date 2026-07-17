# R/genetics.R
#
# Genetic architecture of myrtle rust resistance. Genotypes are stored
# OUTSIDE the main `pop` data frame, in a parallel integer matrix
# `resist_gt` (n_individuals x n_loci, values 0/1/2 = allele dose),
# because the main row-based PlantPopGenFit precedent (jasongbragg/
# PlantPopGenFit, epi branch) keeps demography and genotype in separate
# tables for exactly this reason: n_loci changes shouldn't change the
# shape of the demography table.
#
# Row i of resist_gt always corresponds to row i of pop. Both grow in
# lockstep on every rbind() of new recruits.
#
# Functions here are the only place n_loci is looped over explicitly;
# everything else in the model consumes the single derived
# `resist_score` column and never looks at genotypes directly.

#' Initialise genotypes for a founding population.
#' Loci are independent (Hardy-Weinberg within locus, linkage
#' equilibrium across loci) unless future data justifies otherwise.
init_resist_gt <- function(n_ind, params) {
  n_loci <- length(params$resist_locus_effect)
  gt <- matrix(0L, nrow = n_ind, ncol = n_loci)

  for (j in seq_len(n_loci)) {
    a1 <- rbinom(n_ind, 1, params$resist_freq0[j])
    a2 <- rbinom(n_ind, 1, params$resist_freq0[j])
    gt[, j] <- a1 + a2
  }

  colnames(gt) <- paste0("locus", seq_len(n_loci))
  gt
}

#' Mendelian segregation: for n_rec new recruits, given parental genotype
#' rows, draw one allele from each parent at each locus and sum.
#' Homozygotes deterministically pass their fixed allele; heterozygotes
#' pass either allele with probability 0.5 (vectorised across recruits).
#' Ported from PlantPopGenFit::assign_alleles_to_recruits(), generalised
#' from a fixed 2-genotype loop to an arbitrary n_loci loop.
inherit_alleles <- function(resist_gt, mother_idx, father_idx) {
  n_rec  <- length(mother_idx)
  n_loci <- ncol(resist_gt)
  rec_gt <- matrix(0L, nrow = n_rec, ncol = n_loci)

  for (j in seq_len(n_loci)) {
    mom <- resist_gt[mother_idx, j]
    dad <- resist_gt[father_idx, j]

    allele_m <- ifelse(mom == 2L, 1L,
                 ifelse(mom == 0L, 0L, rbinom(n_rec, 1, 0.5)))
    allele_f <- ifelse(dad == 2L, 1L,
                 ifelse(dad == 0L, 0L, rbinom(n_rec, 1, 0.5)))

    rec_gt[, j] <- allele_m + allele_f
  }

  colnames(rec_gt) <- colnames(resist_gt)
  rec_gt
}

#' Convert genotype matrix -> a single continuous resistance score in
#' [0, 1], interpreted as "fraction of the rust-specific hazard
#' increment removed". Per-locus contributions are additive across loci
#' (extend here if epistasis between resistance loci is ever needed --
#' see PlantPopGenFit::assign_phenotype_quantitative_epistatic() for a
#' template), then the total is clamped at 1.
#'
#' This single clamp is what lets the same function represent:
#'   - one dominant locus fully protective in RR and Rr      (score = 1)
#'   - one recessive locus protective only in RR             (score = 1 in RR, 0 in Rr)
#'   - oligogenic additive partial resistance, no single locus sufficient
resist_score_from_gt <- function(resist_gt, params) {
  n_loci  <- ncol(resist_gt)
  contrib <- matrix(0, nrow = nrow(resist_gt), ncol = n_loci)

  for (j in seq_len(n_loci)) {
    g <- resist_gt[, j]
    contrib[, j] <- ifelse(g == 2L, params$resist_locus_effect[j],
                     ifelse(g == 1L, params$resist_locus_effect[j] *
                                       params$resist_dominance[j], 0))
  }

  pmin(rowSums(contrib), 1)
}

#' Per-locus allele frequency among currently alive individuals.
#' Used by census_year() to report selection trajectories.
allele_freqs <- function(resist_gt, alive) {
  if (sum(alive) == 0) {
    return(rep(NA_real_, ncol(resist_gt)))
  }
  colMeans(resist_gt[alive, , drop = FALSE]) / 2
}
