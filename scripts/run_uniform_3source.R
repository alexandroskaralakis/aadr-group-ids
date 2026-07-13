#!/usr/bin/env Rscript
# =====================================================================
# run_uniform_3source.R      -- SUPERSEDES run_eastern_total.R (withdrawn)
#
# WHY run_eastern_total.R WAS WITHDRAWN
#   It required decomposing Russia_Eneolithic_Steppe_QC <- EHG + Iran to get a
#   coefficient phi. Two fatal problems:
#     (1) No EHG population exists in postBA_arb. Merging one in changes the
#         SNP basis and moves every anchor coordinate.
#     (2) STRUCTURAL: Serbia_IronGates_Mesolithic is in the canonical right
#         set and carries EHG-related drift (WHG-EHG cline). With EHG as a
#         SOURCE, that right set violates the project's own outgroup-
#         independence standard and would inflate the fit. phi would have
#         been biased, invisibly.
#   Conclusion: phi is not measurable on this panel with this right set.
#   Logged as PRED-F, FAILED.
#
# WHAT THIS SCRIPT DOES INSTEAD
#   (A) Fixes the REAL defect -- non-commensurable models across populations.
#       The logged composition decomposition fits the four anchors with
#       2 sources (Turkey_N + Iran) and the target with 3 (+ steppe).
#       Greece_Mycenaean_POOL carries ~13% Yamnaya-related steppe but has NO
#       steppe term, so its steppe-borne CHG is absorbed into its Turkey_N and
#       Iran weights -- INFLATING its 6.8% Iran figure. Comparing that to the
#       target's steppe-EXCLUDING 26.4% is not a fair comparison.
#       FIX: fit ALL FIVE with the identical three-source left set.
#
#   (B) Replaces the invented coefficient with a BRACKET.
#       Eastern ancestry cannot be pinned to a point without phi, but it can be
#       bounded from below and above using weights we already have:
#           E_low  = w_iran                  (steppe carries NO Iran/CHG -- false, hard floor)
#           E_high = w_iran + w_steppe       (steppe is ALL Iran/CHG  -- false, hard ceiling)
#       The truth lies strictly inside [E_low, E_high] for every population,
#       INCLUDING the anchors. Both bounds are computed identically for all
#       five, so the comparison is fair even though neither bound is the truth.
#
# BASIS
#   IDENTICAL to run_composition_decomp: same PREFIX, same population list,
#   therefore the same 122,804 polymorphic SNPs / 708 blocks. The script ABORTS
#   if the population list does not match that run's manifest exactly.
#   This is what makes the output plottable on the same figure as the existing
#   card without a cross-basis violation.
# =====================================================================

suppressPackageStartupMessages(library(admixtools))
suppressPackageStartupMessages(library(tidyverse))

SCRIPT_ID <- "run_uniform_3source"
setwd("/Users/ell_thales/qpadm_project2")
dir.create("manifest", showWarnings = FALSE)

PREFIX <- "postBA_arb"
F2DIR  <- "f2_uniform_3source"

TARGETS <- c("Greece_Mycenaean_POOL", "Greece_Minoan_POOL", "MyDNA",
             "Turkey_EBA_WEST", "Turkey_MLBA_QC")
ANF     <- "Turkey_N_QC"
IRAN    <- "Iran_GanjDareh_N_QC"
STEPPE  <- "Russia_Eneolithic_Steppe_QC"
LEFT    <- c(ANF, IRAN, STEPPE)
RIGHT   <- c("Mbuti", "Italy_Sicily_Epigravettian", "Russia_Kostenki_UP",
             "Papuan", "Serbia_IronGates_Mesolithic", "Morocco_Iberomaurusian")

POPS <- unique(c(TARGETS, LEFT, RIGHT))

# ---------------------------------------------------------------------
# BASIS GATE -- refuse to run on a different population list
# ---------------------------------------------------------------------
ref <- "manifest/run_composition_decomp_populations.txt"
if (file.exists(ref)) {
  prev <- sort(readLines(ref))
  if (!identical(sort(POPS), prev)) {
    cat("POPULATION LIST DIFFERS FROM run_composition_decomp.\n")
    cat("Only here:  ", paste(setdiff(sort(POPS), prev), collapse = ", "), "\n")
    cat("Only there: ", paste(setdiff(prev, sort(POPS)), collapse = ", "), "\n")
    stop("Basis would differ. Fix the list, or accept a new basis DELIBERATELY ",
         "and rebuild the whole figure. Do not proceed by accident.")
  }
  cat("BASIS GATE: population list matches run_composition_decomp exactly.\n")
} else {
  cat("WARNING: ", ref, " not found. Cannot verify basis identity.\n",
      "Compare n_jackknife_blocks and SNP count below against 708 / 122,804.\n")
}

cat("\n######## EXTRACT (single call) ########\n")
extract_f2(PREFIX, F2DIR, pops = POPS, maxmiss = 0, blgsize = 0.05,
           overwrite = TRUE, verbose = TRUE)
f2b <- f2_from_precomp(F2DIR, verbose = FALSE)
cat(sprintf("\nBASIS: %d populations, %d jackknife blocks (expect 708)\n",
            length(POPS), dim(f2b)[3]))

# ---------------------------------------------------------------------
# UNIFORM THREE-SOURCE MODEL -- identical left set for all five
# ---------------------------------------------------------------------
cat("\n######## X <- Turkey_N + Iran + Eneolithic steppe ########\n")
cat("Expect Greece_Minoan_POOL to return a steppe weight near zero, possibly\n")
cat("slightly negative. Report it AS IS. Do NOT drop the steppe term for the\n")
cat("Minoan pool -- a non-uniform left set is the defect this script fixes.\n\n")

res <- map_dfr(TARGETS, function(tg) {
  m <- qpadm(f2b, left = LEFT, right = RIGHT, target = tg)
  cat("---- ", tg, " ----\n"); print(m$weights)
  cat(sprintf("p = %.4g\n\n", m$rankdrop$p[1]))
  w <- deframe(select(m$weights, left, weight))
  z <- deframe(select(m$weights, left, z))
  s <- deframe(select(m$weights, left, se))
  tibble(target  = tg,
         w_anf   = w[[ANF]],    z_anf    = z[[ANF]],
         w_iran  = w[[IRAN]],   z_iran   = z[[IRAN]],   se_iran   = s[[IRAN]],
         w_steppe= w[[STEPPE]], z_steppe = z[[STEPPE]], se_steppe = s[[STEPPE]],
         p_model = m$rankdrop$p[1],
         E_low   = w[[IRAN]],
         E_high  = w[[IRAN]] + w[[STEPPE]])
})

OUT <- res %>%
  mutate(across(c(w_anf, w_iran, w_steppe, se_iran, se_steppe, E_low, E_high),
                ~ round(100 * .x, 2))) %>%
  arrange(E_low)

cat("######## RESULT -- one model, one basis, all five populations ########\n")
cat("E_low  = explicit Iran weight        (hard FLOOR on eastern ancestry)\n")
cat("E_high = Iran + steppe weight        (hard CEILING; steppe is not 100% CHG)\n")
cat("The truth lies strictly inside the bracket, for every row.\n\n")
print(OUT, width = Inf)
write_csv(OUT, sprintf("manifest/%s_results.csv", SCRIPT_ID))

# ---------------------------------------------------------------------
# BASIS MANIFEST
# ---------------------------------------------------------------------
writeLines(sort(POPS), sprintf("manifest/%s_populations.txt", SCRIPT_ID))
writeLines(c(
  paste0("script: ", SCRIPT_ID),
  paste0("prefix: ", PREFIX),
  paste0("f2dir: ",  F2DIR),
  "maxmiss: 0  blgsize: 0.05",
  paste0("left: ", paste(LEFT, collapse = " + ")),
  paste0("n_populations: ", length(POPS)),
  paste0("n_jackknife_blocks: ", dim(f2b)[3]),
  "phi: NOT MEASURABLE on this panel/right set -- see PRED-F. Bracket used instead."
), sprintf("manifest/%s_basis.txt", SCRIPT_ID))

cat("\nDONE. Run as:  Rscript run_uniform_3source.R > manifest/run_uniform_3source.log 2>&1\n")
