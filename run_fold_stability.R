# ============================================================================
# run_fold_stability.R
#
# Paper B v1.0 reported its Test 3 fold across THREE SNP bases (43.8x / 41.9x /
# 41.9x) and argued that the FOLD, not the p-value, is the reportable quantity
# BECAUSE the fold was stable while the p-value was not.
#
# The rebuilt Test 3 (run_test3_controls.R) has been measured on ONE basis:
#     Germany_Esperstedt_CordedWare <- Yamnaya + Czechia_N_GlobularAmphora
#     six refs p = 0.267  ->  six + Kotias p = 0.016   FOLD = 16.6x
#     (898,057 SNPs)
#
# One basis does not support the claim that the fold is the stable quantity.
# This script tests it.
#
# DESIGN
#   Same model, same populations, same right set. ONLY the SNP basis changes:
#     - full basis
#     - random subsamples at ~500k, ~250k, ~120k, THREE replicates each
#     - transversions only (an independent kind of thinning: it removes
#       damage-prone transitions rather than sampling at random)
#
#   ~120k is the point of the exercise: it is the basis a typical 1240K aDNA
#   study actually works at, and it is where Paper B v1.0's numbers lived.
#   If Test 3 only fires at ~900k SNPs, it is not a tool most studies can use,
#   and the paper MUST say so.
#
#   At every basis we also re-run the six legitimate references (leave-one-out
#   and restore). If discrimination collapses at low basis, that is the result.
#
# Run: Rscript run_fold_stability.R 2>&1 | tee manifest/run_fold_stability.log
# ============================================================================

suppressMessages({library(admixtools); library(dplyr); library(tidyr)})

OUT    <- "clean_controls"
PRE    <- file.path(OUT, "controls")
STAMP  <- format(Sys.time(), "%Y%m%d_%H%M%S")
set.seed(20260713)

right  <- readLines(file.path(OUT, "right.txt")); right <- right[nzchar(right)]
KOTIAS <- "Georgia_KotiasKlde_Mesolithic"
TARGET <- "Germany_Esperstedt_CordedWare"
LEFT   <- c("Russia_Samara_EBA_Yamnaya", "Czechia_N_GlobularAmphora")

cat("==========================================================\n")
cat(" run_fold_stability.R  ", format(Sys.time()), "\n")
cat("==========================================================\n\n")
cat("Model (fixed):", TARGET, "<-", paste(LEFT, collapse = " + "), "\n")
cat("Candidate reference:", KOTIAS, "\n")
cat("ONLY the SNP basis changes.\n\n")

cat("PRE-REGISTERED:\n")
cat("  P1. The fold stays > 5x at every basis.                        [~70%]\n")
cat("  P2. The fold varies < 3x across bases while p-values move by\n")
cat("      orders of magnitude. (This is the argument for reporting\n")
cat("      the fold rather than the p-value.)                         [~60%]\n")
cat("  P3. At ~120k SNPs, Kotias may no longer BREAK the model (p stays\n")
cat("      above 0.05) even if the fold persists. If so, Test 3 needs a\n")
cat("      high-basis regime that most published aDNA studies do not\n")
cat("      have, and that is a MAJOR limitation.                      [~50%]\n")
cat("  P4. The six legitimate references keep passing at every basis.  [~75%]\n\n")

adm <- function(m) {
  p <- m$rankdrop$p[1]
  isTRUE(p > 0.05) && all(m$weights$weight > 0) && all(abs(m$weights$z) > 2)
}

snpinfo <- read.table(paste0(PRE, ".snp"), stringsAsFactors = FALSE)
snps    <- snpinfo$V1
# transversions: exclude A<->G and C<->T
ref <- toupper(snpinfo$V5); alt <- toupper(snpinfo$V6)
is_ts <- (ref == "A" & alt == "G") | (ref == "G" & alt == "A") |
         (ref == "C" & alt == "T") | (ref == "T" & alt == "C")
tv <- snps[!is_ts]
cat("panel SNPs:", format(length(snps), big.mark = ","),
    "  transversions:", format(length(tv), big.mark = ","), "\n\n")

runs <- c(
  list(list(tag = "full", keep = NULL, lab = "full basis")),
  list(list(tag = "tv",   keep = tv,   lab = "transversions only")),
  unlist(lapply(c(500000, 250000, 120000), function(n)
    lapply(LETTERS[1:3], function(r)
      list(tag = paste0(n/1000, "k_", r), keep = n,
           lab = sprintf("~%dk subsample (rep %s)", n/1000, r)))),
    recursive = FALSE)
)

one_basis <- function(r) {
  d <- file.path(OUT, sprintf("f2_fold_%s_%s", r$tag, STAMP))
  if (dir.exists(d)) stop("refusing to reuse: ", d)
  dir.create(d, recursive = TRUE)

  keep <- if (is.null(r$keep)) NULL
          else if (is.numeric(r$keep)) sample(snps, min(r$keep, length(snps)))
          else r$keep

  args <- list(PRE, outdir = d, pops = c(TARGET, LEFT, KOTIAS, right),
               maxmiss = 0, auto_only = TRUE, verbose = FALSE)
  if (!is.null(keep)) args$keepsnps <- keep
  do.call(extract_f2, args)
  f2 <- f2_from_precomp(d, verbose = FALSE)
  n  <- tryCatch(sum(count_snps(f2)), error = function(e) NA_integer_)

  m6 <- qpadm(f2, left = LEFT, right = right,            target = TARGET)
  m7 <- qpadm(f2, left = LEFT, right = c(right, KOTIAS), target = TARGET)
  p6 <- m6$rankdrop$p[1]; p7 <- m7$rankdrop$p[1]

  # do the six legitimate references still pass at this basis?
  sixfold <- sapply(right, function(O) {
    a <- qpadm(f2, left = LEFT, right = setdiff(right, O), target = TARGET)
    a$rankdrop$p[1] / p6
  })

  cat(sprintf("  %-26s basis %9s   p6 = %-9.4g  p7 = %-9.4g  FOLD %7.1fx   %s\n",
              r$lab, format(n, big.mark = ","), p6, p7, p6/p7,
              if (adm(m6)) {
                if (!adm(m7)) "control HOLDS, Kotias BREAKS it"
                else          "control holds, Kotias does NOT break it"
              } else "CONTROL DOES NOT HOLD"))
  cat(sprintf("      six legit refs: max fold %.2fx  (%s)\n",
              max(sixfold),
              if (max(sixfold) < 5) "all pass — discriminates"
              else "*** A LEGIT REF ALSO FIRES — DISCRIMINATION LOST ***"))

  tibble(basis_label = r$lab, basis = n,
         p_six = p6, p_kotias = p7, fold = p6/p7,
         control_holds = adm(m6), kotias_breaks = !adm(m7),
         six_max_fold = max(sixfold), discriminates = max(sixfold) < 5)
}

cat("##### the fold across bases #####\n\n")
R <- bind_rows(lapply(runs, one_basis))

cat("\n==========================================================\n")
cat(" SUMMARY\n")
cat("==========================================================\n\n")
print(as.data.frame(R), row.names = FALSE, digits = 4)

# ---- the three questions the paper needs answered -------------------------
held <- R %>% filter(control_holds)
cat("\n----------------------------------------------------------\n")

cat("Q1. Does the fold stay above 5x?\n")
cat(sprintf("    %d of %d bases have fold > 5x. Range: %.1fx - %.1fx.\n",
            sum(R$fold > 5), nrow(R), min(R$fold), max(R$fold)))

cat("\nQ2. Is the FOLD more stable than the P-VALUE?\n")
cat(sprintf("    fold    spans %.1fx  (%.1f - %.1f)\n",
            max(R$fold)/min(R$fold), min(R$fold), max(R$fold)))
cat(sprintf("    p_six   spans %.0f-fold  (%.3g - %.3g)\n",
            max(R$p_six)/min(R$p_six), min(R$p_six), max(R$p_six)))
cat(sprintf("    p_kotias spans %.0f-fold (%.3g - %.3g)\n",
            max(R$p_kotias)/min(R$p_kotias), min(R$p_kotias), max(R$p_kotias)))
cat("    >> If the fold spans much less than the p-values, REPORT THE FOLD.\n")

cat("\nQ3. Does Test 3 still FIRE at aDNA-typical bases (~120k)?\n")
low <- R %>% filter(grepl("120k", basis_label))
print(as.data.frame(low %>% select(basis, p_six, p_kotias, fold,
                                   control_holds, kotias_breaks)),
      row.names = FALSE, digits = 4)
if (any(low$control_holds & !low$kotias_breaks)) {
  cat("\n    *** At ~120k SNPs the control HOLDS but Kotias does NOT break it.\n")
  cat("    *** Test 3 requires a high SNP basis. Most published 1240K aDNA\n")
  cat("    *** studies work well below 900k. THIS IS A MAJOR LIMITATION AND\n")
  cat("    *** MUST BE STATED IN THE ABSTRACT, NOT THE DISCUSSION.\n")
} else if (all(low$control_holds & low$kotias_breaks)) {
  cat("\n    Test 3 fires at ~120k. It is usable at aDNA-typical bases.\n")
} else {
  cat("\n    The control itself does not hold at ~120k. Read the rows above.\n")
}

if (!all(R$discriminates)) {
  cat("\n    *** DISCRIMINATION LOST at one or more bases: a LEGITIMATE\n")
  cat("    *** reference produces a Kotias-sized fold. Test 3 would then be\n")
  cat("    *** reacting to perturbation, not to the violation. Report it.\n")
}
cat("----------------------------------------------------------\n\n")

dir.create("manifest", showWarnings = FALSE)
write.csv(R, "manifest/run_fold_stability.csv", row.names = FALSE)
writeLines(c(
  "script: run_fold_stability",
  paste("date:", format(Sys.time())),
  "purpose: is the Test 3 fold a property of the data, or of the SNP basis?",
  paste0("model: ", TARGET, " <- ", paste(LEFT, collapse = " + ")),
  paste0("candidate reference: ", KOTIAS),
  "",
  paste0("fold range: ", sprintf("%.1fx - %.1fx", min(R$fold), max(R$fold))),
  paste0("p_six range: ", sprintf("%.3g - %.3g", min(R$p_six), max(R$p_six))),
  paste0("bases tested: ", paste(format(R$basis, big.mark=","), collapse=", ")),
  "",
  "Subsamples are random draws from the panel SNP set, three replicates per",
  "  size, seed 20260713. The transversion-only row is an independent kind of",
  "  thinning (it removes damage-prone transitions, not a random subset).",
  "",
  "Every row is a SEPARATE BASIS. The rows are compared to each other ONLY",
  "  through the fold, which is the quantity this script exists to test.",
  "NOTHING here is comparable to 122,891 / 124,362 or any earlier basis."
), "manifest/run_fold_stability.txt")

cat("Manifest: manifest/run_fold_stability.txt\nDone.\n")
