# ============================================================================
# run_chisq_null.R
#
# Paper B §3 claims a diagnostic:
#
#     "chi^2 scales approximately with SNP count. A 17.7% loss therefore
#      predicts chi^2 ~ 12.5. Observed: 7.36 -- a 51% fall."
#     ">> chi^2 falling faster than the SNP count is a diagnostic. It
#      distinguishes 'there is no signal' from 'the signal-bearing sites
#      were dropped.'"
#
# THE INFERENCE HAS NO NULL DISTRIBUTION. Section 3.4 admits this:
#     [OWED: a simulation calibrating the expected chi^2/N ratio under
#      uniform vs footprint-driven SNP loss. Without it this is an
#      observation, not a test with a threshold.]
#
# This script builds that null.
#
# METHOD
#   Take a REJECTED model on the clean, public controls panel. Remove SNPs
#   AT RANDOM -- uniform loss BY CONSTRUCTION, which is precisely the null
#   the diagnostic needs. Measure what chi^2 actually does.
#
#   If a 51% fall in chi^2 is COMMON under random 17.7% SNP loss, then
#   Paper B's exhibit is within noise and the "signature" is not a signature.
#
# THE MODEL
#   Germany_Esperstedt_CordedWare <- Yamnaya + Czechia_N_GlobularAmphora,
#   right = six canonical + Georgia_Kotias.
#   At 898,057 SNPs this REJECTS (p = 0.016), which is what we need: chi^2
#   only scales with N for a model with genuine misfit. It is also a close
#   analogue of the exhibit's models (moderate rejection, similar dof).
#
# A POINT PAPER B NEVER MAKES, AND MUST
#   chi^2 ~ N holds only for a FALSE model. For a TRUE model the misfit
#   vanishes as ~1/N, so chi^2 converges to its null distribution and is
#   BASIS-INVARIANT. (This is why the six-reference control's p-value barely
#   moved from 898k to 108k SNPs.) Section 3 states "chi^2 scales with SNP
#   count" unconditionally. It does not. The script checks both cases.
#
# Run: Rscript run_chisq_null.R 2>&1 | tee manifest/run_chisq_null.log
# ============================================================================

suppressMessages({library(admixtools); library(dplyr)})

OUT   <- "clean_controls"
PRE   <- file.path(OUT, "controls")
STAMP <- format(Sys.time(), "%Y%m%d_%H%M%S")
set.seed(20260714)

right  <- readLines(file.path(OUT, "right.txt")); right <- right[nzchar(right)]
KOTIAS <- "Georgia_KotiasKlde_Mesolithic"
TARGET <- "Germany_Esperstedt_CordedWare"
LEFT   <- c("Russia_Samara_EBA_Yamnaya", "Czechia_N_GlobularAmphora")

RIGHT_REJ <- c(right, KOTIAS)   # the REJECTED model  (p = 0.016 at full basis)
RIGHT_OK  <- right              # the ACCEPTED model  (p = 0.267 at full basis)

NREP <- 20          # replicates per loss fraction
LOSS <- c(0.177, 0.30, 0.50)   # 0.177 = the exhibit's loss fraction

cat("==========================================================\n")
cat(" run_chisq_null.R  ", format(Sys.time()), "\n")
cat("==========================================================\n\n")
cat("Model:", TARGET, "<-", paste(LEFT, collapse = " + "), "\n")
cat("SNP loss is RANDOM by construction. That is the null Section 3 needs.\n\n")

cat("PRE-REGISTERED:\n")
cat("  P1. chi^2 is far noisier than Section 3 assumes. Under RANDOM 17.7%\n")
cat("      loss, chi^2 falls by >= 51% in more than 10% of replicates.  [~65%]\n")
cat("  P2. If P1 holds, the Section 3 exhibit is within noise and the\n")
cat("      'signature' has no discriminating power.       [conditional on P1]\n")
cat("  P3. chi^2 ~ N holds only for a REJECTED model. For an ACCEPTED model\n")
cat("      chi^2 is basis-invariant. Section 3 states the scaling law\n")
cat("      unconditionally, and it is wrong to.                        [~80%]\n")
cat("  P4. If the fall IS anomalous, Section 3 survives and finally gets\n")
cat("      the null distribution it admits it owes.\n\n")

snps <- read.table(paste0(PRE, ".snp"), stringsAsFactors = FALSE)$V1

fit_at <- function(keep, rt, tag) {
  d <- file.path(OUT, sprintf("f2_null_%s_%s", tag, STAMP))
  if (dir.exists(d)) stop("refusing to reuse: ", d)
  dir.create(d, recursive = TRUE)
  args <- list(PRE, outdir = d, pops = c(TARGET, LEFT, KOTIAS, right),
               maxmiss = 0, auto_only = TRUE, verbose = FALSE)
  if (!is.null(keep)) args$keepsnps <- keep
  do.call(extract_f2, args)
  f2 <- f2_from_precomp(d, verbose = FALSE)
  n  <- tryCatch(sum(count_snps(f2)), error = function(e) NA_integer_)
  m  <- qpadm(f2, left = LEFT, right = rt, target = TARGET)
  unlink(d, recursive = TRUE)          # these pile up; the manifest records them
  list(n = n, chisq = m$rankdrop$chisq[1], p = m$rankdrop$p[1],
       dof = m$rankdrop$dof[1])
}

# ---------------------------------------------------------------- baseline --
cat("##### baseline, full basis #####\n")
b_rej <- fit_at(NULL, RIGHT_REJ, "base_rej")
b_ok  <- fit_at(NULL, RIGHT_OK,  "base_ok")
cat(sprintf("  REJECTED model (7 right): N = %s  chi^2 = %.2f  dof %d  p = %.4g\n",
            format(b_rej$n, big.mark=","), b_rej$chisq, b_rej$dof, b_rej$p))
cat(sprintf("  ACCEPTED model (6 right): N = %s  chi^2 = %.2f  dof %d  p = %.4g\n\n",
            format(b_ok$n, big.mark=","), b_ok$chisq, b_ok$dof, b_ok$p))

# ------------------------------------------------- P3: does chi^2 scale? ----
cat("##### P3 -- does chi^2 scale with N? Rejected vs accepted model #####\n")
cat("If Section 3's law is unconditional, BOTH should fall proportionally.\n\n")
p3 <- bind_rows(lapply(c(0.50, 0.75), function(frac) {
  keep <- sample(snps, round(length(snps) * frac))
  r <- fit_at(keep, RIGHT_REJ, sprintf("p3rej_%d", frac*100))
  o <- fit_at(keep, RIGHT_OK,  sprintf("p3ok_%d",  frac*100))
  cat(sprintf("  ~%.0f%% of SNPs (N = %s)\n", frac*100, format(r$n, big.mark=",")))
  cat(sprintf("     REJECTED: chi^2 %.2f -> %.2f  (%.0f%% of baseline; N is %.0f%%)\n",
              b_rej$chisq, r$chisq, 100*r$chisq/b_rej$chisq, 100*r$n/b_rej$n))
  cat(sprintf("     ACCEPTED: chi^2 %.2f -> %.2f  (%.0f%% of baseline)\n",
              b_ok$chisq, o$chisq, 100*o$chisq/b_ok$chisq))
  tibble(frac = frac, n = r$n,
         chisq_rej = r$chisq, chisq_ok = o$chisq,
         ratio_rej = r$chisq/b_rej$chisq, ratio_ok = o$chisq/b_ok$chisq,
         n_ratio = r$n/b_rej$n)
}))
cat("\n")

# ------------------------------------- P1/P2: the null under random loss ----
cat("##### P1/P2 -- the NULL: chi^2 under RANDOM SNP loss #####\n")
cat(sprintf("Rejected model. %d replicates per loss fraction.\n", NREP))
cat("Section 3's exhibit: 17.7%% loss, chi^2 fell 51%% (15.14 -> 7.36).\n\n")

null <- bind_rows(lapply(LOSS, function(L) {
  keep_n <- round(length(snps) * (1 - L))
  reps <- bind_rows(lapply(seq_len(NREP), function(i) {
    r <- fit_at(sample(snps, keep_n), RIGHT_REJ, sprintf("null_%d_%d", L*1000, i))
    tibble(loss = L, rep = i, n = r$n, chisq = r$chisq, p = r$p)
  }))
  reps %>% mutate(chisq_ratio = chisq / b_rej$chisq,
                  n_ratio     = n / b_rej$n)
}))

for (L in LOSS) {
  s <- null %>% filter(loss == L)
  pred <- 1 - L                                  # what "chi^2 ~ N" predicts
  frac51 <- mean(s$chisq_ratio <= 0.49)          # a >=51% fall, as in the exhibit
  cat(sprintf("--- random loss %.1f%%  (predicted chi^2 ratio if chi^2 ~ N: %.2f) ---\n",
              L*100, pred))
  cat(sprintf("    observed chi^2 ratio: median %.2f   range %.2f - %.2f\n",
              median(s$chisq_ratio), min(s$chisq_ratio), max(s$chisq_ratio)))
  cat(sprintf("    observed chi^2:       median %.2f   range %.2f - %.2f\n",
              median(s$chisq), min(s$chisq), max(s$chisq)))
  cat(sprintf("    p-value flipped to non-rejection in %d of %d replicates\n",
              sum(s$p > 0.05), nrow(s)))
  if (abs(L - 0.177) < 1e-9)
    cat(sprintf("    >>> chi^2 fell by >= 51%% (the exhibit's claim) in %.0f%% of replicates <<<\n",
                100*frac51))
  cat("\n")
}

# ------------------------------------------------------------- the verdict --
s177 <- null %>% filter(abs(loss - 0.177) < 1e-9)
frac51 <- mean(s177$chisq_ratio <= 0.49)

cat("==========================================================\n")
cat(" VERDICT\n")
cat("==========================================================\n\n")
cat(sprintf("Under RANDOM 17.7%% SNP loss, chi^2 fell by 51%% or more in\n"))
cat(sprintf("  %.0f%% of %d replicates.\n\n", 100*frac51, nrow(s177)))

if (frac51 > 0.10) {
  cat("  SECTION 3'S SIGNATURE IS NOT A SIGNATURE.\n")
  cat("  A 51% fall in chi^2 after a 17.7% SNP loss is a COMMON outcome of\n")
  cat("  RANDOM loss. The exhibit does not distinguish 'the signal-bearing\n")
  cat("  sites were dropped' from 'chi^2 is a noisy statistic'. The diagnostic\n")
  cat("  has no discriminating power at dof ~4-6, and Section 3 must be\n")
  cat("  RETRACTED or reduced to the mechanical observation that survives:\n")
  cat("  that an extract change reverses a conclusion, which is true, and is\n")
  cat("  about the SNP BASIS, not about chi^2.\n")
} else if (frac51 > 0.05) {
  cat("  MARGINAL. A 51% fall is uncommon but not rare under random loss.\n")
  cat("  Section 3 needs the null reported alongside the exhibit, and cannot\n")
  cat("  claim the fall is diagnostic without a stated false-positive rate.\n")
} else {
  cat("  SECTION 3 SURVIVES. A 51% fall is rare under random loss. The\n")
  cat("  exhibit's chi^2 collapse IS anomalous, and this run is the null\n")
  cat("  distribution Section 3.4 admits it owes. Report both.\n")
}

cat("\n--- and P3, which Section 3 gets wrong either way ---\n")
print(as.data.frame(p3), row.names = FALSE, digits = 3)
cat("\n  If ratio_rej tracks n_ratio but ratio_ok does not, then chi^2 ~ N\n")
cat("  holds ONLY for a rejected model. Section 3 states the law\n")
cat("  unconditionally. That is wrong, and it matters: an analyst applying\n")
cat("  the chi^2/N rule to a model that FITS will conclude non-random SNP\n")
cat("  loss from a statistic that was never going to scale.\n\n")

dir.create("manifest", showWarnings = FALSE)
write.csv(null, "manifest/run_chisq_null.csv", row.names = FALSE)
write.csv(p3,   "manifest/run_chisq_null_scaling.csv", row.names = FALSE)
writeLines(c(
  "script: run_chisq_null",
  paste("date:", format(Sys.time())),
  "purpose: the null distribution of chi^2 under RANDOM SNP loss --",
  "  the calibration Paper B section 3.4 admits it owes",
  paste0("model: ", TARGET, " <- ", paste(LEFT, collapse=" + ")),
  paste0("baseline (7 right, REJECTED): N=", b_rej$n, " chisq=", round(b_rej$chisq,3),
         " dof=", b_rej$dof, " p=", signif(b_rej$p,4)),
  paste0("baseline (6 right, ACCEPTED): N=", b_ok$n,  " chisq=", round(b_ok$chisq,3),
         " dof=", b_ok$dof,  " p=", signif(b_ok$p,4)),
  paste0("replicates per loss fraction: ", NREP),
  paste0("fraction of replicates with a >=51% chi^2 fall at 17.7% random loss: ",
         sprintf("%.2f", frac51)),
  "",
  "Section 3's exhibit reports a 51% chi^2 fall after a 17.7% SNP loss and",
  "  reads it as evidence of NON-RANDOM loss. This script asks how often",
  "  RANDOM loss does the same thing.",
  "",
  "Every population here is an AADR group taken whole. No target genome.",
  "NOTHING here is comparable to 122,891 / 124,362 or any earlier basis."
), "manifest/run_chisq_null.txt")

cat("Manifest: manifest/run_chisq_null.txt\nDone.\n")
