# ============================================================================
# run_chisq_exhibit.R
#
# Paper B §3, rebuilt on public data, with the null it was missing.
#
# THE CLAIM (§3.2): adding a population to `extract_f2` that the MODEL NEVER
# USES collapses the SNP basis and reverses the model's conclusion.
#
# THE ARGUMENT §3 GAVE (§3.3), AND WHY IT IS WRONG:
#     "chi^2 scales approximately with SNP count. A 17.7% loss therefore
#      predicts chi^2 ~ 12.5. Observed: 7.36 -- a 51% fall."
#
#   run_chisq_null.R measured what chi^2 ACTUALLY does under RANDOM SNP loss,
#   on this exact model and panel:
#
#       random loss   "chi^2 ~ N" predicts   OBSERVED chi^2 ratio
#          17.7%             0.82             1.36  (1.15 - 1.52, 20 reps)
#          30%               0.70             1.40
#          50%               0.50             1.14
#
#   chi^2 does not fall with N. IT RISES. Fewer SNPs means fewer jackknife
#   blocks, a noisier f4 covariance matrix, and an inflated inverse -- so
#   chi^2 = E' Q^-1 E is biased UPWARD. This is a small-sample property of
#   the statistic, not of population history.
#
# THEREFORE THE EXHIBIT IS MORE ANOMALOUS THAN §3 CLAIMED, NOT LESS.
#   §3 compared its observed chi^2 against a prediction of ~12.5 (proportional
#   decay). The CORRECT null is ~1.36x the baseline -- an INCREASE. A chi^2
#   that FALLS when the basis shrinks is doing something random loss cannot
#   produce in 20 out of 20 replicates.
#
#   The diagnostic is not "chi^2 fell faster than N."
#   The diagnostic is "chi^2 fell when it should have risen."
#
# THIS SCRIPT tests that on data anyone can download: a thin population is
# added to the extract, the basis collapses, and the resulting chi^2 is
# compared against the RANDOM-LOSS NULL AT MATCHED N.
#
# Run: Rscript run_chisq_exhibit.R 2>&1 | tee manifest/run_chisq_exhibit.log
# ============================================================================

suppressMessages({library(admixtools); library(dplyr)})

OUT   <- "clean_thin"
PRE   <- file.path(OUT, "thin")
STAMP <- format(Sys.time(), "%Y%m%d_%H%M%S")
set.seed(20260714)

right  <- readLines(file.path(OUT, "right.txt")); right <- right[nzchar(right)]
THIN   <- readLines(file.path(OUT, "thin.txt"));  THIN  <- THIN[nzchar(THIN)]
KOTIAS <- "Georgia_KotiasKlde_Mesolithic"
TARGET <- "Germany_Esperstedt_CordedWare"
LEFT   <- c("Russia_Samara_EBA_Yamnaya", "Czechia_N_GlobularAmphora")
RT     <- c(right, KOTIAS)          # the rejected model

NREP <- 20

cat("==========================================================\n")
cat(" run_chisq_exhibit.R  ", format(Sys.time()), "\n")
cat("==========================================================\n\n")
cat("Model (FIXED throughout):", TARGET, "<-", paste(LEFT, collapse=" + "), "\n")
cat("Right (FIXED throughout): six canonical + Kotias\n")
cat("The ONLY thing that changes is which populations are in the EXTRACT.\n\n")

cat("PRE-REGISTERED:\n")
cat("  P1. Each thin population collapses the SNP basis by >10%.        [~85%]\n")
cat("  P2. chi^2 FALLS when a thin population enters the extract, while\n")
cat("      random loss of the same magnitude makes it RISE.             [~70%]\n")
cat("  P3. The observed chi^2 lies BELOW the entire random-loss null at\n")
cat("      matched N -- i.e. outside 20/20 replicates.                  [~60%]\n")
cat("  P4. At least one thin population flips the model from REJECTED to\n")
cat("      NOT REJECTED -- a conclusion reversed by a population the model\n")
cat("      never uses.                                                  [~70%]\n")
cat("  FAILURE MODE: if chi^2 rises (as under random loss), the exhibit is\n")
cat("      just power loss and §3 must be retracted.\n\n")

fit <- function(pops, keep, tag) {
  d <- file.path(OUT, sprintf("f2_%s_%s", tag, STAMP))
  if (dir.exists(d)) stop("refusing to reuse: ", d)
  dir.create(d, recursive = TRUE)
  args <- list(PRE, outdir = d, pops = pops, maxmiss = 0,
               auto_only = TRUE, verbose = FALSE)
  if (!is.null(keep)) args$keepsnps <- keep
  do.call(extract_f2, args)
  f2 <- f2_from_precomp(d, verbose = FALSE)
  n  <- tryCatch(sum(count_snps(f2)), error = function(e) NA_integer_)
  m  <- qpadm(f2, left = LEFT, right = RT, target = TARGET)
  unlink(d, recursive = TRUE)
  list(n = n, chisq = m$rankdrop$chisq[1], p = m$rankdrop$p[1],
       dof = m$rankdrop$dof[1],
       w = paste(sprintf("%.1f%%", 100*m$weights$weight), collapse = " / "))
}

CORE <- c(TARGET, LEFT, KOTIAS, right)
snps <- read.table(paste0(PRE, ".snp"), stringsAsFactors = FALSE)$V1

# ------------------------------------------------------------- the baseline --
cat("##### 1. baseline: the core extract #####\n")
b <- fit(CORE, NULL, "core")
cat(sprintf("  N = %s   chi^2 = %.2f   dof %d   p = %.4g   [%s]\n\n",
            format(b$n, big.mark=","), b$chisq, b$dof, b$p,
            if (b$p < 0.05) "REJECTS" else "does not reject"))
if (b$p >= 0.05) cat("  *** WARNING: baseline does not reject. chi^2 scaling arguments\n",
                     "  *** do not apply to an accepted model. Read the null run.\n\n")

# ---------------------------------------------- 2. the exhibit: thin pops ----
cat("##### 2. THE EXHIBIT: add a thin population to the EXTRACT ONLY #####\n")
cat("The model does not change. The right set does not change. Only the extract.\n\n")

ex <- bind_rows(lapply(THIN, function(T) {
  r <- fit(c(CORE, T), NULL, paste0("thin_", gsub("[^A-Za-z]", "", T)))
  loss <- 1 - r$n / b$n
  cat(sprintf("  + %-26s N = %-9s (%+.1f%%)   chi^2 = %-6.2f  p = %-8.4g  %s\n",
              T, format(r$n, big.mark=","), -100*loss, r$chisq, r$p,
              if (r$p < 0.05) "still rejects" else ">>> NO LONGER REJECTS <<<"))
  tibble(added = T, n = r$n, loss = loss, chisq = r$chisq, p = r$p,
         chisq_ratio = r$chisq / b$chisq, weights = r$w)
}))
cat("\n")

# ------------------------------- 3. the null: RANDOM loss at matched N -------
cat("##### 3. THE NULL: random SNP loss to the SAME basis #####\n")
cat("Uniform loss by construction. This is what chi^2 does when the SNPs that\n")
cat("go missing are NOT ascertained.\n\n")

nulls <- bind_rows(lapply(seq_len(nrow(ex)), function(i) {
  E <- ex[i, ]
  target_n <- E$n
  # match the basis by random subsampling of the CORE extract
  frac <- target_n / b$n
  reps <- bind_rows(lapply(seq_len(NREP), function(j) {
    keep <- sample(snps, round(length(snps) * frac))
    r <- fit(CORE, keep, sprintf("null%d_%d", i, j))
    tibble(added = E$added, rep = j, n = r$n, chisq = r$chisq, p = r$p)
  }))
  lo <- min(reps$chisq); hi <- max(reps$chisq); md <- median(reps$chisq)
  outside <- E$chisq < lo
  cat(sprintf("  matched to '+ %s'  (N ~ %s, %.1f%% loss)\n",
              E$added, format(target_n, big.mark=","), 100*E$loss))
  cat(sprintf("    random-loss chi^2:  median %.2f   range %.2f - %.2f   (%d reps)\n",
              md, lo, hi, NREP))
  cat(sprintf("    \"chi^2 ~ N\" would predict: %.2f\n", b$chisq * frac))
  cat(sprintf("    OBSERVED with the thin pop: %.2f   %s\n\n", E$chisq,
              if (outside) ">>> BELOW THE ENTIRE NULL <<<" else "(inside the null range)"))
  reps %>% mutate(exhibit_chisq = E$chisq, null_min = lo, null_max = hi,
                  null_median = md, outside_null = outside,
                  naive_prediction = b$chisq * frac)
}))

# ---------------------------------------------------------------- verdict ----
cat("==========================================================\n")
cat(" VERDICT\n")
cat("==========================================================\n\n")

summ <- nulls %>% group_by(added) %>%
  summarise(exhibit_chisq = first(exhibit_chisq),
            null_min = first(null_min), null_median = first(null_median),
            null_max = first(null_max), outside = first(outside_null),
            naive = first(naive_prediction), .groups = "drop") %>%
  left_join(ex %>% select(added, n, loss, p), by = "added")
print(as.data.frame(summ), row.names = FALSE, digits = 4)

nout <- sum(summ$outside)
flips <- sum(ex$p > 0.05)

cat("\n")
cat(sprintf("  Thin populations that collapse the basis:  %d of %d\n",
            sum(ex$loss > 0.10), nrow(ex)))
cat(sprintf("  Conclusions REVERSED (reject -> not reject): %d of %d\n", flips, nrow(ex)))
cat(sprintf("  chi^2 BELOW the entire random-loss null:     %d of %d\n\n", nout, nrow(summ)))

if (nout > 0) {
  cat("  §3 SURVIVES, ON PUBLIC DATA, WITH THE CORRECT NULL.\n\n")
  cat("  Under random SNP loss chi^2 RISES. With a thin population in the\n")
  cat("  extract it FALLS -- below every one of 20 random-loss replicates at\n")
  cat("  the same basis. The surviving SNPs are that genome's coverage\n")
  cat("  footprint: an ascertained subset, not a random one.\n\n")
  cat("  RESTATE THE DIAGNOSTIC:\n")
  cat("    NOT  'chi^2 fell faster than the SNP count'   (the premise is false)\n")
  cat("    BUT  'chi^2 fell when random loss says it should have risen'\n")
} else {
  cat("  §3 DOES NOT SURVIVE THE REBUILD.\n")
  cat("  The observed chi^2 lies inside the random-loss null. The exhibit does\n")
  cat("  not distinguish ascertained loss from random loss. RETRACT §3's\n")
  cat("  signature; keep only the mechanical point, which is true and is about\n")
  cat("  the SNP BASIS: a population absent from the model still sets it.\n")
}
cat("\n----------------------------------------------------------\n\n")

dir.create("manifest", showWarnings = FALSE)
write.csv(ex,    "manifest/run_chisq_exhibit.csv", row.names = FALSE)
write.csv(nulls, "manifest/run_chisq_exhibit_null.csv", row.names = FALSE)
writeLines(c(
  "script: run_chisq_exhibit",
  paste("date:", format(Sys.time())),
  "purpose: Paper B section 3's exhibit, rebuilt on public data, with a null",
  paste0("model: ", TARGET, " <- ", paste(LEFT, collapse=" + ")),
  "right: six canonical + Georgia_KotiasKlde_Mesolithic",
  paste0("baseline: N=", b$n, " chisq=", round(b$chisq,3), " dof=", b$dof,
         " p=", signif(b$p,4)),
  "",
  "thin populations added to the EXTRACT ONLY (never in the model):",
  paste0("  ", ex$added, ": N=", ex$n, " (", sprintf("%.1f%%", -100*ex$loss),
         ") chisq=", round(ex$chisq,3), " p=", signif(ex$p,4)),
  "",
  paste0("random-loss null: ", NREP, " replicates per matched basis"),
  paste0("chi^2 below the entire null: ", nout, " of ", nrow(summ)),
  "",
  "THE PREMISE OF SECTION 3.3 IS FALSE. Under random SNP loss chi^2 does not",
  "  fall with N; it RISES (1.36x at 17.7% loss, 20/20 reps; run_chisq_null).",
  "  Fewer SNPs -> fewer jackknife blocks -> noisier f4 covariance -> inflated",
  "  inverse -> chi^2 = E'Q^-1 E biased upward. The exhibit is therefore MORE",
  "  anomalous than section 3 claimed, and the diagnostic must be restated:",
  "  chi^2 fell when it should have RISEN.",
  "",
  "Every population is an AADR group taken whole. No target genome, no pool.",
  "NOTHING here is comparable to 124,362 / 102,314 or any earlier basis."
), "manifest/run_chisq_exhibit.txt")

cat("Manifest: manifest/run_chisq_exhibit.txt\nDone.\n")
