# ============================================================================
# run_myc_power.R
#
# Myc_PELOPONNESE <- Crete_Lasithi_EMBA + Yamnaya rejects at p = 8.6e-07 on
# 970,660 SNPs. TWO EXPLANATIONS, OPPOSITE CONSEQUENCES:
#
#   (a) POWER.  At ~1M SNPs qpAdm rejects models that are approximately true.
#       If NOTHING fits at this basis, the rejection is uninformative about the
#       descent relationship, and Test 3 survives with a stated basis regime.
#
#   (b) MISFIT. If some other two-source model DOES fit while Lasithi+Yamnaya
#       does not, the descent relationship is wrong and Test 3 has no control.
#
# PART A -- the power curve. Same model, SNPs randomly subsampled to 120k /
#           250k / 500k / full. Two replicates at 120k. If p climbs above 0.05
#           as the basis shrinks toward Paper B's original 122,891, then the
#           original 'p = 0.5113' was a power artifact and so is its fold.
#
# PART B -- the competitor grid. Every two-source model of Myc_PELOPONNESE on
#           ONE full basis. If every one of them rejects, that is (a). If one
#           fits, that is (b).
#
# NOTE ON THE BUG THIS SCRIPT EXISTS TO CORRECT:
#   run_test3_clean.R's adm() checked weights and Z ONLY. It printed
#   [ADMISSIBLE] for a model qpAdm rejects at p = 8.6e-07, and its verdict
#   logic then reported a 111,374x fold between two dead models. Admissibility
#   REQUIRES the model p-value. Fixed here, and logged.
#
# Run: Rscript run_myc_power.R 2>&1 | tee manifest/run_myc_power.log
# ============================================================================

suppressMessages({library(admixtools); library(dplyr); library(tidyr)})

OUT   <- "clean_paperB_wide"
PRE   <- file.path(OUT, "wide")
STAMP <- format(Sys.time(), "%Y%m%d_%H%M%S")
set.seed(20260713)

right  <- readLines(file.path(OUT, "right.txt")); right <- right[nzchar(right)]
TARGET <- "Myc_PELOPONNESE"

cat("==========================================================\n")
cat(" run_myc_power.R  ", format(Sys.time()), "\n")
cat("==========================================================\n\n")

cat("PRE-REGISTERED:\n")
cat("  P1. At ~120k SNPs, Lasithi+Yamnaya does NOT reject (p > 0.05).      [~65%]\n")
cat("  P2. At full basis, EVERY two-source model rejects.                  [~60%]\n")
cat("  P3. If P2 holds, the rejection is a power regime, not a misfit, and\n")
cat("      Test 3 is recoverable with a declared basis.\n")
cat("  P4. If some model fits at full basis and Lasithi+Yamnaya does not,\n")
cat("      the descent relationship is WRONG. Test 3 loses its control and\n")
cat("      Paper B loses its headline. Report it.\n\n")

# --- admissibility, done properly this time --------------------------------
#   A model is admissible ONLY if: p > 0.05 AND all weights > 0 AND all |Z| > 2.
adm <- function(m) {
  p <- m$rankdrop$p[1]
  isTRUE(p > 0.05) && all(m$weights$weight > 0) && all(abs(m$weights$z) > 2)
}
wtxt <- function(m) paste(sprintf("%s %.1f%% (Z %.2f)", m$weights$left,
                                  100*m$weights$weight, m$weights$z), collapse = "  ")

extract <- function(pops, tag, keep = NULL) {
  d <- file.path(OUT, sprintf("f2_%s_%s", tag, STAMP))
  if (dir.exists(d)) stop("refusing to reuse: ", d)
  dir.create(d, recursive = TRUE)
  args <- list(PRE, outdir = d, pops = pops, maxmiss = 0,
               auto_only = TRUE, verbose = FALSE)
  if (!is.null(keep)) args$keepsnps <- keep
  do.call(extract_f2, args)
  f2 <- f2_from_precomp(d, verbose = FALSE)
  n  <- tryCatch(sum(count_snps(f2)), error = function(e) NA_integer_)
  list(f2 = f2, n = n)
}

# ============================== PART A ======================================
cat("##### PART A -- the power curve #####\n")
cat("Model fixed: ", TARGET, "<- Crete_Lasithi_EMBA + Yamnaya_Samara\n")
cat("Only the SNP basis changes.\n\n")

snps <- read.table(paste0(PRE, ".snp"), stringsAsFactors = FALSE)$V1
cat("panel SNPs available:", format(length(snps), big.mark = ","), "\n\n")

POPS_A <- c(TARGET, "Crete_Lasithi_EMBA", "Yamnaya_Samara", right)

runs <- list(
  list(tag = "full",   keep = NULL,   lab = "full basis"),
  list(tag = "500k",   keep = 500000, lab = "~500k subsample"),
  list(tag = "250k",   keep = 250000, lab = "~250k subsample"),
  list(tag = "120k_a", keep = 120000, lab = "~120k subsample (rep A)"),
  list(tag = "120k_b", keep = 120000, lab = "~120k subsample (rep B)")
)

partA <- bind_rows(lapply(runs, function(r) {
  keep <- if (is.null(r$keep)) NULL else sample(snps, min(r$keep, length(snps)))
  X <- extract(POPS_A, paste0("A_", r$tag), keep)
  m <- qpadm(X$f2, left = c("Crete_Lasithi_EMBA", "Yamnaya_Samara"),
             right = right, target = TARGET)
  p <- m$rankdrop$p[1]
  cat(sprintf("  %-24s basis %9s   p = %-11.4g  %s\n      %s\n",
              r$lab, format(X$n, big.mark = ","), p,
              if (adm(m)) "[ADMISSIBLE]" else "[REJECTED]", wtxt(m)))
  tibble(run = r$lab, basis = X$n, p = p,
         w_lasithi = 100*m$weights$weight[1], w_yamnaya = 100*m$weights$weight[2],
         admissible = adm(m))
}))

cat("\n")
print(as.data.frame(partA), row.names = FALSE)
cat("\n>> Paper B's original basis was 122,891 SNPs and reported p = 0.5113.\n")
cat(">> If the ~120k rows above sit near that, the original fit was a POWER\n")
cat(">> artifact -- and so was the 41.9x fold computed from it.\n\n")

# ============================== PART B ======================================
cat("##### PART B -- does ANYTHING fit at full basis? #####\n\n")

SRC <- c("Crete_Lasithi_EMBA", "Crete_Chania_LBA", "Yamnaya_Samara",
         "Khvalynsk_Eneo", "Turkey_N_BIN", "Iran_GanjDareh")

XB <- extract(c(TARGET, SRC, right), "B_grid")
cat("ONE basis for the whole grid:", format(XB$n, big.mark = ","), "SNPs\n")
cat("Every row below is comparable to every other row.\n\n")

pairs <- t(combn(SRC, 2))
partB <- bind_rows(lapply(seq_len(nrow(pairs)), function(i) {
  L <- pairs[i, ]
  m <- tryCatch(qpadm(XB$f2, left = L, right = right, target = TARGET),
                error = function(e) NULL)
  if (is.null(m)) return(tibble())
  p <- m$rankdrop$p[1]
  tibble(source1 = L[1], source2 = L[2], p = p,
         w1 = 100*m$weights$weight[1], z1 = m$weights$z[1],
         w2 = 100*m$weights$weight[2], z2 = m$weights$z[2],
         admissible = adm(m))
})) %>% arrange(desc(p))

print(as.data.frame(partB), row.names = FALSE, digits = 4)

nfit <- sum(partB$admissible)
cat(sprintf("\n>> %d of %d two-source models are ADMISSIBLE at %s SNPs.\n",
            nfit, nrow(partB), format(XB$n, big.mark = ",")))

cat("\n----------------------------------------------------------\n")
if (nfit == 0) {
  cat("VERDICT: POWER REGIME.\n")
  cat("  NOTHING fits Myc_PELOPONNESE at ~1M SNPs. The rejection of\n")
  cat("  Lasithi+Yamnaya says nothing about that relationship in particular.\n")
  cat("  Test 3 is recoverable, but ONLY with a declared basis regime, and\n")
  cat("  Paper B must say that its positive control is basis-conditional.\n")
} else {
  best <- partB[partB$admissible, ][1, ]
  cat("VERDICT: MISFIT.\n")
  cat(sprintf("  %s + %s FITS (p = %.4g) where Lasithi+Yamnaya does not.\n",
              best$source1, best$source2, best$p))
  cat("  The descent relationship Test 3 relies on is NOT the best account of\n")
  cat("  Myc_PELOPONNESE. Test 3 has no validated positive control.\n")
  cat("  PAPER B'S HEADLINE DOES NOT SURVIVE. Report it.\n")
}
cat("----------------------------------------------------------\n\n")

dir.create("manifest", showWarnings = FALSE)
write.csv(partA, "manifest/run_myc_power_partA.csv", row.names = FALSE)
write.csv(partB, "manifest/run_myc_power_partB.csv", row.names = FALSE)
writeLines(c(
  "script: run_myc_power",
  paste("date:", format(Sys.time())),
  "purpose: is the rejection of Myc <- Lasithi + Yamnaya power, or misfit?",
  paste0("partB_basis: ", XB$n),
  paste0("partB_admissible_models: ", nfit, " of ", nrow(partB)),
  "",
  "admissibility = p > 0.05 AND all weights > 0 AND all |Z| > 2.",
  "  run_test3_clean.R omitted the p-value from this test and reported",
  "  [ADMISSIBLE] for a model rejected at p = 8.6e-07. That run is VOID.",
  "",
  "NOTHING here is comparable to 122,891 / 124,362 or any earlier basis."
), "manifest/run_myc_power.txt")

cat("Manifests written.\nDone.\n")
