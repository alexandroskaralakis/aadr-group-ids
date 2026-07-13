# ============================================================================
# run_test3_clean.R
#
# Paper B's Test 3, rebuilt on populations that are populations.
#
# STAGE 1 (gate).  Is the Mycenaean pool one population?
#     qpWave rank-0: Myc_PELOPONNESE (18) vs Myc_ISLANDS (8).
#     A pool is a hypothesis. This one is tested BEFORE it is used.
#
# STAGE 2 (branch, PRE-REGISTERED -- the rule is fixed before the answer):
#     if rank-0 p > 0.05  -> the pool holds. Mycenaean = all 26.
#     if rank-0 p <= 0.05 -> it is a bin. Mycenaean = Myc_PELOPONNESE (18).
#                            The islands are dropped and that is declared.
#
# STAGE 3.  Test 3 on the surviving definition:
#     Mycenaean <- Crete_Lasithi_EMBA + Yamnaya
#       right = six canonical            -> p_base
#       right = six + Georgia_Kotias      -> p_kotias
#       fold  = p_base / p_kotias
#     BOTH models on ONE basis (Kotias is in the extract either way), so the
#     fold is a like-for-like comparison. That is the whole point of the fold.
#
# Run: Rscript run_test3_clean.R 2>&1 | tee manifest/run_test3_clean.log
# ============================================================================

suppressMessages(library(admixtools))
suppressMessages(library(dplyr))

OUT   <- "clean_paperB"
STAMP <- format(Sys.time(), "%Y%m%d_%H%M%S")
PRE   <- file.path(OUT, "panelB")

right   <- readLines(file.path(OUT, "right.txt"));   right   <- right[nzchar(right)]
YAMNAYA <- readLines(file.path(OUT, "yamnaya.txt"))[1]
KOTIAS  <- readLines(file.path(OUT, "kotias.txt"))[1]
SOURCE  <- "Crete_Lasithi_EMBA"

cat("==========================================================\n")
cat(" run_test3_clean.R\n ", format(Sys.time()), "\n")
cat(" admixtools:", as.character(packageVersion("admixtools")), "\n")
cat("==========================================================\n\n")

cat("PRE-REGISTERED, before any model is run:\n")
cat("  P1. The Mycenaean pool is a bin and REJECTS rank 0.        [~55%]\n")
cat("  P2. Mycenaean <- Lasithi + Yamnaya stays ADMISSIBLE, weights\n")
cat("      near the published 86.9 / 13.1.                        [~70%]\n")
cat("  P3. Kotias still degrades it -- fold > 5x.                 [~70%]\n")
cat("  P4. If the fold collapses, Test 3's headline was an artifact of a\n")
cat("      bin, and Paper B loses its exhibit. That is reportable.\n")
cat("  BRANCH RULE (fixed now): rank-0 p > 0.05 -> use all 26.\n")
cat("                           rank-0 p <= 0.05 -> use Peloponnese 18.\n\n")

extract <- function(pops, tag) {
  d <- file.path(OUT, sprintf("f2_%s_%s", tag, STAMP))
  if (dir.exists(d)) stop("refusing to reuse an extract dir: ", d)
  dir.create(d, recursive = TRUE)
  cat("  extract_f2 pops (these, and only these, set the basis):\n")
  for (p in pops) cat("    -", p, "\n")
  extract_f2(PRE, outdir = d, pops = pops, maxmiss = 0,
             auto_only = TRUE, verbose = FALSE)
  f2 <- f2_from_precomp(d, verbose = FALSE)
  n  <- tryCatch(sum(count_snps(f2)), error = function(e) NA_integer_)
  cat("  basis:", format(n, big.mark = ","), "SNPs\n\n")
  list(f2 = f2, n = n, dir = d)
}

adm <- function(m) all(m$weights$weight > 0) && all(abs(m$weights$z) > 2)
wtxt <- function(m) paste(sprintf("%s %.1f%% (Z %.2f)", m$weights$left,
                                  100 * m$weights$weight, m$weights$z),
                          collapse = "   ")

# ---------------------------------------------------------------- STAGE 1 ---
cat("##### STAGE 1 -- is the Mycenaean pool one population? #####\n\n")
X1 <- extract(c("Myc_PELOPONNESE", "Myc_ISLANDS", right), "myc_coherence")
q  <- qpwave(X1$f2, left = c("Myc_PELOPONNESE", "Myc_ISLANDS"), right = right)
r0 <- q$rankdrop[q$rankdrop$f4rank == 0, ]
print(as.data.frame(q$rankdrop))
cat(sprintf("\n  RANK 0: chisq = %.2f  dof = %d  p = %.5g  ->  %s\n\n",
            r0$chisq, r0$dof, r0$p,
            if (r0$p <= 0.05) "REJECTED — THE MYCENAEAN POOL IS A BIN TOO"
            else "not rejected — the pool holds"))

# ---------------------------------------------------------------- STAGE 2 ---
if (r0$p > 0.05) {
  MYC <- "Greece_Mycenaean_ALL"; myc_n <- 26
  cat(">> BRANCH A: pool holds. Mycenaean = all 26 (Peloponnese + islands).\n\n")
  relabel <- TRUE
} else {
  MYC <- "Myc_PELOPONNESE"; myc_n <- 18
  cat(">> BRANCH B: pool rejects. Mycenaean = Peloponnese 18.\n")
  cat(">> The 8 island individuals (Aegina, Paros) are DROPPED and declared.\n\n")
  relabel <- FALSE
}

if (relabel) {
  ind <- read.table(paste0(PRE, ".ind"), stringsAsFactors = FALSE)
  ind$V3[ind$V3 %in% c("Myc_PELOPONNESE", "Myc_ISLANDS")] <- MYC
  P2 <- file.path(OUT, "panelB_merged")
  write.table(ind, paste0(P2, ".ind"), sep = "\t",
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  for (e in c(".geno", ".snp"))
    if (!file.exists(paste0(P2, e)))
      file.symlink(normalizePath(paste0(PRE, e)), paste0(P2, e))
  PRE <- P2
}

# ---------------------------------------------------------------- STAGE 3 ---
cat("##### STAGE 3 -- Test 3 #####\n\n")
cat("  target:", MYC, sprintf("(n = %d)\n", myc_n))
cat("  model :", MYC, "<-", SOURCE, "+", YAMNAYA, "\n\n")

X <- extract(c(MYC, SOURCE, YAMNAYA, KOTIAS, right), "test3")

fit <- function(rt, lab) {
  m <- qpadm(X$f2, left = c(SOURCE, YAMNAYA), right = rt, target = MYC)
  p <- m$rankdrop$p[1]
  cat(sprintf("  %-28s p = %-10.5g  %s\n      %s\n",
              lab, p, if (adm(m)) "[ADMISSIBLE]" else "[inadmissible]", wtxt(m)))
  list(p = p, adm = adm(m), m = m)
}

cat("  Both models below run on the ONE basis printed above. Kotias sits in the\n")
cat("  extract for both, so the only thing that changes is the RIGHT set.\n\n")

base <- fit(right,                  "six canonical refs")
kot  <- fit(c(right, KOTIAS),       "six + Georgia_Kotias")

fold <- base$p / kot$p

cat("\n----------------------------------------------------------\n")
cat(sprintf("  FOLD (p_base / p_kotias) = %.1fx\n", fold))
cat("----------------------------------------------------------\n\n")

verdict <- if (!base$adm) {
  "VOID — the descent relationship is not admissible on this basis.\n       Test 3 cannot run on a control that is not a control."
} else if (fold > 5) {
  "HOLDS — Kotias degrades a true, admissible descent relationship.\n       Test 3 survives the rebuild on coherent populations."
} else {
  "COLLAPSES — the 41.9x fold does not survive coherent populations.\n       Paper B's headline exhibit was an artifact of a bin. Report it."
}
cat("  VERDICT:", verdict, "\n\n")

cat("  Paper B v1.0 reported: p 0.5113 -> 0.0122, fold 41.9x, on\n")
cat("  Greece_Minoan_POOL (6 group IDs) -> Greece_Mycenaean_POOL (5 group IDs),\n")
cat("  at 122,891 SNPs. That basis and those populations are superseded.\n")
cat("  DO NOT put the old and new numbers in one table.\n\n")

# ---- manifest --------------------------------------------------------------
dir.create("manifest", showWarnings = FALSE)
writeLines(c(
  "script: run_test3_clean",
  paste("date:", format(Sys.time())),
  "purpose: Paper B Test 3 rebuilt on coherent AADR-labelled populations",
  "",
  "stage1_mycenaean_coherence:",
  paste0("  left: Myc_PELOPONNESE (18) vs Myc_ISLANDS (8)"),
  paste0("  basis: ", X1$n),
  paste0("  chisq: ", r0$chisq),
  paste0("  dof: ", r0$dof),
  paste0("  rank0_p: ", r0$p),
  paste0("  verdict: ", if (r0$p <= 0.05) "BIN" else "coherent"),
  "",
  "stage3_test3:",
  paste0("  target: ", MYC, " (n=", myc_n, ")"),
  paste0("  model: ", MYC, " <- ", SOURCE, " + ", YAMNAYA),
  paste0("  basis: ", X$n),
  paste0("  populations_in_extract: ", paste(c(MYC, SOURCE, YAMNAYA, KOTIAS, right), collapse=", ")),
  paste0("  p_six_refs: ", base$p, "  admissible: ", base$adm),
  paste0("  p_plus_kotias: ", kot$p, "  admissible: ", kot$adm),
  paste0("  fold: ", fold),
  "",
  "SUPERSEDES Paper B v1.0's 41.9x on 122,891 SNPs, which used",
  "  Greece_Minoan_POOL (not a population, p=3.2e-08) as the SOURCE and",
  "  Greece_Mycenaean_POOL (5 AADR group IDs) as the DESCENDANT.",
  "NOTHING here is comparable to 122,891 / 124,362 or any earlier basis."
), "manifest/run_test3_clean.txt")

cat("Manifest: manifest/run_test3_clean.txt\nDone.\n")
