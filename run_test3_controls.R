# ============================================================================
# run_test3_controls.R
#
# Test 3, on positive controls that are actually controls.
#
# STAGE 1 -- GATE. Does each candidate descent relationship HOLD with the six
#            canonical references? Admissibility here means:
#                p > 0.05  AND  all weights > 0  AND  all |Z| > 2
#            A relationship that does not hold is not a control, and Test 3
#            cannot be run on it. This gate is why the previous run was void:
#            run_test3_clean.R's adm() omitted the p-value, called a model
#            rejected at p = 8.6e-07 "ADMISSIBLE", and then computed a
#            111,374x fold between two dead models.
#
# STAGE 2 -- TEST 3 proper, on whichever controls passed the gate.
#            fold = p(six refs) / p(six + Georgia_Kotias)
#            Both fits on ONE basis; Kotias sits in the extract either way, so
#            the ONLY thing that changes is the right set.
#
# STAGE 3 -- THE SIX MUST PASS. Paper B claims all six canonical references
#            leave a true relationship intact. Verified by leave-one-out and
#            restore: right = five, then right = six. If any of the six shows
#            a fold comparable to Kotias's, Test 3 does not discriminate and
#            the diagnostic is worthless. (The 5->6 vs 6->7 asymmetry is
#            declared, not hidden.)
#
# Run: Rscript run_test3_controls.R 2>&1 | tee manifest/run_test3_controls.log
# ============================================================================

suppressMessages({library(admixtools); library(dplyr)})

OUT   <- "clean_controls"
PRE   <- file.path(OUT, "controls")
STAMP <- format(Sys.time(), "%Y%m%d_%H%M%S")

right  <- readLines(file.path(OUT, "right.txt")); right <- right[nzchar(right)]
KOTIAS <- "Georgia_KotiasKlde_Mesolithic"

cat("==========================================================\n")
cat(" run_test3_controls.R  ", format(Sys.time()), "\n")
cat(" admixtools:", as.character(packageVersion("admixtools")), "\n")
cat("==========================================================\n\n")

cat("PRE-REGISTERED, before any model runs:\n")
cat("  P1. Afanasievo <- Yamnaya (one-source) HOLDS.                [~70%]\n")
cat("  P2. Corded Ware <- Yamnaya + GlobularAmphora HOLDS, ~75/25.  [~65%]\n")
cat("  P3. LBK <- Turkey_N (one-source) REJECTS -- LBK carries a few\n")
cat("      percent WHG. If so, the natural second source is a WHG\n")
cat("      population, but Serbia_IronGates_Mesolithic is IN OUR RIGHT SET.\n")
cat("      That is a reference-into-source conflict, and it is exactly what\n")
cat("      Test 3 exists to detect. Reported, not resolved.            [~60%]\n")
cat("  P4. Kotias degrades whatever holds -- fold > 5x.              [~60%]\n")
cat("  P5. All six canonical references leave it intact.             [~70%]\n")
cat("  FAILURE MODE: if NO control holds, Test 3 cannot be demonstrated on\n")
cat("      real data with this archive, and that is the paper.\n\n")

# --- admissibility, WITH the p-value this time ------------------------------
adm <- function(m) {
  p <- m$rankdrop$p[1]
  isTRUE(p > 0.05) && all(m$weights$weight > 0) && all(abs(m$weights$z) > 2)
}
wtxt <- function(m) paste(sprintf("%s %.1f%% (Z %.2f)", m$weights$left,
                                  100*m$weights$weight, m$weights$z), collapse = "  ")

CONTROLS <- list(
  list(id = "afanasievo",
       target = "Russia_Khakassia_Afanasievo",
       left   = c("Russia_Samara_EBA_Yamnaya"),
       screens = "STEPPE", lit = "Narasimhan 2019; near-clade"),
  list(id = "cordedware",
       target = "Germany_Esperstedt_CordedWare",
       left   = c("Russia_Samara_EBA_Yamnaya", "Czechia_N_GlobularAmphora"),
       screens = "STEPPE", lit = "Haak 2015; canonical ~75/25"),
  list(id = "lbk_slovakia",
       target = "Slovakia_N_LBK",
       left   = c("Turkey_N"),
       screens = "ANATOLIAN", lit = "European farmers <- Anatolian farmers"),
  list(id = "lbk_halberstadt",
       target = "Germany_HalberstadtSonntagsfeld_EN_LBK",
       left   = c("Turkey_N"),
       screens = "ANATOLIAN", lit = "European farmers <- Anatolian farmers")
)

extract <- function(pops, tag) {
  d <- file.path(OUT, sprintf("f2_%s_%s", tag, STAMP))
  if (dir.exists(d)) stop("refusing to reuse: ", d)
  dir.create(d, recursive = TRUE)
  extract_f2(PRE, outdir = d, pops = pops, maxmiss = 0,
             auto_only = TRUE, verbose = FALSE)
  f2 <- f2_from_precomp(d, verbose = FALSE)
  n  <- tryCatch(sum(count_snps(f2)), error = function(e) NA_integer_)
  list(f2 = f2, n = n)
}

results <- list()

for (C in CONTROLS) {
  cat("==========================================================\n")
  cat(sprintf(" %s   [screens %s]\n", C$id, C$screens))
  cat(sprintf(" %s <- %s\n", C$target, paste(C$left, collapse = " + ")))
  cat(sprintf(" literature: %s\n", C$lit))
  cat("==========================================================\n")

  X <- extract(c(C$target, C$left, KOTIAS, right), C$id)
  cat("  basis:", format(X$n, big.mark = ","), "SNPs\n")
  cat("  (Kotias is in the extract for BOTH fits below. Only the right set changes.)\n\n")

  m6 <- qpadm(X$f2, left = C$left, right = right, target = C$target)
  p6 <- m6$rankdrop$p[1]; ok6 <- adm(m6)
  cat(sprintf("  six canonical refs      p = %-11.4g %s\n      %s\n",
              p6, if (ok6) "[HOLDS]" else "[does NOT hold]", wtxt(m6)))

  if (!ok6) {
    cat("\n  >> GATE FAILED. Not a control. Test 3 is not run on it.\n\n")
    results[[C$id]] <- tibble(control = C$id, target = C$target,
                              screens = C$screens, basis = X$n,
                              p_six = p6, holds = FALSE,
                              p_kotias = NA_real_, fold = NA_real_)
    next
  }

  m7 <- qpadm(X$f2, left = C$left, right = c(right, KOTIAS), target = C$target)
  p7 <- m7$rankdrop$p[1]
  fold <- p6 / p7
  cat(sprintf("  six + Georgia_Kotias    p = %-11.4g %s\n      %s\n",
              p7, if (adm(m7)) "[still holds]" else "[BROKEN]", wtxt(m7)))
  cat(sprintf("\n  >> FOLD = %.1fx   %s\n\n", fold,
              if (fold > 5) "-- Kotias degrades a TRUE relationship. Test 3 fires."
              else "-- Kotias does NOT degrade it. Test 3 does not fire."))

  # ---- STAGE 3: do the six themselves pass? ----
  cat("  --- the six canonical references, leave-one-out then restore ---\n")
  six <- bind_rows(lapply(right, function(O) {
    five <- setdiff(right, O)
    a <- qpadm(X$f2, left = C$left, right = five,        target = C$target)
    b <- qpadm(X$f2, left = C$left, right = right,       target = C$target)
    f <- a$rankdrop$p[1] / b$rankdrop$p[1]
    cat(sprintf("    %-30s fold %8.2fx  %s\n", O, f,
                if (f > 5) "<<< FAILS -- behaves like Kotias" else "passes"))
    tibble(control = C$id, ref = O, fold = f, passes = f <= 5)
  }))
  cat("\n")

  results[[C$id]] <- tibble(control = C$id, target = C$target,
                            screens = C$screens, basis = X$n,
                            p_six = p6, holds = TRUE,
                            p_kotias = p7, fold = fold)
  attr(results[[C$id]], "six") <- six
}

# ------------------------------------------------------------------ SUMMARY --
cat("==========================================================\n")
cat(" SUMMARY\n")
cat("==========================================================\n\n")
S <- bind_rows(results)
print(as.data.frame(S), row.names = FALSE, digits = 4)

held <- S %>% filter(holds)
cat("\n----------------------------------------------------------\n")
if (nrow(held) == 0) {
  cat("VERDICT: NO POSITIVE CONTROL HOLDS.\n")
  cat("  Test 3 cannot be demonstrated on real data with this archive.\n")
  cat("  The diagnostic may still be sound in principle; we cannot show it.\n")
  cat("  THAT IS THE PAPER. Report it, and do not claim the method works.\n")
} else {
  fired <- held %>% filter(fold > 5)
  cat(sprintf("VERDICT: %d of %d candidate controls HOLD.\n", nrow(held), nrow(S)))
  if (nrow(fired) > 0) {
    cat(sprintf("  Kotias degrades %d of them (fold > 5x).\n", nrow(fired)))
    cat("  TEST 3 SURVIVES -- on controls established outside the model, with\n")
    cat("  coherent populations, at a full public-data basis.\n")
    if (any(held$screens == "ANATOLIAN" & held$fold > 5))
      cat("  AND it now screens Turkey_N -- the hole Paper B admitted it had.\n")
  } else {
    cat("  But Kotias degrades NONE of them. Test 3 does not fire on any true\n")
    cat("  relationship. The 41.9x was an artifact of the bins. REPORT IT.\n")
  }
}
cat("----------------------------------------------------------\n\n")

dir.create("manifest", showWarnings = FALSE)
write.csv(S, "manifest/run_test3_controls_summary.csv", row.names = FALSE)
sixall <- bind_rows(lapply(results, function(r) attr(r, "six")))
if (nrow(sixall)) write.csv(sixall, "manifest/run_test3_controls_sixrefs.csv", row.names = FALSE)

writeLines(c(
  "script: run_test3_controls",
  paste("date:", format(Sys.time())),
  "purpose: Test 3 on literature-established, coherent positive controls",
  "",
  "admissibility = p > 0.05 AND all weights > 0 AND all |Z| > 2.",
  "  The p-value is REQUIRED. run_test3_clean.R omitted it and reported a",
  "  111,374x fold between two models both rejected at p < 1e-06. VOID.",
  "",
  "SUPERSEDES Paper B v1.0's Test 3 exhibit (41.9x on Mycenaean <- Minoan +",
  "  Yamnaya, 122,891 SNPs). That control is void: both pools are bins, the",
  "  coherent version rejects at every basis, and the pool's apparent fit came",
  "  from post-conquest Chania individuals already carrying steppe ancestry.",
  "",
  "NO usable descendant of Iran_GanjDareh exists in v66.p1. The Iran/CHG",
  "  source cannot be Test-3 screened. Structural and permanent.",
  "",
  "Every population here is an AADR group taken whole. No private filter.",
  "NOTHING here is comparable to 122,891 / 124,362 or any earlier basis."
), "manifest/run_test3_controls.txt")

cat("Manifests written.\nDone.\n")
