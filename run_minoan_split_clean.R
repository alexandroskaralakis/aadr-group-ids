# ============================================================================
# run_minoan_split_clean.R
#
# Paper A, §2.5 — rebuilt on an ancients-only basis.
#
# Two qpWave rank-0 tests on Greece_Minoan_POOL sub-pools, using panels that
# contain NO target genome. The SNP basis is therefore set by AADR alone and
# is reproducible by anyone with the public archive.
#
# Sub-pools defined from AADR's own Locality and Date fields. Never from
# genetic similarity. Date threshold fixed before any f-statistic was computed.
#
# One-basis rule: each test gets its OWN fresh extract_f2 directory. No cache
# is reused. The two tests are compared to nothing but themselves.
#
# Run:  Rscript run_minoan_split_clean.R 2>&1 | tee manifest/run_minoan_split_clean.log
# ============================================================================

suppressMessages(library(admixtools))

OUT <- "clean_minoan"
STAMP <- format(Sys.time(), "%Y%m%d_%H%M%S")

cat("==========================================================\n")
cat(" run_minoan_split_clean.R\n")
cat(" ", format(Sys.time()), "\n")
cat(" admixtools:", as.character(packageVersion("admixtools")), "\n")
cat("==========================================================\n\n")

right <- readLines(file.path(OUT, "right.txt"))
right <- right[nzchar(right)]

cat("Right set (", length(right), "):\n", sep = "")
for (r in right) cat("  -", r, "\n")
cat("\nExpected dof for a 2-left rank-0 test:", (2 - 1) * (length(right) - 1), "\n\n")

# --- pre-registration, restated in the log --------------------------------
cat("PRE-REGISTERED PREDICTIONS (logged before this script was run):\n")
cat("  1. Basis rises to 400k-900k SNPs (chip ceiling removed).\n")
cat("  2. Date split still rejects, and harder than p = 0.0014.\n")
cat("  3. Site split still rejects.\n")
cat("  4. If the date split WEAKENS, the original rejection was partly a\n")
cat("     property of the 23andMe v5 SNP set. That is reportable either way.\n\n")

# --- the runner -----------------------------------------------------------
run_split <- function(prefix, left, tag) {

  cat("----------------------------------------------------------\n")
  cat(" TEST:", tag, "\n")
  cat(" left :", paste(left, collapse = " vs "), "\n")
  cat(" panel:", prefix, "\n")
  cat("----------------------------------------------------------\n")

  # FRESH directory. Stale f2 cache contamination is a documented failure
  # mode in this project; never reuse an outdir across model changes.
  f2dir <- file.path(OUT, sprintf("f2_%s_%s", tag, STAMP))
  if (dir.exists(f2dir)) stop("outdir already exists — refusing to reuse: ", f2dir)
  dir.create(f2dir, recursive = TRUE)

  pops <- c(left, right)
  cat("\nextract_f2 populations (THESE, and only these, set the basis):\n")
  for (p in pops) cat("  -", p, "\n")
  cat("\n")

  extract_f2(prefix,
             outdir   = f2dir,
             pops     = pops,
             maxmiss  = 0,
             auto_only = TRUE,
             verbose  = TRUE)

  f2 <- f2_from_precomp(f2dir, verbose = FALSE)

  # SNP basis. extract_f2 prints it above; capture it here too if available.
  nsnp <- tryCatch(sum(count_snps(f2)), error = function(e) NA_integer_)
  cat("\nSNP basis (sum of block counts):", nsnp, "\n")
  if (is.na(nsnp))
    cat("  (count_snps unavailable — read the basis from the extract_f2 line above)\n")

  res <- qpwave(f2, left = left, right = right)

  cat("\nqpWave rank table:\n")
  print(as.data.frame(res$rankdrop))

  r0 <- res$rankdrop[res$rankdrop$f4rank == 0, ]
  cat("\n  RANK 0:  chisq =", signif(r0$chisq, 5),
      "  dof =", r0$dof,
      "  p =", signif(r0$p, 5),
      if (r0$p < 0.05) "  -> REJECTED\n" else "  -> not rejected\n")
  cat("\n")

  list(tag = tag, left = paste(left, collapse = " | "), pops = pops,
       f2dir = f2dir, nsnp = nsnp,
       chisq = r0$chisq, dof = r0$dof, p = r0$p)
}

# --- the two tests --------------------------------------------------------
site <- run_split(file.path(OUT, "panel_site"),
                  c("Crete_Lasithi_EMBA", "Crete_Chania_LBA"),
                  "site")

date <- run_split(file.path(OUT, "panel_date"),
                  c("Crete_EARLY", "Crete_LATE"),
                  "date")

# --- summary --------------------------------------------------------------
cat("==========================================================\n")
cat(" SUMMARY — ancients-only basis, no target genome\n")
cat("==========================================================\n\n")

summ <- data.frame(
  split = c("Lasithi vs Chania (AADR group IDs)", "Early vs Late (date, 3600 BP)"),
  chisq = c(site$chisq, date$chisq),
  dof   = c(site$dof,   date$dof),
  p     = c(site$p,     date$p),
  basis = c(site$nsnp,  date$nsnp),
  old_p     = c(0.0073636, 0.0013873),
  old_basis = c(88776, 121732)
)
print(summ, row.names = FALSE)

cat("\nThe old_p / old_basis columns are the values from run_anchor_integrity.R.\n")
cat("They are shown for PROVENANCE ONLY, never as a like-for-like comparison:\n")
cat("  - the old panel was MyDNA-merged, so its SNP set was intersected down\n")
cat("    to a 23andMe v5 chip before any extract ran;\n")
cat("  - the old site extract additionally carried Minoan_OTHER, a population\n")
cat("    its own model never used, which cost ~33,000 SNPs;\n")
cat("  - the old site arm defined sub-pools by Locality substring; this one\n")
cat("    uses AADR's own Group IDs.\n")
cat("Different basis, different sub-pool definition, different test.\n")
cat("DO NOT place both in one table in the paper.\n\n")

for (r in list(site, date)) {
  verdict <- if (r$p < 0.05) "REJECTS" else "DOES NOT REJECT"
  cat(sprintf("  %-6s : %s at p = %.5g on %s SNPs\n",
              r$tag, verdict, r$p,
              ifelse(is.na(r$nsnp), "?", format(r$nsnp, big.mark = ","))))
}

# --- manifest -------------------------------------------------------------
dir.create("manifest", showWarnings = FALSE)
mf <- file.path("manifest", "run_minoan_split_clean.txt")
writeLines(c(
  "script: run_minoan_split_clean",
  paste("date:", format(Sys.time())),
  paste("admixtools:", as.character(packageVersion("admixtools"))),
  "purpose: Paper A 2.5 rebuilt on ancients-only AADR basis; no target genome",
  "panels: clean_minoan/panel_site, clean_minoan/panel_date (EIGENSTRAT, from convertf)",
  "extract: maxmiss=0, fresh outdir per test, no cache reuse",
  "",
  "right_set:",
  paste0("  - ", right),
  "",
  "site_test:",
  paste0("  left: ", site$left),
  paste0("  populations_in_extract: ", paste(site$pops, collapse = ", ")),
  paste0("  basis: ", site$nsnp),
  paste0("  chisq: ", site$chisq),
  paste0("  dof: ", site$dof),
  paste0("  rank0_p: ", site$p),
  paste0("  f2dir: ", site$f2dir),
  "",
  "date_test:",
  paste0("  left: ", date$left),
  paste0("  populations_in_extract: ", paste(date$pops, collapse = ", ")),
  paste0("  basis: ", date$nsnp),
  paste0("  chisq: ", date$chisq),
  paste0("  dof: ", date$dof),
  paste0("  rank0_p: ", date$p),
  paste0("  f2dir: ", date$f2dir),
  "",
  "NOTHING here is comparable to 88,776 / 121,732 / 122,704 / 122,702.",
  "Those bases were set by an extract containing the target genome.",
  "This one is not. The two are separate bases and separate tables."
), mf)

cat("\nManifest:", mf, "\n")
cat("Done.\n")
