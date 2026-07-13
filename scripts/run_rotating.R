# run_rotating.R — base vs rotating reference-set design
#
# Rotating set: every candidate source not currently in use serves as a reference.
# The six base outgroups are always references.
# Harney et al. (2021) recommend this; this project used the base design throughout.

library(admixtools)
setwd("/Users/ell_thales/qpadm_project2")

PREFIX <- "postBA_arb"   # = postBA_west + Turkey_MLBA_QC
TARGET <- "MyDNA"

OG6 <- c("Mbuti","Italy_Sicily_Epigravettian","Russia_Kostenki_UP",
         "Serbia_IronGates_Mesolithic","Papuan","Morocco_Iberomaurusian")

ROT <- c("Greece_Minoan_POOL","Turkey_EBA_WEST","Turkey_N_QC",
         "Iran_GanjDareh_N_QC","Russia_Eneolithic_Steppe_QC","Russia_Yamnaya_QC")

CTRL <- c("Greece_Mycenaean_POOL","Turkey_MLBA_QC")

POPS <- unique(c(TARGET, OG6, ROT, CTRL))

f2dir <- "f2_rotating_audit"
unlink(f2dir, recursive = TRUE); dir.create(f2dir)
extract_f2(PREFIX, f2dir, pops = POPS, maxmiss = 0, blgsize = 0.05, overwrite = TRUE)
f2b <- f2_from_precomp(f2dir)

dir.create("manifest", showWarnings = FALSE)
writeLines(sort(POPS), "manifest/rotating_populations.txt")
cat("\nBASIS: populations =", length(POPS), " blocks =", dim(f2b)[3], "\n")

fit <- function(tgt, L, R, tag) {
  r <- tryCatch(qpadm(f2b, left = L, right = R, target = tgt),
                error = function(e) NULL)
  if (is.null(r)) { cat("  ", tag, ": FAILED\n"); return(invisible(NULL)) }
  w <- paste(sprintf("%s %.1f%% (Z %.2f)", r$weights$left,
             100*r$weights$weight, r$weights$z), collapse = " | ")
  cat(sprintf("  %-46s p = %-11.3g  %s\n", tag, r$rankdrop$p[1], w))
  invisible(r)
}

MODELS <- list(
  c("Turkey_EBA_WEST","Russia_Yamnaya_QC"),
  c("Turkey_EBA_WEST","Russia_Eneolithic_Steppe_QC"),
  c("Turkey_EBA_WEST","Iran_GanjDareh_N_QC","Russia_Yamnaya_QC"),
  c("Greece_Minoan_POOL","Iran_GanjDareh_N_QC","Russia_Eneolithic_Steppe_QC"),
  c("Greece_Minoan_POOL","Iran_GanjDareh_N_QC","Russia_Yamnaya_QC")
)

cat("\n######## BASE DESIGN: six fixed outgroups ########\n\n")
for (L in MODELS) fit(TARGET, L, OG6, paste(L, collapse = " + "))

cat("\n######## ROTATING DESIGN: unused candidates go right ########\n\n")
for (L in MODELS) {
  R <- unique(c(OG6, setdiff(ROT, L)))
  fit(TARGET, L, R, paste(L, collapse = " + "))
}

cat("\n######## ROTATION POWER CONTROL ########\n")
cat("Mycenaean from its OWN ancestors. Independent audit: 86.9 / 13.1.\n")
cat("Rotation must RECOVER this, and must ANNIHILATE the wrong source.\n\n")
for (L in list(c("Greece_Minoan_POOL","Russia_Yamnaya_QC"),
               c("Turkey_EBA_WEST","Russia_Yamnaya_QC"))) {
  R <- unique(c(OG6, setdiff(ROT, L)))
  fit("Greece_Mycenaean_POOL", L, R, paste(L, collapse = " + "))
}

cat("\n######## ROTATION REJECTS A TRUE RELATIONSHIP ########\n")
cat("Turkey_MLBA from its OWN Bronze Age Anatolian ancestors.\n")
cat("If rotation rejects this, it cannot adjudicate a disputed case.\n\n")
MLBA <- list(
  c("Turkey_EBA_WEST","Iran_GanjDareh_N_QC","Russia_Yamnaya_QC"),
  c("Turkey_EBA_WEST","Iran_GanjDareh_N_QC"),
  c("Turkey_EBA_WEST")
)
for (L in MLBA) {
  R <- unique(c(OG6, setdiff(ROT, L)))
  fit("Turkey_MLBA_QC", L, R, paste(L, collapse = " + "))
}
cat("\n  -- same, with Turkey_N removed from the rotation --\n")
for (L in MLBA) {
  R <- unique(c(OG6, setdiff(ROT, c(L, "Turkey_N_QC"))))
  fit("Turkey_MLBA_QC", L, R, paste(L, collapse = " + "))
}

cat("\n######## MATCHED DEGREES OF FREEDOM ########\n")
cat("At equal source counts, does the ranking of the two anchors flip?\n\n")
cat("2-source:\n")
for (L in list(c("Greece_Minoan_POOL","Russia_Yamnaya_QC"),
               c("Turkey_EBA_WEST","Russia_Yamnaya_QC"))) {
  R <- unique(c(OG6, setdiff(ROT, L)))
  fit(TARGET, L, R, paste(L, collapse = " + "))
}
cat("3-source:\n")
for (L in list(c("Greece_Minoan_POOL","Iran_GanjDareh_N_QC","Russia_Yamnaya_QC"),
               c("Turkey_EBA_WEST","Iran_GanjDareh_N_QC","Russia_Yamnaya_QC"))) {
  R <- unique(c(OG6, setdiff(ROT, L)))
  fit(TARGET, L, R, paste(L, collapse = " + "))
}
