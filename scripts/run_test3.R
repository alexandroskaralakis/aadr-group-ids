# run_test3.R — the source-descendant test
#
# For each candidate outgroup O: fit a population independently known to
# descend from source S, with and without O in the right set.
# If O BREAKS a true descent relationship, reject O as a reference.

library(admixtools)
setwd("/Users/ell_thales/qpadm_project2")

PREFIX <- "postBA_west"
DESC   <- "Greece_Mycenaean_POOL"
SRC    <- c("Greece_Minoan_POOL", "Russia_Yamnaya_QC")
TSRC   <- c("Greece_Minoan_POOL", "Iran_GanjDareh_N_QC", "Russia_Eneolithic_Steppe_QC")
OG6    <- c("Mbuti","Italy_Sicily_Epigravettian","Russia_Kostenki_UP",
            "Serbia_IronGates_Mesolithic","Papuan","Morocco_Iberomaurusian")
POPS   <- unique(c(DESC, SRC, TSRC, OG6, "Georgia_Kotias_QC", "MyDNA"))

f2dir <- "f2_test3_audit"
unlink(f2dir, recursive = TRUE); dir.create(f2dir)
extract_f2(PREFIX, f2dir, pops = POPS, maxmiss = 0, blgsize = 0.05, overwrite = TRUE)
f2b <- f2_from_precomp(f2dir)

dir.create("manifest", showWarnings = FALSE)
writeLines(sort(POPS), "manifest/test3_populations.txt")
cat("\nBASIS: populations =", length(POPS), " blocks =", dim(f2b)[3], "\n")

cat("\n######## BASELINE: the true relationship ########\n")
base <- qpadm(f2b, left = SRC, right = OG6, target = DESC)
print(base$weights); print(base$rankdrop[1, ])

cat("\n######## TEST 3: add Kotias to the right set ########\n")
k <- qpadm(f2b, left = SRC, right = c(OG6, "Georgia_Kotias_QC"), target = DESC)
print(k$weights); print(k$rankdrop[1, ])

cat("\n######## Six canonical outgroups: leave-one-out ########\n")
for (O in OG6) {
  r <- qpadm(f2b, left = SRC, right = setdiff(OG6, O), target = DESC)
  cat("--- without", O, "---  p =", signif(r$rankdrop$p[1], 3), "\n")
}

cat("\n######## STANDARD SCREEN 1: Kotias cladal with a SOURCE? ########\n")
cat("Non-rejection = the screen PASSES Kotias.\n\n")
for (S in SRC) {
  w <- qpwave(f2b, left = c("Georgia_Kotias_QC", S), right = OG6)
  cat("  Kotias vs", S, ": p =", signif(w$rankdrop$p[1], 3), "\n")
}

cat("\n######## STANDARD SCREEN 2: nested delta-chi2 on the TARGET ########\n")
cat("On the target CANONICAL 3-source model.\n\n")
a <- qpadm(f2b, left = TSRC, right = OG6, target = "MyDNA")
b <- qpadm(f2b, left = TSRC, right = c(OG6, "Georgia_Kotias_QC"), target = "MyDNA")
cat("  MyDNA without Kotias: p =", signif(a$rankdrop$p[1], 3), "\n")
print(a$weights)
cat("  MyDNA with Kotias:    p =", signif(b$rankdrop$p[1], 3), "\n")
print(b$weights)

cat("\n######## SIDE RESULT: two-source models reject (section 8) ########\n")
for (L in list(c("Greece_Minoan_POOL","Russia_Eneolithic_Steppe_QC"),
               c("Greece_Minoan_POOL","Russia_Yamnaya_QC"))) {
  r <- qpadm(f2b, left = L, right = OG6, target = "MyDNA")
  cat("  ", paste(L, collapse = " + "), ": p =", signif(r$rankdrop$p[1], 3), "\n")
}

cat("\n######## THRESHOLD SENSITIVITY ########\n")
cat("  baseline p =", signif(base$rankdrop$p[1], 4), "\n")
cat("  Test 3   p =", signif(k$rankdrop$p[1], 4), "\n")
cat("  fold drop  =", signif(base$rankdrop$p[1] / k$rankdrop$p[1], 3), "x\n")
cat("  Rejects at alpha = 0.05. Flegontov et al. (2025) use alpha = 0.01.\n")
cat("  This p is BASIS-DEPENDENT: 0.0146 at 124,252 SNPs; 0.0113 at 122,969.\n")
cat("  Adding two populations to the extract moved it 23%.\n")
