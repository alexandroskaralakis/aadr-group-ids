# run_composition_decomp.R
#
# The single axis separating the two anchors is Iran/CHG content.
# Decompose every relevant population on ONE basis as Turkey_N + Iran.
# NOTE: several of these models are REJECTED. Their point estimates are
# indicative only and must be labelled as such wherever quoted.

library(admixtools)
setwd("/Users/ell_thales/qpadm_project2")

PREFIX <- "postBA_arb"   # carries Turkey_MLBA_QC
OG6 <- c("Mbuti","Italy_Sicily_Epigravettian","Russia_Kostenki_UP",
         "Serbia_IronGates_Mesolithic","Papuan","Morocco_Iberomaurusian")

SRC2 <- c("Turkey_N_QC","Iran_GanjDareh_N_QC")
SRC3 <- c("Turkey_N_QC","Iran_GanjDareh_N_QC","Russia_Eneolithic_Steppe_QC")

TARGETS2 <- c("Greece_Minoan_POOL","Greece_Mycenaean_POOL",
              "Turkey_EBA_WEST","Turkey_MLBA_QC")

POPS <- unique(c(TARGETS2, "MyDNA", SRC3, OG6))

f2dir <- "f2_composition"
unlink(f2dir, recursive = TRUE); dir.create(f2dir)
extract_f2(PREFIX, f2dir, pops = POPS, maxmiss = 0, blgsize = 0.05, overwrite = TRUE)
f2b <- f2_from_precomp(f2dir)

dir.create("manifest", showWarnings = FALSE)
writeLines(sort(POPS), "manifest/composition_populations.txt")
cat("\nBASIS: populations =", length(POPS), " blocks =", dim(f2b)[3], "\n")

show <- function(tgt, L) {
  r <- tryCatch(qpadm(f2b, left = L, right = OG6, target = tgt),
                error = function(e) NULL)
  if (is.null(r)) { cat(sprintf("  %-26s FAILED\n", tgt)); return(invisible(NULL)) }
  w <- paste(sprintf("%s %.1f%% (Z %.2f)", r$weights$left,
             100*r$weights$weight, r$weights$z), collapse = " | ")
  p <- r$rankdrop$p[1]
  cat(sprintf("  %-26s p = %-10.3g %s%s\n", tgt, p, w,
      ifelse(p < 0.05, "   <<< MODEL REJECTED - indicative only", "")))
  invisible(r)
}

cat("\n######## TWO-SOURCE: Turkey_N + Iran ########\n\n")
for (t in TARGETS2) show(t, SRC2)

cat("\n######## THREE-SOURCE (target needs steppe) ########\n\n")
show("MyDNA", SRC3)

cat("\n######## THE POINT ########\n")
cat("The anchors differ chiefly in ONE dimension: internal Iran/CHG fraction.\n")
cat("That single axis is the whole story. The six-outgroup panel is blind on it.\n")
