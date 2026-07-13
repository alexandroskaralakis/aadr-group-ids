# run_pvalue_pathology.R
#
# Two artefacts, one basis:
#   (a) the HIGHEST-p Bronze Age model in the study requires a NEGATIVE anchor weight
#   (b) a model posts p = 0.09 on a -200% weight and sails past a 0.05 threshold
#
# Model p is a diagnostic of misfit, NOT a criterion of selection.
# Judge on weights and Z.

library(admixtools)
setwd("/Users/ell_thales/qpadm_project2")

PREFIX <- "postBA_west"
TARGET <- "MyDNA"
OG6 <- c("Mbuti","Italy_Sicily_Epigravettian","Russia_Kostenki_UP",
         "Serbia_IronGates_Mesolithic","Papuan","Morocco_Iberomaurusian")

MODELS <- list(
  "A_negative_anchor_highest_p" =
    c("Greece_Minoan_POOL","Turkey_EBA_WEST","Russia_Eneolithic_Steppe_QC"),
  "B_minus200_percent_passes_threshold" =
    c("Greece_Minoan_POOL","Iran_GanjDareh_N_QC","Turkey_PPN_QC"),
  "C_reference_admissible_model" =
    c("Greece_Minoan_POOL","Iran_GanjDareh_N_QC","Russia_Eneolithic_Steppe_QC")
)

POPS <- unique(c(TARGET, OG6, unlist(MODELS)))

f2dir <- "f2_pvalue_pathology"
unlink(f2dir, recursive = TRUE); dir.create(f2dir)
extract_f2(PREFIX, f2dir, pops = POPS, maxmiss = 0, blgsize = 0.05, overwrite = TRUE)
f2b <- f2_from_precomp(f2dir)

dir.create("manifest", showWarnings = FALSE)
writeLines(sort(POPS), "manifest/pvalue_pathology_populations.txt")
cat("\nBASIS: populations =", length(POPS), " blocks =", dim(f2b)[3], "\n")

for (nm in names(MODELS)) {
  L <- MODELS[[nm]]
  cat("\n########", nm, "########\n")
  cat("  ", paste(L, collapse = " + "), "\n")
  r <- tryCatch(qpadm(f2b, left = L, right = OG6, target = TARGET),
                error = function(e) { cat("  FAILED:", conditionMessage(e), "\n"); NULL })
  if (is.null(r)) next
  print(r$weights)
  cat("  model p =", signif(r$rankdrop$p[1], 4), "\n")
  bad <- r$weights$weight < 0 | r$weights$weight > 1
  if (any(bad)) {
    cat("  >>> WEIGHT OUTSIDE [0,1]:",
        paste(sprintf("%s = %.1f%%", r$weights$left[bad], 100*r$weights$weight[bad]),
              collapse = ", "), "\n")
    cat("  >>> A p-value threshold of 0.05 would ADMIT this model.\n")
  }
}

cat("\n######## THE POINT ########\n")
cat("If A has the highest p in the table and an impossible weight,\n")
cat("then p-value ranking is not model selection.\n")
cat("Per Flegontov et al. (2025): out-of-[0,1] weights on a landscape are\n")
cat("SYMMETRY statements -- the target lies beyond a source on a cline.\n")
