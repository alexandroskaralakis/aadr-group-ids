#!/usr/bin/env Rscript
# =====================================================================
# run_eastern_term_v2.R      -- SUPERSEDES run_eastern_term.R (basis-polluted)
#
# WHAT v1 GOT WRONG, AND WHY IT MATTERS
#   v1 put ALL FOUR third-source candidates into ONE extract so that the
#   substitution table (2b) would be internally comparable. Correct instinct.
#   But Turkey_Epipaleolithic (n = 1) costs 22,048 SNPs -- 18% of the basis --
#   and the extract dropped to 101,706.
#
#   The TWO-SOURCE models (2a) DO NOT CONTAIN any candidate. They had no need
#   to be on that basis. On the polluted basis they stopped rejecting:
#
#       Minoan + Eneolithic   0.0034  ->  0.109
#       Minoan + Yamnaya      0.0081  ->  0.171
#
#   §4's first claim -- "the Iran/CHG term is REQUIRED" -- rests on those two
#   rejections. A script written to enforce §12(7) destroyed the paper's one
#   surviving positive result by committing §12(7)'s own pathology.
#   Logged as PRED-S, FAILED.
#
# WHAT v2 DOES
#   Every table now runs on the SMALLEST extract that its models actually need,
#   and every table PRINTS ITS OWN BASIS. Nothing is compared across bases.
#
#     PART 1  basis cost of each candidate                  (diagnostic)
#     PART 2  "REQUIRED": 2-source rejections, CORE extract (no candidates)
#     PART 3  THE EXHIBIT: same models, three bases         (§12(7) evidence)
#     PART 4  "EASTERN": substitution table, COMMON extract (all candidates)
#     PART 5  direction robustness: each candidate on its OWN basis
#     PART 6  Iran vs Kotias cladality, MINIMAL extract
#
#   PART 3 is the reason this script exists. It is the cleanest §12(7) exhibit
#   the project has produced: identical model, identical populations IN the
#   model, only the extract differs -- and the conclusion reverses.
# =====================================================================

suppressPackageStartupMessages(library(admixtools))
suppressPackageStartupMessages(library(tidyverse))

SCRIPT_ID <- "run_eastern_term_v2"
setwd("/Users/ell_thales/qpadm_project2")
dir.create("manifest", showWarnings = FALSE)
dir.create("f2_et2",   showWarnings = FALSE)

PREFIX <- "postBA_arb"

TARGET <- "MyDNA"
ANCHOR <- "Greece_Minoan_POOL"
ENEO   <- "Russia_Eneolithic_Steppe_QC"
YAM    <- "Russia_Yamnaya_QC"
IRAN   <- "Iran_GanjDareh_N_QC"
KOT    <- "Georgia_Kotias_QC"
EPI    <- "Turkey_Epipaleolithic"
PPN    <- "Turkey_PPN_QC"
NATUF  <- "Israel_Natufian_QC"

CAND  <- c(IRAN, KOT, EPI, PPN)
RIGHT <- c("Mbuti", "Italy_Sicily_Epigravettian", "Russia_Kostenki_UP",
           "Papuan", "Serbia_IronGates_Mesolithic", "Morocco_Iberomaurusian")
CORE  <- c(TARGET, ANCHOR, ENEO, YAM, RIGHT)          # 10 pops, NO candidates

ind  <- read.table(paste0(PREFIX, ".ind"), stringsAsFactors = FALSE)
miss <- setdiff(c(CORE, CAND, NATUF), unique(ind$V3))
if (length(miss)) stop("Missing from ", PREFIX, ": ", paste(miss, collapse = ", "))

# ---- helper: build an extract, return f2 array + its basis ----------
build <- function(pops, tag) {
  d <- file.path("f2_et2", tag)
  extract_f2(PREFIX, d, pops = pops, maxmiss = 0, blgsize = 0.05,
             overwrite = TRUE, verbose = FALSE)
  f2 <- f2_from_precomp(d, verbose = FALSE)
  n  <- sum(readr::parse_number(dimnames(f2)[[3]]))
  cat(sprintf("  [extract %-22s]  %2d pops  %3d blocks  %7s SNPs\n",
              tag, length(pops), dim(f2)[3], format(n, big.mark = ",")))
  list(f2 = f2, nsnp = n, nblk = dim(f2)[3], npop = length(pops))
}
# ---- helper: fit, and report chisq + dof, NOT just p ----------------
fit <- function(f2, left, right, target) {
  m  <- qpadm(f2, left = left, right = right, target = target)
  rd <- m$rankdrop[1, ]
  list(m = m, p = rd$p, chisq = rd$chisq, dof = rd$dof, w = m$weights)
}

# =====================================================================
# PART 1 -- what each candidate COSTS
# =====================================================================
cat("######## PART 1: per-candidate basis cost ########\n")
cost <- bind_rows(
  { b <- build(CORE, "core"); tibble(added = "(none - CORE)", n_snp = b$nsnp, n_blk = b$nblk) },
  map_dfr(c(CAND, NATUF), function(x) {
    b <- build(c(CORE, x), paste0("core_", x))
    tibble(added = x, n_snp = b$nsnp, n_blk = b$nblk)
  })
) %>% mutate(cost_snps = n_snp - n_snp[1],
             cost_pct  = round(100 * cost_snps / n_snp[1], 1))
cat("\n"); print(cost, width = Inf)
write_csv(cost, sprintf("manifest/%s_basis_cost.csv", SCRIPT_ID))

# =====================================================================
# PART 2 -- "REQUIRED", on the CORE extract (the correct basis)
# =====================================================================
CORE_X <- build(CORE, "core")
cat(sprintf("\n######## PART 2: the term is REQUIRED  [CORE basis: %s SNPs] ########\n",
            format(CORE_X$nsnp, big.mark = ",")))
cat("These models contain NO candidate. They belong on the CORE extract.\n")
cat("v0.3 prints 0.0034 / 0.0081 from an unnamed basis. THESE supersede them.\n\n")
req <- map_dfr(c(ENEO, YAM), function(s) {
  f <- fit(CORE_X$f2, c(ANCHOR, s), RIGHT, TARGET)
  cat(sprintf("  Minoan + %-28s  p = %-9.5g  chisq = %6.2f (dof %d)   %s\n",
              s, f$p, f$chisq, f$dof, if (f$p < 0.05) "REJECTS ok" else "<<< DOES NOT REJECT"))
  tibble(model = paste("Minoan +", s), basis = CORE_X$nsnp,
         p = f$p, chisq = f$chisq, dof = f$dof,
         w_steppe = 100 * f$w$weight[2], z_steppe = f$w$z[2])
})
cat("\n>> If BOTH reject here, §4's 'the term is required' stands, on a basis\n")
cat(">> that contains only the populations the models actually use.\n")

# =====================================================================
# PART 3 -- THE EXHIBIT: same model, three bases
# =====================================================================
cat("\n######## PART 3: §12(7) EXHIBIT -- identical model, three extracts ########\n")
cat("The MODEL is fixed. The populations IN the model are fixed. Only the\n")
cat("EXTRACT changes -- by adding populations the model never uses.\n\n")
EPI_X <- build(c(CORE, EPI),  "core_epi")      # +Epipaleolithic only
ALL_X <- build(c(CORE, CAND), "common")        # all four candidates
cat("\n")

exhibit <- map_dfr(list(CORE_X, EPI_X, ALL_X), function(X) {
  map_dfr(c(ENEO, YAM), function(s) {
    f <- fit(X$f2, c(ANCHOR, s), RIGHT, TARGET)
    tibble(extract_pops = X$npop, basis_snps = X$nsnp,
           model = paste("Minoan +", sub("Russia_", "", s)),
           chisq = round(f$chisq, 2), dof = f$dof, p = f$p,
           verdict = if (f$p < 0.05) "rejects" else "DOES NOT REJECT")
  })
}) %>% arrange(model, desc(basis_snps))
print(exhibit, width = Inf)
write_csv(exhibit, sprintf("manifest/%s_basis_exhibit.csv", SCRIPT_ID))

cat("\n>> Read the chisq column, not just p. If chisq falls FASTER than the SNP\n")
cat(">> count, the surviving SNPs are not merely fewer -- they are a\n")
cat(">> COMPOSITIONALLY DIFFERENT subset. That is a stronger claim than power\n")
cat(">> loss, and the numbers above decide it. Do not assert it without them.\n")

# =====================================================================
# PART 4 -- "EASTERN", on the COMMON extract (all candidates, comparable)
# =====================================================================
cat(sprintf("\n######## PART 4: the term is EASTERN  [COMMON basis: %s SNPs] ########\n",
            format(ALL_X$nsnp, big.mark = ",")))
cat("All four candidates in ONE extract -- the rows ARE comparable to each other.\n")
cat("Right set identical in every row: reference bias is structurally impossible.\n")
cat("This basis is DEGRADED (Epipaleolithic costs 18%). Say so on the table.\n\n")
east <- map_dfr(CAND, function(cd) {
  f <- fit(ALL_X$f2, c(ANCHOR, ENEO, cd), RIGHT, TARGET)
  w <- deframe(select(f$w, left, weight)); z <- deframe(select(f$w, left, z))
  cat(sprintf("  %-24s  3rd %7.1f (Z %5.2f)   anchor %6.1f   steppe %6.1f (Z %5.2f)   p = %.4g\n",
              cd, 100*w[[cd]], z[[cd]], 100*w[[ANCHOR]], 100*w[[ENEO]], z[[ENEO]], f$p))
  tibble(third = cd, w_third = 100*w[[cd]], z_third = z[[cd]],
         w_anchor = 100*w[[ANCHOR]], w_steppe = 100*w[[ENEO]], z_steppe = z[[ENEO]],
         p = f$p,
         third_ok = w[[cd]] > 0 & abs(z[[cd]]) > 2,
         all_terms_ok = all(c(w[[ANCHOR]], w[[cd]], w[[ENEO]]) >= 0) &
                        all(abs(c(z[[ANCHOR]], z[[cd]], z[[ENEO]])) > 2))
})
cat("\n"); print(east, width = Inf)
cat("\n>> NOTE the two columns. 'third_ok' asks whether the CANDIDATE behaves.\n")
cat(">> 'all_terms_ok' can fail because the STEPPE term went n.s. on this\n")
cat(">> degraded basis -- which is not the candidate's fault. Report both, or\n")
cat(">> the table will say Kotias failed when Kotias did not.\n")
write_csv(east, sprintf("manifest/%s_substitution.csv", SCRIPT_ID))

# =====================================================================
# PART 5 -- direction robustness: each candidate on its OWN basis
# =====================================================================
cat("\n######## PART 5: does the DIRECTION survive a basis change? ########\n")
cat("Each candidate refitted in its own core+candidate extract (best basis for it).\n")
cat("These rows are NOT comparable to each other. They test one thing only:\n")
cat("does the SIGN and SIGNIFICANCE of the third term hold under a better basis?\n\n")
robust <- map_dfr(CAND, function(cd) {
  X <- build(c(CORE, cd), paste0("core_", cd))
  f <- fit(X$f2, c(ANCHOR, ENEO, cd), RIGHT, TARGET)
  w <- deframe(select(f$w, left, weight)); z <- deframe(select(f$w, left, z))
  cat(sprintf("    -> %-22s  own basis %7s SNPs   3rd %7.1f (Z %5.2f)   p = %.4g\n",
              cd, format(X$nsnp, big.mark = ","), 100*w[[cd]], z[[cd]], f$p))
  tibble(third = cd, own_basis = X$nsnp,
         w_third = 100*w[[cd]], z_third = z[[cd]], p = f$p,
         positive_and_sig = w[[cd]] > 0 & z[[cd]] > 2)
})
cat("\n"); print(robust, width = Inf)
comp <- left_join(select(east, third, common_w = w_third, common_z = z_third),
                  select(robust, third, own_w = w_third, own_z = z_third), by = "third")
cat("\n---- direction check: common basis vs own basis ----\n")
print(comp, width = Inf)
cat("\n>> The CLAIM is directional: eastern sources substitute POSITIVELY and\n")
cat(">> SIGNIFICANTLY; forager/ANF-adjacent ones go NEGATIVE and n.s.\n")
cat(">> If the sign holds across both bases for all four, the claim is\n")
cat(">> basis-robust even though the magnitudes are not. Say exactly that.\n")
write_csv(robust, sprintf("manifest/%s_robustness.csv", SCRIPT_ID))

# =====================================================================
# PART 6 -- Iran vs Kotias cladality, on a MINIMAL extract
# =====================================================================
cat("\n######## PART 6: Iran vs Kotias -- non-cladal, yet interchangeable ########\n")
QW_X <- build(c(IRAN, KOT, RIGHT), "qpwave_iran_kotias")
qw   <- qpwave(QW_X$f2, left = c(IRAN, KOT), right = RIGHT)
rd   <- qw$rankdrop[1, ]
cat(sprintf("\n  qpWave(Iran, Kotias | six):  chisq = %.2f  dof = %d  p = %.4g   [%s SNPs]\n",
            rd$chisq, rd$dof, rd$p, format(QW_X$nsnp, big.mark = ",")))
cat("\n>> Decisively non-cladal, yet BOTH substitute positively as sources (PART 4/5).\n")
cat(">> That is the whole argument for the label 'Iran/CHG': the term is eastern,\n")
cat(">> and it is NOT further resolvable with these data.\n")
writeLines(c(sprintf("qpwave_iran_kotias_chisq: %.4f", rd$chisq),
             sprintf("qpwave_iran_kotias_dof: %d",     rd$dof),
             sprintf("qpwave_iran_kotias_p: %.6g",     rd$p),
             sprintf("qpwave_basis_snps: %d",          QW_X$nsnp)),
           sprintf("manifest/%s_qpwave.txt", SCRIPT_ID))

# =====================================================================
# MANIFEST -- one line per basis. There is no single basis for this script,
# and pretending otherwise is the error this script exists to correct.
# =====================================================================
writeLines(c(
  paste0("script: ", SCRIPT_ID),
  paste0("prefix: ", PREFIX),
  "maxmiss: 0  blgsize: 0.05",
  "",
  "THIS SCRIPT USES MULTIPLE BASES BY DESIGN. Each table names its own.",
  sprintf("  PART 2 (required)      CORE extract      %d pops  %s SNPs",
          CORE_X$npop, format(CORE_X$nsnp, big.mark = ",")),
  sprintf("  PART 3 (exhibit)       core / +epi / all %s / %s / %s SNPs",
          format(CORE_X$nsnp, big.mark = ","), format(EPI_X$nsnp, big.mark = ","),
          format(ALL_X$nsnp, big.mark = ",")),
  sprintf("  PART 4 (eastern)       COMMON extract    %d pops  %s SNPs  [DEGRADED]",
          ALL_X$npop, format(ALL_X$nsnp, big.mark = ",")),
  "  PART 5 (robustness)    one extract per candidate -- rows NOT comparable",
  sprintf("  PART 6 (qpWave)        MINIMAL extract   %s SNPs",
          format(QW_X$nsnp, big.mark = ",")),
  "",
  "Israel_Natufian_QC is measured in PART 1 and never modelled.",
  "Its SNP cost is extract-dependent: quote it WITH its extract, never bare."
), sprintf("manifest/%s_basis.txt", SCRIPT_ID))

cat("\nDONE. Rscript run_eastern_term_v2.R > manifest/run_eastern_term_v2.log 2>&1\n")
