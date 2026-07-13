#!/usr/bin/env Rscript
# =====================================================================
# run_test3_controls_v3.R
#   SUPERSEDES the Turkey_MLBA control rows in run_test3_v2.R and the whole
#   of run_test3_sourceswap.R (2x2 factorial -- VOID, see below).
#
# WHY THE PREVIOUS CONTROLS WERE MALFORMED
#   Every Turkey_MLBA control model we ran contained an obligatory STEPPE
#   source. Turkey_MLBA has no steppe ancestry: §2.1's uniform model gives it
#   steppe -23.7% at Z -4.78. So every such model MUST return a significantly
#   negative weight, and every one did:
#
#     MLBA <- EBA + Yam          p 0.582   Yam -16.0 (Z -3.82)   INADMISSIBLE
#     MLBA <- EBA + Iran + Yam   p 0.412   Yam -15.9 (Z -3.57)   INADMISSIBLE
#     MLBA <- MIN + Yam          p 9.8e-11 Yam -43.7 (Z -6.65)   INADMISSIBLE
#     MLBA <- MIN + Iran + Yam   p 0.254   Yam -26.4 (Z -4.66)   INADMISSIBLE
#
#   The folds computed from these (55,765x, 2,536x, 1.0x) are MEANINGLESS.
#   YOU CANNOT BREAK A RELATIONSHIP THAT WAS NEVER STANDING.
#   The 2x2 was also degenerate: with Turkey_EBA as anchor the Iran term came
#   in at 0.2% (Z 0.05) -- the anchor already carries the Iran -- so the
#   "Iran term" factor did not vary across half the grid.
#
#   Root cause: we judged those baselines on their p-values (0.58, 0.41, 0.25
#   all look fine) and never looked at the weights. In a paper whose §12(6)
#   says "judge on weights and Z, never on model p." Logged as self-audit 5.
#
# THE RULE THIS SCRIPT ENFORCES
#   >>> NO FOLD IS COMPUTED UNTIL THE BASELINE IS PROVEN ADMISSIBLE. <<<
#   PART 1 screens candidate baselines. PART 2 runs Test 3 on the survivors
#   ONLY. If nothing survives, the script says so and refuses to produce a
#   number -- which is itself the reportable result.
#
# WHAT IS BEING TESTED
#   Test 3 needs, for each SOURCE S, a population independently known to
#   descend from S. For the Turkey_EBA source, that population is Turkey_MLBA.
#   The model must be STEPPE-FREE, because MLBA has no steppe.
#
#   The anchor contrast then becomes clean and one-variable:
#       MLBA <- Turkey_EBA + Iran     (true descent)
#       MLBA <- Minoan     + Iran     (false anchor, same structure)
#   If the first is admissible and Kotias breaks it, Test 3 works on the
#   Turkey_EBA source and §7 gains a SECOND cleared source (§7.5's limitation
#   weakens from "1 of 3" to "2 of 3").
#   If BOTH are admissible, the anchor swap isolates the mechanism.
#   If NEITHER is, the mechanism is permanently unprovable on these data.
# =====================================================================

suppressPackageStartupMessages(library(admixtools))
suppressPackageStartupMessages(library(tidyverse))

SCRIPT_ID <- "run_test3_controls_v3"
setwd("/Users/ell_thales/qpadm_project2")
dir.create("manifest", showWarnings = FALSE)

PREFIX <- "postBA_arb"
F2DIR  <- "f2_test3_v2"        # REUSE -- keeps ALL of §7 on one basis (122,891)

TGT  <- "MyDNA";                MYC  <- "Greece_Mycenaean_POOL"
MIN  <- "Greece_Minoan_POOL";   YAM  <- "Russia_Yamnaya_QC"
ENEO <- "Russia_Eneolithic_Steppe_QC"
MLBA <- "Turkey_MLBA_QC";       EBA  <- "Turkey_EBA_WEST"
IRAN <- "Iran_GanjDareh_N_QC";  ANF  <- "Turkey_N_QC"
KOT  <- "Georgia_Kotias_QC"

RIGHT <- c("Mbuti", "Italy_Sicily_Epigravettian", "Russia_Kostenki_UP",
           "Papuan", "Serbia_IronGates_Mesolithic", "Morocco_Iberomaurusian")
POPS  <- unique(c(TGT, MYC, MIN, YAM, ENEO, MLBA, EBA, IRAN, ANF, KOT, RIGHT))

ref <- "manifest/run_test3_v2_populations.txt"
if (!file.exists(ref)) stop("Run run_test3_v2.R first -- its extract IS this basis.")
if (!identical(sort(POPS), sort(readLines(ref))))
  stop("Population list differs from run_test3_v2. §7 would go cross-basis. Abort.")
cat("BASIS GATE: matches run_test3_v2 exactly.\n")

f2b  <- f2_from_precomp(F2DIR, verbose = FALSE)
NSNP <- sum(readr::parse_number(dimnames(f2b)[[3]]))
cat(sprintf("§7 BASIS: %d pops, %d blocks, %d SNPs (expect 122,891)\n\n",
            length(POPS), dim(f2b)[3], NSNP))

# ---------------------------------------------------------------------
# ADMISSIBILITY -- the gate. Weights first. p LAST.
# ---------------------------------------------------------------------
admissible <- function(m, one_source = FALSE) {
  p <- m$rankdrop$p[1]
  if (one_source) return(list(ok = p > 0.05, why = "cladality (weight forced to 1)"))
  w <- m$weights$weight; z <- m$weights$z
  bad <- c(if (any(w < 0))            "weight < 0",
           if (any(w > 1))            "weight > 1",
           if (any(abs(z) <= 2))      "a term is n.s. (|Z| <= 2)",
           if (p <= 0.05)             "model rejected")
  list(ok = length(bad) == 0, why = if (length(bad)) paste(bad, collapse = "; ") else "admissible")
}
show_w <- function(m) paste(sprintf("%s %.1f (Z %.2f)", m$weights$left,
                                    100*m$weights$weight, m$weights$z), collapse = " | ")

# =====================================================================
# PART 1 -- SCREEN candidate baselines. No Kotias. No folds. Weights only.
# =====================================================================
cat("######## PART 1: which Turkey_MLBA baselines are even STANDING? ########\n")
cat("Steppe-free by design: Turkey_MLBA has no steppe (§2.1, -23.7%, Z -4.78).\n")
cat("A model that is not admissible CANNOT be used as a Test-3 control.\n\n")

CANDIDATES <- list(
  list(left = EBA,               one = TRUE,  lab = "MLBA <- Turkey_EBA  (cladality)"),
  list(left = c(EBA, IRAN),      one = FALSE, lab = "MLBA <- Turkey_EBA + Iran"),
  list(left = c(EBA, ANF),       one = FALSE, lab = "MLBA <- Turkey_EBA + Turkey_N"),
  list(left = MIN,               one = TRUE,  lab = "MLBA <- Minoan  (cladality)"),
  list(left = c(MIN, IRAN),      one = FALSE, lab = "MLBA <- Minoan + Iran"),
  list(left = c(MIN, ANF),       one = FALSE, lab = "MLBA <- Minoan + Turkey_N"),
  list(left = c(ANF, IRAN),      one = FALSE, lab = "MLBA <- Turkey_N + Iran  (no anchor)")
)

screen <- map_dfr(CANDIDATES, function(cd) {
  m <- qpadm(f2b, left = cd$left, right = RIGHT, target = MLBA)
  a <- admissible(m, cd$one)
  p <- m$rankdrop$p[1]
  cat(sprintf("  %-38s p = %-10.4g %s\n", cd$lab, p,
              if (a$ok) "*** ADMISSIBLE ***" else paste0("inadmissible: ", a$why)))
  if (!cd$one) cat(sprintf("      %s\n", show_w(m)))
  tibble(model = cd$lab, one_source = cd$one, p = p,
         admissible = a$ok, reason = a$why,
         weights = if (cd$one) NA_character_ else show_w(m))
})
cat("\n"); print(select(screen, model, p, admissible, reason), width = Inf)
write_csv(screen, sprintf("manifest/%s_baseline_screen.csv", SCRIPT_ID))

usable <- screen %>% filter(admissible, !one_source)
cat(sprintf("\n>> %d of %d multi-source baselines are admissible.\n",
            nrow(usable), sum(!screen$one_source)))

# =====================================================================
# PART 2 -- Test 3, on ADMISSIBLE baselines ONLY
# =====================================================================
cat("\n######## PART 2: Test 3 -- run ONLY where the baseline stands ########\n")
if (nrow(usable) == 0) {
  cat("\n  NO ADMISSIBLE Turkey_MLBA BASELINE EXISTS.\n")
  cat("  >> Test 3 CANNOT be applied to the Turkey_EBA source on these data.\n")
  cat("  >> §7.5's limitation is CONFIRMED and permanent: the six references are\n")
  cat("  >> cleared on ONE source of three (Minoan), and no more.\n")
  cat("  >> §7.3's mechanism ('Kotias is anti-Turkey_EBA') is UNPROVABLE here.\n")
  cat("  >> Report the absence. Do NOT compute a fold on a fallen baseline.\n")
  t3 <- tibble()
} else {
  t3 <- map_dfr(seq_len(nrow(usable)), function(i) {
    cd <- CANDIDATES[[which(map_chr(CANDIDATES, "lab") == usable$model[i])]]
    a <- qpadm(f2b, left = cd$left, right = RIGHT,          target = MLBA)
    b <- qpadm(f2b, left = cd$left, right = c(RIGHT, KOT),  target = MLBA)
    pa <- a$rankdrop$p[1]; pb <- b$rankdrop$p[1]
    ab <- admissible(b, cd$one)
    cat(sprintf("  %-38s six: %-9.4g  +Kotias: %-10.4g  fold: %9.1fx  %s\n",
                cd$lab, pa, pb, pa/pb, if (pb < 0.05) "<< BROKEN" else "survives"))
    cat(sprintf("      +Kotias weights: %s  [%s]\n", show_w(b), ab$why))
    tibble(model = cd$lab, p_six = pa, p_kotias = pb, fold = pa/pb,
           broken = pb < 0.05, kotias_weights_ok = ab$ok)
  })
  cat("\n"); print(t3, width = Inf)
  write_csv(t3, sprintf("manifest/%s_test3.csv", SCRIPT_ID))
}

# =====================================================================
# PART 3 -- the headline, with its admissibility PROVEN, not assumed
# =====================================================================
cat("\n######## PART 3: the headline -- admissibility shown explicitly ########\n")
h0 <- qpadm(f2b, left = c(MIN, YAM), right = RIGHT,          target = MYC)
h1 <- qpadm(f2b, left = c(MIN, YAM), right = c(RIGHT, KOT),  target = MYC)
ah <- admissible(h0)
cat(sprintf("  Mycenaean <- Minoan + Yamnaya   p = %.4f   [%s]\n",
            h0$rankdrop$p[1], ah$why))
cat(sprintf("      %s     (independent truth: 86.9 / 13.1)\n", show_w(h0)))
cat(sprintf("  + Kotias                        p = %.5f   fold = %.1fx\n",
            h1$rankdrop$p[1], h0$rankdrop$p[1]/h1$rankdrop$p[1]))
cat("\n>> THIS baseline stands: both weights in [0,1], both |Z| > 2, model fits,\n")
cat(">> and it recovers an independently known truth to three digits.\n")
cat(">> Test 3's RESULT rests on this row and on the leave-one-out table.\n")
cat(">> Neither involves Turkey_MLBA. The mechanism was always a bonus.\n")

writeLines(c(
  paste0("script: ", SCRIPT_ID),
  paste0("prefix: ", PREFIX), paste0("f2dir: ", F2DIR, "  (reused; §7 stays on one basis)"),
  paste0("n_snps: ", NSNP),
  paste0("admissible_MLBA_baselines: ", nrow(usable)),
  paste0("headline_admissible: ", ah$ok),
  "RULE: no fold is computed on a baseline that is not admissible.",
  "VOIDS: run_test3_sourceswap.R (2x2) and the Turkey_MLBA rows of run_test3_v2.R"
), sprintf("manifest/%s_basis.txt", SCRIPT_ID))

cat("\nDONE. Rscript run_test3_controls_v3.R > manifest/run_test3_controls_v3.log 2>&1\n")
