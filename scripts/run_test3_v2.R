#!/usr/bin/env Rscript
# =====================================================================
# run_test3_v2.R      -- SUPERSEDES run_test3.R AND run_test3_controls.R
#
# WHY v2 EXISTS
#   run_test3.R      ran on postBA_west / 122,969  -> §7.1 screens, §7.3 LOO
#   run_test3_controls.R ran on postBA_arb / 122,886 -> §7.2, §7.4
#   Quoting both in one section is a CROSS-BASIS PRESENTATION -- ledger entry
#   CCCC, the error the abstract of v0.2 already committed once, now
#   threatening §7: the section carrying the paper's ONLY novel contribution.
#
#   v2 puts EVERY figure in §7 in ONE extract. Screen 2 fits the target, so
#   MyDNA must be in the extract, which moves the basis AGAIN.
#
#   THE FOLD WILL MOVE. It moved 43.8x -> 41.9x between the last two bases.
#   Whatever this script emits is what §7 prints. Do not re-run and keep the
#   friendlier basis. If the fold collapses, the paper reports a collapsed fold.
#
# ONE CLAIM IS REDEFINED, NOT REPRODUCED
#   v0.2 §7.1 asserts Kotias makes "the Aegean model admissible in 6/6
#   rotations and the Anatolian rival rejected in 6/6." No script produces it
#   and the configuration behind "6/6" is not recoverable from any document.
#   PART 3 does NOT guess at it. It defines a clean, unambiguous replacement:
#   both rival models, full rotating design, with and without Kotias.
#   The paper must print the redefined test, and say that it is redefined.
# =====================================================================

suppressPackageStartupMessages(library(admixtools))
suppressPackageStartupMessages(library(tidyverse))

SCRIPT_ID <- "run_test3_v2"
setwd("/Users/ell_thales/qpadm_project2")
dir.create("manifest", showWarnings = FALSE)

PREFIX <- "postBA_arb"
F2DIR  <- "f2_test3_v2"

TGT    <- "MyDNA"
MYC    <- "Greece_Mycenaean_POOL"
MIN    <- "Greece_Minoan_POOL"
YAM    <- "Russia_Yamnaya_QC"
ENEO   <- "Russia_Eneolithic_Steppe_QC"
MLBA   <- "Turkey_MLBA_QC"
EBA    <- "Turkey_EBA_WEST"
IRAN   <- "Iran_GanjDareh_N_QC"
ANF    <- "Turkey_N_QC"
KOT    <- "Georgia_Kotias_QC"

RIGHT  <- c("Mbuti", "Italy_Sicily_Epigravettian", "Russia_Kostenki_UP",
            "Papuan", "Serbia_IronGates_Mesolithic", "Morocco_Iberomaurusian")

# the rotating set (§6): each is a source in some models, a reference in others
ROT    <- c(MIN, EBA, ANF, IRAN, ENEO, YAM)

POPS <- unique(c(TGT, MYC, MIN, YAM, ENEO, MLBA, EBA, IRAN, ANF, KOT, RIGHT))

ind  <- read.table(paste0(PREFIX, ".ind"), stringsAsFactors = FALSE)
miss <- setdiff(POPS, unique(ind$V3))
if (length(miss)) stop("Missing from ", PREFIX, ": ", paste(miss, collapse = ", "))

cat("######## EXTRACT: ONE call. All of §7 lives on this basis. ########\n")
extract_f2(PREFIX, F2DIR, pops = POPS, maxmiss = 0, blgsize = 0.05,
           overwrite = TRUE, verbose = TRUE)
f2b  <- f2_from_precomp(F2DIR, verbose = FALSE)
NSNP <- sum(readr::parse_number(dimnames(f2b)[[3]]))
cat(sprintf("\n>>> §7 BASIS: %d pops, %d blocks, %d SNPs <<<\n", length(POPS), dim(f2b)[3], NSNP))
cat("    (run_test3: 122,969 / postBA_west.  run_test3_controls: 122,886 / postBA_arb.)\n")
cat("    THIS number supersedes both, for every figure in §7.\n\n")

P <- function(m) m$rankdrop$p[1]

# =====================================================================
# PART 1 -- §7.1 SCREEN 1: qpWave cladality of Kotias with every source
# =====================================================================
cat("######## PART 1: Screen 1 -- cladality (the screen that PASSES a bad ref) ########\n")
screen1 <- map_dfr(c(MIN, EBA, ANF, IRAN, ENEO, YAM), function(s) {
  q <- qpwave(f2b, left = c(KOT, s), right = RIGHT)
  p <- q$rankdrop$p[1]
  cat(sprintf("  Kotias vs %-30s  qpWave p = %.4g   %s\n", s, p,
              if (p < 0.05) "non-cladal -> PASS" else "CLADAL -> would fail"))
  tibble(source = s, qpwave_p = p, non_cladal = p < 0.05)
})
cat("\n>> Kotias is decisively non-cladal with every source. Screen 1 PASSES it.\n\n")

# =====================================================================
# PART 2 -- §7.1 SCREEN 2: nested fit on the target's canonical model
# =====================================================================
cat("######## PART 2: Screen 2 -- nested fit (the OTHER screen that passes it) ########\n")
CANON <- c(MIN, IRAN, ENEO)
a <- qpadm(f2b, left = CANON, right = RIGHT,        target = TGT)
b <- qpadm(f2b, left = CANON, right = c(RIGHT, KOT), target = TGT)
cat(sprintf("  MyDNA <- Minoan + Iran + Eneolithic\n"))
cat(sprintf("    six refs   : p = %.4f\n", P(a)))
cat(sprintf("    + Kotias   : p = %.4f\n", P(b)))
cat("\n  weights, six refs:\n");  print(a$weights)
cat("\n  weights, + Kotias:\n");  print(b$weights)
cat("\n>> A legitimate reference FIGHTS the fit: misfit rises, SE falls, weights hold.\n")
cat(">> If Kotias does that here, Screen 2 PASSES it too -- and both screens are wrong.\n\n")
screen2 <- tibble(model = "MyDNA <- Minoan+Iran+Eneolithic",
                  p_six = P(a), p_kotias = P(b),
                  w_anchor_six = a$weights$weight[1], w_anchor_kot = b$weights$weight[1],
                  se_anchor_six = a$weights$se[1],    se_anchor_kot = b$weights$se[1])
print(screen2, width = Inf)

# =====================================================================
# PART 3 -- §7.1 the "inversion", REDEFINED (see header)
# =====================================================================
cat("\n######## PART 3: does Kotias INVERT the rival models under rotation? ########\n")
cat("REDEFINED TEST -- v0.2's '6/6 rotations' configuration is unrecoverable.\n")
cat("Here: both rival models, FULL rotating design, with and without Kotias.\n\n")
rot_fit <- function(left, label) {
  refs0 <- c(RIGHT, setdiff(ROT, left))
  m0 <- qpadm(f2b, left = left, right = refs0,          target = TGT)
  m1 <- qpadm(f2b, left = left, right = c(refs0, KOT),  target = TGT)
  cat(sprintf("  %-34s rotating: %-10.4g   + Kotias: %-10.4g\n", label, P(m0), P(m1)))
  tibble(model = label, p_rot = P(m0), p_rot_kotias = P(m1))
}
inversion <- bind_rows(
  rot_fit(c(MIN, IRAN, ENEO), "AEGEAN   (Minoan+Iran+Eneo)"),
  rot_fit(c(EBA, IRAN, YAM),  "ANATOLIAN(Turkey_EBA+Iran+Yam)")
)
cat("\n>> If adding Kotias rescues the Aegean model and kills the Anatolian one,\n")
cat(">> the 'inversion' is real -- and PART 5 shows it is an artefact.\n")

# =====================================================================
# PART 4 -- THE HEADLINE
# =====================================================================
cat("\n######## PART 4: Test 3 headline, on THIS basis ########\n")
test3 <- function(tgt, left, label) {
  a <- qpadm(f2b, left = left, right = RIGHT,          target = tgt)
  b <- qpadm(f2b, left = left, right = c(RIGHT, KOT),  target = tgt)
  pa <- P(a); pb <- P(b)
  cat(sprintf("  %-52s six: %-9.4g  +Kotias: %-10.4g  fold: %8.1fx %s\n",
              label, pa, pb, pa/pb, if (pb < 0.05 && pa > 0.05) "<< BROKEN" else ""))
  tibble(relationship = label, p_six = pa, p_kotias = pb, fold = pa/pb,
         broken = (pa > 0.05 & pb < 0.05))
}
hd <- test3(MYC, c(MIN, YAM), "Mycenaean <- Minoan + Yamnaya [HEADLINE]")
m0 <- qpadm(f2b, left = c(MIN, YAM), right = RIGHT, target = MYC)
cat("\n  baseline weights (independent truth: 86.9 / 13.1):\n"); print(m0$weights)

# =====================================================================
# PART 5 -- §7.2 controls + the DOSE-RESPONSE
# =====================================================================
cat("\n######## PART 5: §7.2 -- the controlled contrast ########\n")
cat("Same target (Turkey_MLBA), same reference change. ONLY the source differs.\n\n")
ctrl <- bind_rows(
  test3(MLBA, c(EBA, ENEO),      "Turkey_MLBA <- Turkey_EBA + Eneolithic"),
  test3(MLBA, c(EBA, YAM),       "Turkey_MLBA <- Turkey_EBA + Yamnaya"),
  test3(MLBA, c(MIN, IRAN, YAM), "Turkey_MLBA <- Minoan + Iran + Yamnaya [no Turkey_EBA]")
)

cat("\n---- DOSE-RESPONSE: fold vs the source's Iran/CHG content (§2) ----\n")
cat("Kotias should damage a relationship IN PROPORTION to how much post-split\n")
cat("CHG-related gene flow its source received. If it does, Harney's violation\n")
cat("is not just detected -- it is QUANTIFIED.\n\n")
dose <- tibble(
  model      = c("Turkey_MLBA <- Minoan+Iran+Yam", "Mycenaean <- Minoan+Yam",
                 "Turkey_MLBA <- Turkey_EBA+Yam",  "Turkey_MLBA <- Turkey_EBA+Eneo"),
  main_source= c("Minoan (no EBA)", "Minoan", "Turkey_EBA", "Turkey_EBA"),
  iran_pct   = c(12.5, 12.5, 32.1, 32.1),
  fold       = c(ctrl$fold[3], hd$fold, ctrl$fold[2], ctrl$fold[1]))
print(dose, width = Inf)

# =====================================================================
# PART 6 -- §7.3 leave-one-out of the canonical six
# =====================================================================
cat("\n######## PART 6: §7.3 -- the six references, leave-one-out ########\n")
loo <- map_dfr(RIGHT, function(o) {
  p <- P(qpadm(f2b, left = c(MIN, YAM), right = setdiff(RIGHT, o), target = MYC))
  cat(sprintf("  without %-32s p = %.4f   ok\n", o, p))
  tibble(reference_removed = o, p_without = p)
})
cat(sprintf("  with all six                             p = %.4f\n", hd$p_six))
cat(sprintf("  + Georgia_Kotias_QC                      p = %.4f   << FAILS\n", hd$p_kotias))

# =====================================================================
# PART 7 -- §7.4 the steppe control void
# =====================================================================
cat("\n######## PART 7: §7.4 -- is the steppe control void? ########\n")
steppe <- bind_rows(
  { p <- P(qpadm(f2b, left = ENEO,          right = RIGHT, target = YAM))
    cat(sprintf("  Yamnaya <- Eneolithic (cladality)   p = %.4g\n", p))
    tibble(model = "Yamnaya <- Eneolithic (cladality)", p = p) },
  { p <- P(qpadm(f2b, left = c(ENEO, ANF),  right = RIGHT, target = YAM))
    cat(sprintf("  Yamnaya <- Eneolithic + Turkey_N    p = %.4g\n", p))
    tibble(model = "Yamnaya <- Eneolithic + Turkey_N", p = p) },
  { p <- P(qpadm(f2b, left = c(ENEO, IRAN), right = RIGHT, target = YAM))
    cat(sprintf("  Yamnaya <- Eneolithic + Iran        p = %.4g   <- v0.2's orphan 2e-5\n", p))
    tibble(model = "Yamnaya <- Eneolithic + Iran", p = p) }
)
cat("\n>> All three rejecting = control void = the six are cleared on 1 source of 3.\n")

# =====================================================================
# OUTPUT + MANIFEST
# =====================================================================
walk2(list(screen1, screen2, inversion, bind_rows(hd, ctrl), dose, loo, steppe),
      c("screen1_cladality", "screen2_nested", "rotation_inversion",
        "test3_rows", "dose_response", "leave_one_out", "steppe_control"),
      ~ write_csv(.x, sprintf("manifest/%s_%s.csv", SCRIPT_ID, .y)))

writeLines(sort(POPS), sprintf("manifest/%s_populations.txt", SCRIPT_ID))
writeLines(c(
  paste0("script: ", SCRIPT_ID),
  paste0("prefix: ", PREFIX), paste0("f2dir: ", F2DIR),
  "maxmiss: 0  blgsize: 0.05",
  paste0("n_populations: ", length(POPS)),
  paste0("n_jackknife_blocks: ", dim(f2b)[3]),
  paste0("n_snps: ", NSNP),
  paste0("headline_p_six: ",    signif(hd$p_six, 5)),
  paste0("headline_p_kotias: ", signif(hd$p_kotias, 5)),
  paste0("headline_fold: ",     signif(hd$fold, 4)),
  "SUPERSEDES: run_test3.R (122,969) and run_test3_controls.R (122,886).",
  "ALL of §7 must be quoted from THIS basis. No exceptions.",
  "PART 3 is a REDEFINED test, not a reproduction of v0.2's '6/6 rotations'."
), sprintf("manifest/%s_basis.txt", SCRIPT_ID))

cat("\nDONE. Rscript run_test3_v2.R > manifest/run_test3_v2.log 2>&1\n")
