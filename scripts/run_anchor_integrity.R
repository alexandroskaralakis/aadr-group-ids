#!/usr/bin/env Rscript
# =====================================================================
# run_anchor_integrity.R
#
# The AADR audit turned inward. Three tests on the populations THIS PROJECT
# actually used, after the audit showed that Turkey_N (6 flags), Turkey_PPN
# (5 flags) and Serbia_IronGates_Mesolithic (5 flags) are AADR bins.
#
# TEST 1 -- IS THE MINOAN POOL A POPULATION?
#   Greece_Minoan_POOL = 6 sites, 1,042-year span, and it STRADDLES THE
#   MYCENAEAN TAKEOVER OF CRETE (Knossos falls ~1450 BCE = ~3400 BP):
#       Lasithi   n=26   3933-4177 BP   (~2230-1985 BCE, Early/Middle Minoan)
#       Chania    n=23   3175-3300 BP   (~1350-1225 BCE, LATE MINOAN III)
#       Armenoi   n=1    3300 BP        (LM III cemetery)
#   The project recorded "homogeneous, qpWave rank 0, p = 0.09."
#   p = 0.09 IS NOT HOMOGENEITY. IT IS ONE STEP FROM REJECTION.
#   If the Lasithi/Chania split rejects rank 0, the LAST coherent anchor in
#   this project is a bin spanning a conquest, and the 62.3% Minoan weight is
#   a fit to a bin.
#
# TEST 2 -- IS THE HEADLINE NUMBER CONTAMINATED?
#   Turkey_N_QC contains Kumtepe (1 individual, 5615 BP) among 33 others at
#   7870-8559 BP. It is ~2,400 years younger and LATE CHALCOLITHIC, in a pool
#   called Neolithic. It is ulu117 again, inside the PRIMARY SOURCE of the
#   canonical model (60.0 / 26.4 / 13.6, p = 0.752) -- the one number this
#   project still trusts. Drop it and refit.
#   (Also: the project's own notes say "Turkey_N_QC = Barcin." It is 21/34
#    Barcin. Nine individuals are from Konya and Nigde, 450-550 km away.)
#
# TEST 3 -- IS ONE OF THE SIX OUTGROUPS A SOURCE?
#   Serbia_IronGates_Mesolithic spans 11,465 -> 7,803 BP. Farming reaches the
#   Danube Gorges ~8,000 BP. Lepenski Vir (7,825) and several Vlasac
#   individuals (7,803-7,867) POSTDATE IT. If the late fraction carries
#   Anatolian-farmer ancestry, then a REFERENCE shares drift with a SOURCE
#   (Turkey_N) -- the exact Harney violation -- and Test 3 (the diagnostic)
#   IS BLIND TO IT, because Test 3 cleared the six only against the MINOAN
#   source. Turkey_N has no Test-3 control in this panel.
#
# PRE-REGISTERED, BEFORE RUNNING:
#   PRED-AF  Minoan Lasithi vs Chania REJECTS rank 0 (p < 0.05). The pool is
#            a bin. [If this holds, the project has no coherent anchor at all.]
#   PRED-AG  Dropping Kumtepe moves the canonical weights by < 2 points.
#            One individual in 34 should not matter. [If it DOES matter, the
#            headline number is contaminated and Test 2 is the paper.]
#   PRED-AH  Late Iron Gates carries a significant Turkey_N-related component
#            (weight > 0, Z > 2). One of the six references is a source.
#
#   Sub-pools are defined from AADR's OWN Locality and Date fields.
#   NEVER from genetic similarity -- that would be circular.
#   The 8,000 BP Iron Gates threshold is the published arrival of farming in
#   the Danube Gorges, fixed BEFORE looking at any f-statistic.
#
# BASIS: relabelling creates a NEW PANEL and a NEW BASIS. Two panels are
#   needed (Kumtepe in / out). BOTH SNP counts are printed. Nothing from
#   this script may be quoted alongside 122,804 / 122,891 / 124,362.
# =====================================================================

suppressPackageStartupMessages({library(admixtools); library(tidyverse)})

SCRIPT_ID <- "run_anchor_integrity"
setwd("/Users/ell_thales/qpadm_project2")
dir.create("manifest", showWarnings = FALSE)

PREFIX <- "postBA_arb"
ANNO   <- "v66.p1_1240K.aadr.PUB.anno"

RIGHT <- c("Mbuti", "Italy_Sicily_Epigravettian", "Russia_Kostenki_UP",
           "Papuan", "Serbia_IronGates_Mesolithic", "Morocco_Iberomaurusian")
IRAN  <- "Iran_GanjDareh_N_QC"; ENEO <- "Russia_Eneolithic_Steppe_QC"
ANF   <- "Turkey_N_QC";         TGT  <- "MyDNA"

# ---------------------------------------------------------------------
# BUILD SUB-POOL LABELS FROM THE ANNOTATION. No guessed IDs.
# ---------------------------------------------------------------------
ind <- as_tibble(read.table(paste0(PREFIX, ".ind"), stringsAsFactors = FALSE,
                            col.names = c("id", "sex", "pop")))
# as_tibble() is load-bearing: read.table() returns a base data.frame, and
# dplyr verbs on it stay base data.frames, so print(n = 20) silently matches
# na.print and errors. Assumed a class instead of checking it. Same species of
# error as every other identifier bug in this project.
a   <- read_tsv(ANNO, show_col_types = FALSE, guess_max = 20000, name_repair = "unique")
nm  <- names(a)
gf  <- function(p) { h <- setdiff(nm[str_detect(nm, regex(p, ignore_case = TRUE))], nm[1]); h[1] }
ann <- a %>% transmute(id  = .data[[nm[1]]],
                       loc = .data[[gf("^Locality$")]],
                       bp  = suppressWarnings(as.numeric(.data[[gf("Date mean in BP")]])))

j <- ind %>% left_join(ann, by = "id")
cat(sprintf("Joined annotation for %d / %d individuals.\n\n",
            sum(!is.na(j$loc)), nrow(j)))

MINOAN_SPLIT_BP <- 3600     # between Lasithi (>=3933) and Chania (<=3300)
IG_SPLIT_BP     <- 8000     # farming arrives in the Danube Gorges

# MINOAN + KUMTEPE splits only. Serbia_IronGates_Mesolithic stays WHOLE here,
# because it is a member of RIGHT and must exist as a reference in Tests 1-2.
# (v1 split it globally, so RIGHT named a population that no longer existed.)
j <- j %>% mutate(newpop = case_when(
  pop == "Greece_Minoan_POOL" & str_detect(loc, "Lasithi")  ~ "Minoan_LASITHI",
  pop == "Greece_Minoan_POOL" & str_detect(loc, "Chania")   ~ "Minoan_CHANIA",
  pop == "Greece_Minoan_POOL"                               ~ "Minoan_OTHER",
  pop == ANF & str_detect(loc, "Kumtepe")                   ~ "Turkey_KUMTEPE",
  TRUE ~ pop))

# The Iron Gates split lives on its OWN column, used ONLY by TEST 3.
j <- j %>% mutate(newpop_ig = case_when(
  pop == "Serbia_IronGates_Mesolithic" & bp >= IG_SPLIT_BP  ~ "IronGates_EARLY",
  pop == "Serbia_IronGates_Mesolithic" & bp <  IG_SPLIT_BP  ~ "IronGates_LATE",
  TRUE ~ pop))

cat("SUB-POOLS BUILT FROM THE ANNOTATION:\n")
j %>% filter(newpop != pop) %>% count(pop, newpop) %>% print(n = 20)
cat("\n")

# also an EARLY/LATE Minoan split by DATE, to separate site from period
j <- j %>% mutate(newpop2 = case_when(
  pop == "Greece_Minoan_POOL" & bp >= MINOAN_SPLIT_BP ~ "Minoan_EARLY",
  pop == "Greece_Minoan_POOL" & bp <  MINOAN_SPLIT_BP ~ "Minoan_LATE",
  TRUE ~ newpop))

write_panel <- function(popcol, tag) {
  p <- paste0("panel_", tag)
  out <- j %>% transmute(id, sex, pop = .data[[popcol]])
  write.table(out, paste0(p, ".ind"), quote = FALSE, sep = "\t",
              row.names = FALSE, col.names = FALSE)
  for (e in c(".snp", ".geno"))
    if (!file.exists(paste0(p, e))) file.symlink(paste0(PREFIX, e), paste0(p, e))
  p
}
P_OUT <- write_panel("newpop2", "kumOUT")     # Kumtepe = its own pop -> ANF has 33
j2    <- j %>% mutate(newpop3 = ifelse(newpop2 == "Turkey_KUMTEPE", ANF, newpop2))
j$newpop3 <- j2$newpop3
P_IN  <- write_panel("newpop3", "kumIN")      # Kumtepe folded back -> ANF has 34

bld <- function(prefix, pops, dir) {
  extract_f2(prefix, dir, pops = pops, maxmiss = 0, blgsize = 0.05,
             overwrite = TRUE, verbose = FALSE)
  f2 <- f2_from_precomp(dir, verbose = FALSE)
  n  <- sum(readr::parse_number(dimnames(f2)[[3]]))
  cat(sprintf("  [%s]  %d pops  %d blocks  %s SNPs\n", dir, length(pops),
              dim(f2)[3], format(n, big.mark = ",")))
  list(f2 = f2, n = n)
}
adm <- function(m) { w <- m$weights$weight; z <- m$weights$z
  all(w >= 0 & w <= 1) & all(abs(z) > 2) & m$rankdrop$p[1] > 0.05 }
wshow <- function(m) paste(sprintf("%s %.1f (Z %.2f)", m$weights$left,
                                   100*m$weights$weight, m$weights$z), collapse = " | ")

# =====================================================================
# TEST 1 -- IS THE MINOAN POOL A POPULATION?
# =====================================================================
cat("\n######## TEST 1: Greece_Minoan_POOL -- one population, or a bin? ########\n")
POPS1 <- c("Minoan_LASITHI", "Minoan_CHANIA", "Minoan_OTHER",
           "Minoan_EARLY", "Minoan_LATE", RIGHT)
POPS1 <- intersect(unique(c("Minoan_LASITHI","Minoan_CHANIA","Minoan_OTHER", RIGHT)),
                   unique(j$newpop))
X1 <- bld(write_panel("newpop", "sites"), POPS1, "f2_ai_sites")

cat("\n>> POWER WARNING: splitting the Minoan pool into 3 sub-pools at maxmiss=0\n")
cat(">> costs a large share of the basis (each sub-pool needs complete data).\n")
cat(">> A NON-rejection here is therefore NOT evidence of homogeneity -- it is a\n")
cat(">> weaker test on fewer SNPs. Only a REJECTION is conclusive.\n")
cat("\n-- by SITE: Lasithi (Early/Middle Minoan) vs Chania (Late Minoan III) --\n")
q1 <- qpwave(X1$f2, left = c("Minoan_LASITHI", "Minoan_CHANIA"), right = RIGHT)
print(q1$rankdrop)
p1 <- q1$rankdrop$p[q1$rankdrop$f4rank == 0]
cat(sprintf("\n  rank 0 p = %.4g  ->  %s\n", p1,
            if (p1 > 0.05) "CLADAL: the pool survives"
            else "REJECTED: THE MINOAN POOL IS A BIN"))

cat("\n-- by DATE: pre- vs post-3600 BP (Mycenaean takeover ~3400 BP) --\n")
POPS1b <- intersect(unique(c("Minoan_EARLY","Minoan_LATE", RIGHT)), unique(j$newpop2))
X1b <- bld(P_OUT, POPS1b, "f2_ai_dates")
q1b <- qpwave(X1b$f2, left = c("Minoan_EARLY", "Minoan_LATE"), right = RIGHT)
print(q1b$rankdrop)
p1b <- q1b$rankdrop$p[q1b$rankdrop$f4rank == 0]
cat(sprintf("\n  rank 0 p = %.4g  ->  %s\n", p1b,
            if (p1b > 0.05) "CLADAL" else "REJECTED: the split is TEMPORAL, not just spatial"))

cat("\n-- if it IS a bin: does the LATE fraction carry steppe? --\n")
X1c <- bld(P_OUT, c("Minoan_LATE", "Minoan_EARLY", "Russia_Yamnaya_QC", RIGHT), "f2_ai_lm3")
m1c <- qpadm(X1c$f2, left = c("Minoan_EARLY", "Russia_Yamnaya_QC"),
             right = RIGHT, target = "Minoan_LATE")
cat(sprintf("  Minoan_LATE <- Minoan_EARLY + Yamnaya    p = %.4g\n", m1c$rankdrop$p[1]))
cat(sprintf("    %s   [%s]\n", wshow(m1c), if (adm(m1c)) "ADMISSIBLE" else "inadmissible"))
cat("  >> A significant Yamnaya term = the 'Minoan' pool contains MYCENAEANS.\n")

# =====================================================================
# TEST 2 -- DOES ONE CHALCOLITHIC INDIVIDUAL MOVE THE HEADLINE?
# =====================================================================
cat("\n######## TEST 2: Turkey_N_QC -- with and without Kumtepe (5615 BP) ########\n")
cat("Two panels, therefore TWO BASES. Both printed. This is unavoidable and declared.\n\n")
CANON <- function(X, lab, n) {
  m <- qpadm(X$f2, left = c(ANF, IRAN, ENEO), right = RIGHT, target = TGT)
  w <- deframe(select(m$weights, left, weight)); z <- deframe(select(m$weights, left, z))
  cat(sprintf("  %-28s [%s SNPs]  %.1f / %.1f / %.1f   p = %.4f\n",
              lab, format(X$n, big.mark = ","),
              100*w[[ANF]], 100*w[[IRAN]], 100*w[[ENEO]], m$rankdrop$p[1]))
  cat(sprintf("      Z: %.2f / %.2f / %.2f   [%s]\n",
              z[[ANF]], z[[IRAN]], z[[ENEO]], if (adm(m)) "ADMISSIBLE" else "inadmissible"))
  tibble(panel = lab, n_anf = n, basis = X$n, w_anf = 100*w[[ANF]],
         w_iran = 100*w[[IRAN]], w_steppe = 100*w[[ENEO]], p = m$rankdrop$p[1])
}
POPS2 <- c(TGT, ANF, IRAN, ENEO, RIGHT)
XI <- bld(P_IN,  POPS2, "f2_ai_kumIN")
XO <- bld(P_OUT, POPS2, "f2_ai_kumOUT")
cat("\n")
t2 <- bind_rows(CANON(XI, "Turkey_N_QC (34, Kumtepe IN)",  34),
                CANON(XO, "Turkey_N_QC (33, Kumtepe OUT)", 33))
cat(sprintf("\n  WEIGHT SHIFT: ANF %+.1f   Iran %+.1f   steppe %+.1f  points\n",
            diff(t2$w_anf), diff(t2$w_iran), diff(t2$w_steppe)))
cat("  >> The logged canonical result is 60.0 / 26.4 / 13.6, p = 0.752\n")
cat("  >> on a DIFFERENT basis (122,804). Compare the two rows ABOVE to each\n")
cat("  >> other, never to that. If they differ by >2 points, ONE CHALCOLITHIC\n")
cat("  >> INDIVIDUAL IS MOVING THE PROJECT'S HEADLINE NUMBER.\n")

# =====================================================================
# TEST 3 -- IS AN OUTGROUP ACTUALLY A SOURCE?
# =====================================================================
cat("\n######## TEST 3: is Serbia_IronGates_Mesolithic a reference, or a source? ########\n")
cat("Farming reaches the Danube Gorges ~8,000 BP. The pool spans 11,465-7,803 BP.\n\n")
RIGHT_NOIG <- setdiff(RIGHT, "Serbia_IronGates_Mesolithic")
POPS3 <- c("IronGates_EARLY", "IronGates_LATE", ANF, RIGHT_NOIG)
P_IG  <- write_panel("newpop_ig", "ig")
X3 <- bld(P_IG, POPS3, "f2_ai_ig")

cat("\n-- (a) f4: does the LATE fraction lean toward Anatolian farmers? --\n")
f4a <- f4(X3$f2, "Mbuti", ANF, "IronGates_EARLY", "IronGates_LATE")
print(f4a %>% select(pop1, pop2, pop3, pop4, est, se, z, p))
cat("  >> Z << 0 means Turkey_N shares MORE drift with LATE Iron Gates.\n")

cat("\n-- (b) qpAdm: model the late fraction from the early + farmers --\n")
m3 <- qpadm(X3$f2, left = c("IronGates_EARLY", ANF),
            right = RIGHT_NOIG, target = "IronGates_LATE")
cat(sprintf("  IronGates_LATE <- IronGates_EARLY + Turkey_N   p = %.4g\n", m3$rankdrop$p[1]))
cat(sprintf("    %s   [%s]\n", wshow(m3), if (adm(m3)) "ADMISSIBLE" else "inadmissible"))
cat("\n  >> A POSITIVE, SIGNIFICANT Turkey_N weight means one of the six\n")
cat("  >> REFERENCES carries ancestry from a SOURCE. That is the Harney\n")
cat("  >> violation -- and Test 3 (the diagnostic) NEVER SAW IT, because it\n")
cat("  >> cleared the six only against the MINOAN source. Turkey_N has no\n")
cat("  >> Test-3 control in this panel. Paper B must say so.\n")

# ---------------------------------------------------------------------
write_csv(t2, sprintf("manifest/%s_kumtepe.csv", SCRIPT_ID))
writeLines(c(
  paste0("script: ", SCRIPT_ID),
  "PANELS: panel_sites / panel_kumIN / panel_kumOUT (RELABELLED -- NEW BASES)",
  sprintf("basis_kumIN: %d", XI$n), sprintf("basis_kumOUT: %d", XO$n),
  sprintf("minoan_site_rank0_p: %.5g", p1),
  sprintf("minoan_date_rank0_p: %.5g", p1b),
  "Sub-pools from AADR Locality and Date fields ONLY. Never from genetic similarity.",
  "IronGates threshold 8000 BP = published arrival of farming, fixed BEFORE any f-stat.",
  "NOTHING here is comparable to 122,804 / 122,891 / 124,362."
), sprintf("manifest/%s_basis.txt", SCRIPT_ID))

cat("\nDONE. Rscript run_anchor_integrity.R > manifest/run_anchor_integrity.log 2>&1\n")
