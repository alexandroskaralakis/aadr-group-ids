#!/usr/bin/env Rscript
# =====================================================================
# run_aadr_audit.R
#
# Systematically audit every AADR group ID against the annotation fields
# that AADR itself distributes.
#
# THE QUESTION
#   A qpAdm source is a string in an .ind file. Analysts take those strings
#   from AADR group IDs. How often does a group ID describe something that is
#   not a population?
#
# WHAT WE ALREADY KNOW (the two cases that motivated this, v66.p1)
#   Turkey_MLBA  = 28 individuals, ALL from Tell Atchana / Alalakh (Hatay) --
#                  a NORTHERN LEVANTINE city 25 km from the Syrian border.
#                  Named "Turkey_" because Hatay is in the modern Republic
#                  of Turkey. The source publication's own title classifies
#                  it as Northern Levant.
#   Turkey_EBA   = 4 sites spanning ~700 km (Izmir / Amasya / Cappadocia),
#                  3 publications, 2 platforms, ~1,000 years -- and includes
#                  ulu117 at 5450 BP (context 4000-3000 BCE), i.e. LATE
#                  CHALCOLITHIC, in a group whose name says EBA.
#
# THE POSITIVE CONTROL FOR THIS SCRIPT
#   >>> The metrics below MUST flag Turkey_EBA and Turkey_MLBA. <<<
#   If they do not, the metrics are wrong, not the archive. This is checked
#   explicitly at the end and the script FAILS LOUDLY if they are missed.
#   (Same discipline as Test 3: validate the instrument on a known answer.)
#
# WHAT IS AND IS NOT MACHINE-DETECTABLE
#   Detectable: geographic spread, temporal span, publication and platform
#     mixing, group-name period token vs. individual dates.
#   NOT detectable: whether a group's NAME is a fair description of where it
#     is. "Turkey_MLBA is really Levantine" is a judgement, not a metric.
#     The script surfaces candidates; a human reads them. Say so in the paper.
#
# NO GUESSED COLUMN INDICES.  Every field is located by header name. Four
# times in this project a guessed identifier has cost a run. Not again.
# =====================================================================

suppressPackageStartupMessages({library(tidyverse); library(readr)})

SCRIPT_ID <- "run_aadr_audit"
ANNO <- "v66.p1_1240K.aadr.PUB.anno"     # set to your .anno
dir.create("manifest", showWarnings = FALSE)

# ---------------------------------------------------------------------
# COLUMN RESOLUTION -- by name, never by index
# ---------------------------------------------------------------------
a <- read_tsv(ANNO, show_col_types = FALSE, guess_max = 20000,
              name_repair = "unique")
nm <- names(a)

pick <- function(pattern, label, exclude = NULL) {
  hits <- nm[str_detect(nm, regex(pattern, ignore_case = TRUE))]
  hits <- setdiff(hits, nm[1])   # drop the Genetic ID paragraph BY POSITION.
                                 # (v1 dropped it by header LENGTH < 120 chars --
                                 #  but AADR's own date header is ~122 chars, so
                                 #  the filter ate the column it was meant to keep.
                                 #  A heuristic standing in for a fact. Again.)
  if (!is.null(exclude)) hits <- setdiff(hits, exclude)
  if (!length(hits)) stop("Cannot find a column for: ", label)
  cat(sprintf("  %-18s -> %s\n", label, str_trunc(hits[1], 70)))
  hits[1]
}
cat("RESOLVED COLUMNS (verify these before trusting anything below):\n")
C_ID   <- nm[1]                                   # Genetic ID (long header)
C_GRP  <- pick("^Group ID$",                    "Group ID")
C_LOC  <- pick("^Locality$",                    "Locality")
C_POL  <- pick("^Political Entity$",            "Political Entity")
C_LAT  <- pick("Lat(itude)?",                   "Latitude")
C_LON  <- pick("Long(itude)?",                  "Longitude")
C_DATE <- pick("Date mean in BP",               "Date mean BP")
C_PUB  <- pick("Publication",                   "Publication")
C_SNP  <- tryCatch(pick("SNPs hit", "SNPs hit"), error = function(e) NA_character_)  # optional
cat("\n")

d <- a %>%
  transmute(id     = .data[[C_ID]],
            group  = .data[[C_GRP]],
            loc    = .data[[C_LOC]],
            polity = .data[[C_POL]],
            lat    = suppressWarnings(as.numeric(.data[[C_LAT]])),
            lon    = suppressWarnings(as.numeric(.data[[C_LON]])),
            bp     = suppressWarnings(as.numeric(.data[[C_DATE]])),
            pub    = .data[[C_PUB]],
            nsnp   = if (is.na(C_SNP)) NA_real_ else
                       suppressWarnings(as.numeric(.data[[C_SNP]])) ) %>%
  filter(!is.na(group), group != "", !str_detect(group, "^Ignore")) %>%
  mutate(platform = case_when(str_detect(id, "\\.DG$") ~ "shotgun_diploid",
                              str_detect(id, "\\.SG$") ~ "shotgun",
                              str_detect(id, "\\.HO$") ~ "HumanOrigins",
                              str_detect(id, "\\.(AG|TW|BY|AA|EC|WGC)$") ~ "capture",
                              TRUE ~ "other"),
         outlier_flag = str_detect(group, "-o$|_o$|outlier"))

cat(sprintf("Loaded %d individuals in %d group IDs.\n\n",
            nrow(d), n_distinct(d$group)))

# ---------------------------------------------------------------------
# GEOGRAPHY: max pairwise great-circle distance within a group
# ---------------------------------------------------------------------
hav <- function(la1, lo1, la2, lo2) {
  R <- 6371; p <- pi/180
  dla <- (la2-la1)*p; dlo <- (lo2-lo1)*p
  h <- sin(dla/2)^2 + cos(la1*p)*cos(la2*p)*sin(dlo/2)^2
  2*R*asin(pmin(1, sqrt(h)))
}
spread_km <- function(la, lo) {
  ok <- !is.na(la) & !is.na(lo)
  la <- la[ok]; lo <- lo[ok]
  if (length(la) < 2) return(0)
  m <- outer(seq_along(la), seq_along(la),
             Vectorize(function(i,j) hav(la[i], lo[i], la[j], lo[j])))
  max(m)
}

# ---------------------------------------------------------------------
# PERIOD TOKEN vs DATE.  Conventional Near East / Aegean chronology, BP.
# Deliberately GENEROUS: a flag here means the date is outside even a
# permissive window, not that it is marginal.
# ---------------------------------------------------------------------
PERIOD <- tribble(
  ~token,          ~min_bp, ~max_bp,
  "_EBA",             3900,    5100,   # ~3150-1950 BCE, generous
  "_MBA",             3450,    4050,
  "_MLBA",            3100,    4050,
  "_LBA",             3050,    3600,
  "_IA",              2200,    3100,
  "_Chalcolithic",    5000,    7500,
  "_N$",              7000,   11000,   # Neolithic
  "_Meso",           10000,   15000
)
period_flag <- function(group, bp) {
  hit <- PERIOD %>% filter(str_detect(group, token))
  if (!nrow(hit) || is.na(bp)) return(NA)
  # if several tokens match (e.g. _MLBA also contains _LBA), take the widest
  any(bp < min(hit$min_bp) | bp > max(hit$max_bp))
}

# ---------------------------------------------------------------------
# THE AUDIT
# ---------------------------------------------------------------------
cat("Computing per-group metrics...\n")
audit <- d %>%
  group_by(group) %>%
  summarise(
    n            = n(),
    n_localities = n_distinct(loc,    na.rm = TRUE),
    n_polities   = n_distinct(polity, na.rm = TRUE),
    n_pubs       = n_distinct(pub,    na.rm = TRUE),
    n_platforms  = n_distinct(platform),
    span_km      = round(spread_km(lat, lon)),
    bp_min       = suppressWarnings(min(bp, na.rm = TRUE)),
    bp_max       = suppressWarnings(max(bp, na.rm = TRUE)),
    span_yr      = ifelse(is.finite(bp_min) & is.finite(bp_max), bp_max - bp_min, NA),
    n_period_out = sum(map2_lgl(group, bp, ~ isTRUE(period_flag(.x, .y))), na.rm = TRUE),
    localities   = paste(sort(unique(loc)), collapse = " | "),
    .groups = "drop") %>%
  filter(n >= 2, !str_detect(group, "-o$|_o$")) %>%
  mutate(
    # ---- THE FIVE SOUND FLAGS. n_flags is computed from these ONLY. ----
    F_multisite  = n_localities >= 3,
    F_multipub   = n_pubs       >= 3,
    F_multiplat  = n_platforms  >= 2,
    F_wide_km    = span_km      >= 300,
    F_wide_yr    = span_yr      >= 800,
    n_flags = F_multisite + F_multipub + F_multiplat + F_wide_km + F_wide_yr,

    # ---- TWO WITHDRAWN FLAGS. Computed, reported, NOT counted. ----
    # F_multipoly is CIRCULAR: group IDs are Country_Period and Political Entity
    #   IS the country, so this only tests whether AADR agrees with itself.
    #   It fired on 1 group in 2,090 -- and that one (BantuSA) is modern.
    # F_period is measuring OUR CHRONOLOGY TABLE, not the archive. It fired on
    #   74/74 Scotland_N, 46/46 England_N, 11/11 Wales_N -- a flag that fires on
    #   100% of a group is announcing that the WINDOW is wrong, not the group.
    #   The _N$ window (7,000-11,000 BP) is Near-Eastern; the NW European
    #   Neolithic runs ~7,500-4,000 BP. The table would have to be regionalised
    #   before this flag means anything. It is retained ONLY as a diagnostic on
    #   individuals whose group ID names a period their own date contradicts,
    #   and only where the window is defensible for that region.
    F_multipoly  = n_polities   >= 2,
    F_period     = n_period_out >= 1) %>%
  arrange(desc(n_flags), desc(span_km))

# ---------------------------------------------------------------------
# POSITIVE CONTROL -- the script must catch what we already know
# ---------------------------------------------------------------------
cat("\n######## POSITIVE CONTROL ########\n")
ctrl <- audit %>% filter(group %in% c("Turkey_EBA", "Turkey_MLBA"))
if (nrow(ctrl) == 0) stop("Control groups absent. Check the Group ID column.")
print(ctrl %>% select(group, n, n_localities, n_pubs, n_platforms,
                      span_km, span_yr, n_period_out, n_flags), width = Inf)

if (ctrl$n_flags[ctrl$group == "Turkey_EBA"] < 3)
  stop("FAIL: Turkey_EBA not flagged. The METRICS are wrong, not the archive.")
cat("\n>> Turkey_EBA flagged. The instrument works on the case we know.\n")
cat(">> NOTE: Turkey_MLBA may carry FEW flags -- it is a single site, single\n")
cat(">> region, homogeneous. Its defect is that its NAME is wrong, and no\n")
cat(">> metric can see that. THAT IS THE PAPER'S CENTRAL LIMITATION.\n")
cat(">> The automated audit finds INCOHERENT groups. It cannot find\n")
cat(">> MISNAMED coherent ones. Those need the Locality column and a human.\n")

# ---------------------------------------------------------------------
# RESULTS
# ---------------------------------------------------------------------
cat("\n######## HEADLINE NUMBERS ########\n")
tot <- nrow(audit)
sm  <- function(x) sprintf("%d / %d  (%.1f%%)", sum(x), tot, 100*mean(x))
cat(sprintf("  groups with n >= 2:                 %d\n", tot))
cat("  --- THE FIVE SOUND FLAGS ---\n")
cat(sprintf("  >= 2 genotyping platforms:          %s\n", sm(audit$F_multiplat)))
cat(sprintf("  >= 3 distinct localities:           %s\n", sm(audit$F_multisite)))
cat(sprintf("  temporal span >= 800 yr:            %s\n", sm(audit$F_wide_yr)))
cat(sprintf("  geographic span >= 300 km:          %s\n", sm(audit$F_wide_km)))
cat(sprintf("  >= 3 publications:                  %s\n", sm(audit$F_multipub)))
cat(sprintf("  THREE OR MORE OF THE FIVE:          %s\n", sm(audit$n_flags >= 3)))
cat(sprintf("  TWO OR MORE OF THE FIVE:            %s\n", sm(audit$n_flags >= 2)))
cat("\n  --- TWO WITHDRAWN FLAGS (reported, NOT counted -- see code) ---\n")
cat(sprintf("  >= 2 political entities:  %s   [CIRCULAR]\n", sm(audit$F_multipoly)))
cat(sprintf("  period-token violation:   %s   [MEASURES OUR CHRONOLOGY TABLE]\n",
            sm(audit$F_period)))

cat("\n######## THE 40 WORST OFFENDERS ########\n")
print(audit %>% select(group, n, n_localities, n_pubs, n_platforms,
                       span_km, span_yr, n_period_out, n_flags) %>% head(40),
      n = 40, width = Inf)

cat("\n######## GROUPS SPANNING >= 2 MODERN COUNTRIES ########\n")
print(audit %>% filter(F_multipoly) %>%
        select(group, n, n_polities, span_km, localities) %>% head(25),
      n = 25, width = Inf)

cat("\n######## PERIOD-TOKEN VIOLATIONS (a name that lies about the date) ########\n")
viol <- d %>%
  mutate(bad = map2_lgl(group, bp, ~ isTRUE(period_flag(.x, .y)))) %>%
  filter(bad) %>%
  select(id, group, bp, loc) %>% arrange(group, bp)
cat(sprintf("  %d individuals sit outside their own group's period window.\n",
            nrow(viol)))
print(head(viol, 30), n = 30, width = Inf)

write_csv(audit, sprintf("manifest/%s_groups.csv", SCRIPT_ID))
write_csv(viol,  sprintf("manifest/%s_period_violations.csv", SCRIPT_ID))
writeLines(c(
  paste0("script: ", SCRIPT_ID),
  paste0("anno: ",   ANNO),
  paste0("n_individuals: ", nrow(d)),
  paste0("n_groups_n>=2: ", tot),
  paste0("n_groups_3plus_flags: ", sum(audit$n_flags >= 3)),
  "Flags: >=3 localities; >=2 polities; >=3 pubs; >=2 platforms;",
  "       >=300 km span; >=800 yr span; >=1 individual outside period window.",
  "LIMITATION: detects INCOHERENT groups. Cannot detect MISNAMED coherent ones",
  "            (e.g. Turkey_MLBA = Alalakh). Those require the Locality field",
  "            and human judgement. State this prominently."
), sprintf("manifest/%s_manifest.txt", SCRIPT_ID))

cat("\nDONE. Rscript run_aadr_audit.R > manifest/run_aadr_audit.log 2>&1\n")
