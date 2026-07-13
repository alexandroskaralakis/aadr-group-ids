# ============================================================================
# basis_census.R
#
# Scans manifest/ and logs/ for every SNP basis this project has produced, and
# reports them in one table.
#
# WHY THIS EXISTS
#   The project's own rule is ONE AUDITED SNP BASIS PER TABLE. No number from
#   one basis may sit in a table with a number from another. That rule is only
#   enforceable if the bases are enumerated -- and until now they were not.
#   The original census was run interactively and never written to a file.
#   A project about unexamined provenance ran its provenance audit off-script.
#
# WHAT IT DOES
#   Greps every manifest and log for the SNP counts that scripts print, and
#   groups them. It does NOT judge which comparisons are legitimate -- that is
#   a human reading. It surfaces candidates, exactly like run_aadr_audit.R.
#
# Run:  Rscript basis_census.R 2>&1 | tee manifest/basis_census.log
# ============================================================================

suppressMessages({library(tidyverse)})

cat("==========================================================\n")
cat(" basis_census.R\n")
cat(" ", format(Sys.time()), "\n")
cat("==========================================================\n\n")

files <- c(list.files("manifest", full.names = TRUE, pattern = "\\.(txt|log|csv)$"),
           list.files("logs",     full.names = TRUE, pattern = "\\.log$"))
files <- files[file.exists(files)]
cat("Scanning", length(files), "files in manifest/ and logs/\n\n")

# SNP-count patterns the project's scripts actually emit:
#   "  [f2_ai_sites]  9 pops  706 blocks  88,776 SNPs"
#   "basis_kumIN: 122704"
#   "  basis: 950383"
#   "SNP basis (sum of block counts): 950383"
#   "! 1038102 SNPs remain after filtering. 950383 are polymorphic."
pat <- "([0-9][0-9,]{4,})\\s*(SNPs|snps)|basis[_a-zA-Z]*:?\\s*([0-9][0-9,]{4,})|polymorphic"

rows <- map_dfr(files, function(f) {
  ln <- tryCatch(readLines(f, warn = FALSE), error = function(e) character(0))
  hit <- str_which(ln, regex(pat))
  if (!length(hit)) return(tibble())
  tibble(file = basename(f), line = hit, text = str_trim(ln[hit]))
}) %>%
  mutate(n = map(text, ~ str_extract_all(.x, "[0-9][0-9,]{4,}")[[1]])) %>%
  unnest(n) %>%
  mutate(basis = suppressWarnings(as.numeric(str_remove_all(n, ",")))) %>%
  filter(!is.na(basis), basis >= 10000, basis <= 1300000) %>%
  select(file, line, basis, text)

if (!nrow(rows)) { cat("No SNP counts found. Check the patterns.\n"); quit(save = "no") }

cat("######## EVERY DISTINCT BASIS IN THE PROJECT ########\n\n")
census <- rows %>%
  group_by(basis) %>%
  summarise(n_mentions = n(),
            files = paste(sort(unique(file)), collapse = ", "),
            .groups = "drop") %>%
  arrange(desc(basis))
print(census, n = Inf, width = Inf)

cat(sprintf("\n>> %d DISTINCT SNP BASES across %d files.\n",
            nrow(census), n_distinct(rows$file)))
cat(">> The one-basis rule says: no table may mix any two of these.\n")

# --- the clean ancients-only bases, called out ------------------------------
cat("\n######## THE ANCIENTS-ONLY BASES (Paper A 2.5) ########\n")
clean <- census %>% filter(basis > 500000)
if (nrow(clean)) {
  print(clean, n = Inf, width = Inf)
  cat("\n>> These come from panels containing NO target genome, built from AADR\n")
  cat(">> alone. They are the only bases in this project reproducible from\n")
  cat(">> public data. They are NOT comparable to anything below ~124,000.\n")
} else {
  cat("  none found -- has run_minoan_split_clean.R been run?\n")
}

# --- the chip-ceilinged bases ----------------------------------------------
cat("\n######## THE CHIP-CEILINGED BASES (target genome in the panel) ########\n")
chip <- census %>% filter(basis <= 130000)
print(chip, n = Inf, width = Inf)
cat("\n>> Every one of these sits under the ~123,871 SNP ceiling imposed by a\n")
cat(">> 23andMe v5 chip. Any panel built by merging the target genome carries\n")
cat(">> that ceiling -- INCLUDING tests that never use the target.\n")
cat(">> That is the defect corrected in Paper A v1.2 2.5.\n")

write_csv(rows,   "manifest/basis_census_mentions.csv")
write_csv(census, "manifest/basis_census.csv")

cat("\nWritten: manifest/basis_census.csv, manifest/basis_census_mentions.csv\n")
cat("Done.\n")
