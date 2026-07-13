#!/usr/bin/env bash
# ============================================================================
# build_thin_panel.sh
#
# Rebuilds Paper B §3's exhibit on PUBLIC DATA.
#
# WHAT §3 CLAIMS
#   Adding a thin population to `extract_f2` -- a population the MODEL NEVER
#   USES -- collapses the SNP basis and reverses the model's conclusion.
#   §3 demonstrated this with Turkey_Epipaleolithic (n=1) against a model of
#   the author's own genome from Greece_Minoan_POOL. Neither is reproducible:
#   the target is private, and the source is not a population (companion
#   paper, §2.5).
#
# WHAT THIS BUILDS
#   The same demonstration with nothing private and no bins:
#     model:  Germany_Esperstedt_CordedWare <- Yamnaya + Czechia_N_GlobularAmphora
#     right:  six canonical + Georgia_Kotias   (REJECTS: chisq 13.93, dof 5, p=0.016)
#   then add a thin population to the EXTRACT ONLY and watch the basis fall.
#
# WHY IT MATTERS MORE THAN §3 REALISED
#   run_chisq_null.R showed that under RANDOM SNP loss, chi^2 does not fall
#   with N -- it RISES (~1.36x at 17.7% loss; 20/20 replicates). §3's premise
#   ("chi^2 scales approximately with SNP count") is FALSE. But its exhibit is
#   therefore MORE anomalous, not less: a chi^2 that FALLS when the basis
#   shrinks is doing something random loss cannot produce.
#
#   This panel lets us test that on data anyone can download.
#
# THIN POPULATIONS ARE DISCOVERED, NOT GUESSED.
#   AADR annotation column 27 is "SNPs hit on autosomal targets (1240k snpset)".
#   We read it, and pick the genuinely low-coverage groups. No pattern-matching
#   on names, no assumptions.
# ============================================================================

set -euo pipefail
cd ~/qpadm_project2

AADR_PREFIX="/Users/ell_thales/qpadm_project2/v66.p1_1240K.aadr.patch.PUB"
ANNO="/Users/ell_thales/qpadm_project2/v66.p1_1240K.aadr.PUB.anno"
CONVERTF="/Users/ell_thales/miniconda3/pkgs/admixtools-8.0.2-ha8dfbab_0/bin/convertf"

RIGHT=("Mbuti" "Italy_Sicily_Epigravettian" "Russia_Kostenki_UP" "Papuan" \
       "Serbia_IronGates_Mesolithic" "Morocco_Iberomaurusian")

CORE=(
  "Germany_Esperstedt_CordedWare"
  "Russia_Samara_EBA_Yamnaya"
  "Czechia_N_GlobularAmphora"
  "Georgia_KotiasKlde_Mesolithic"
)

OUT="clean_thin"
mkdir -p "$OUT" manifest

echo "=========================================================="
echo " build_thin_panel.sh   $(date)"
echo "=========================================================="

# --- columns, verified -----------------------------------------------------
COL_GRP=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if ($i ~ /^Group ID$/) {print i; exit}}' "$ANNO")
COL_COV=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if ($i ~ /SNPs hit on autosomal targets .*1240k/) {print i; exit}}' "$ANNO")
[ "$COL_GRP" = "15" ] || { echo "STOP: Group ID column moved (got $COL_GRP)"; exit 1; }
[ -n "$COL_COV" ]     || { echo "STOP: could not find the 1240k coverage column"; exit 1; }
echo ""
echo "[1] Group ID = col $COL_GRP;  1240k SNP coverage = col $COL_COV"
echo "    $(head -1 "$ANNO" | cut -f"$COL_COV" | cut -c1-70)"

# --- discover THIN groups: small n, low coverage ---------------------------
echo ""
echo "[2] The thinnest small groups in AADR (n<=3, by median 1240k coverage):"
echo "    These are DISCOVERED from the coverage column, not guessed from names."
awk -F'\t' -v cg="$COL_GRP" -v cc="$COL_COV" '
  FNR==1 { next }
  $cg ~ /_o$|-o$|Ignore|_dup/ { next }
  {
    c = $cc + 0
    if (c > 0) { n[$cg]++; sum[$cg]+=c; if (mn[$cg]==0 || c<mn[$cg]) mn[$cg]=c }
  }
  END {
    for (g in n)
      if (n[g] >= 1 && n[g] <= 3 && sum[g]/n[g] < 120000)
        printf "  %8.0f  n=%-2d  %s\n", sum[g]/n[g], n[g], g
  }' "$ANNO" | sort -n | head -25

echo ""
echo "    (left column = mean SNPs hit on the 1240K panel; the whole panel is 1.23M)"

# --- pick the thin populations by coverage, and assert they exist ----------
# Two, spanning a range of thinness. Chosen on coverage, before any statistic.
THIN=(
  "Turkey_Epipaleolithic"
  "Israel_Natufian"
)

echo ""
echo "[3] thin populations selected for the exhibit:"
for T in "${THIN[@]}"; do
  n=$(awk -F'\t' -v cg="$COL_GRP" -v g="$T" 'FNR>1 && $cg==g' "$ANNO" | wc -l | tr -d ' ')
  cov=$(awk -F'\t' -v cg="$COL_GRP" -v cc="$COL_COV" -v g="$T" \
        'FNR>1 && $cg==g {s+=$cc; k++} END{if(k) printf "%.0f", s/k; else print "0"}' "$ANNO")
  printf "    %-30s n=%-3s  mean 1240k coverage = %s\n" "$T" "$n" "$cov"
  [ "$n" -gt 0 ] || { echo "    STOP: $T not found under that exact Group ID."; \
                      echo "    Pick one from the list in [2] and edit THIN=()."; exit 1; }
done

# --- poplist ---------------------------------------------------------------
{ printf '%s\n' "${CORE[@]}"; printf '%s\n' "${RIGHT[@]}"; printf '%s\n' "${THIN[@]}"; } \
  > "$OUT/poplist.txt"
printf '%s\n' "${RIGHT[@]}" > "$OUT/right.txt"
printf '%s\n' "${THIN[@]}"  > "$OUT/thin.txt"

echo ""
echo "[4] poplist ($(wc -l < "$OUT/poplist.txt" | tr -d ' ') pops):"
sed 's/^/    /' "$OUT/poplist.txt"

# --- convertf --------------------------------------------------------------
cat > "$OUT/par.thin" <<PAR
genotypename:    $AADR_PREFIX.geno
snpname:         $AADR_PREFIX.snp
indivname:       $AADR_PREFIX.ind
outputformat:    EIGENSTRAT
genotypeoutname: $OUT/thin.geno
snpoutname:      $OUT/thin.snp
indivoutname:    $OUT/thin.ind
poplistname:     $OUT/poplist.txt
hashcheck:       NO
PAR

echo ""
echo "[5] convertf (10-20 min)..."
"$CONVERTF" -p "$OUT/par.thin" > "$OUT/convertf.log" 2>&1 || {
  echo "    STOP: convertf failed. See $OUT/convertf.log"; exit 1; }

echo "    individuals: $(wc -l < "$OUT/thin.ind" | tr -d ' ')   panel SNPs: $(wc -l < "$OUT/thin.snp" | tr -d ' ')"
awk '{c[$3]++} END{for(k in c) printf "      %-34s %d\n", k, c[k]}' "$OUT/thin.ind" | sort

miss=0
while read -r p; do
  [ -n "$p" ] || continue
  awk -v g="$p" '$3==g{f=1} END{exit !f}' "$OUT/thin.ind" || { echo "    *** MISSING: $p ***"; miss=1; }
done < "$OUT/poplist.txt"
[ "$miss" -eq 0 ] || { echo "    STOP: a poplist population is absent."; exit 1; }

echo ""
echo "[6] target-genome check (by population label):"
N=$(awk '$3=="MyDNA"' "$OUT/thin.ind" | wc -l | tr -d ' ')
[ "$N" -eq 0 ] && echo "    CLEAN — MyDNA absent." || { echo "    *** FAIL ***"; exit 1; }

{
  echo "script: build_thin_panel"
  echo "date: $(date)"
  echo "purpose: rebuild Paper B section 3's exhibit on public data"
  echo ""
  echo "model (REJECTS at full basis: chisq 13.93, dof 5, p = 0.016):"
  echo "  Germany_Esperstedt_CordedWare <- Russia_Samara_EBA_Yamnaya + Czechia_N_GlobularAmphora"
  echo "  right = six canonical + Georgia_KotiasKlde_Mesolithic"
  echo ""
  echo "thin populations added to the EXTRACT ONLY (never in the model):"
  printf '  - %s\n' "${THIN[@]}"
  echo "  Selected on AADR column $COL_COV (SNPs hit, 1240k panel), before any statistic."
  echo ""
  echo "Every population is an AADR group taken whole. No target genome, no pool."
  echo "SUPERSEDES section 3's exhibit, which targeted the author's genome and used"
  echo "  Greece_Minoan_POOL (not a population, rank-0 p = 3.2e-08) as a source."
  echo "NOTHING here is comparable to 124,362 / 102,314 or any earlier basis."
} > manifest/build_thin_panel.txt

echo ""
echo " DONE. Next: Rscript run_chisq_exhibit.R"
