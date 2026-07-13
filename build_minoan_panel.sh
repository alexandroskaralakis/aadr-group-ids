#!/usr/bin/env bash
# ============================================================================
# build_minoan_panel.sh   (v2 -- sub-pools from AADR GROUP IDs, not heuristics)
#
# Builds TWO ancients-only EIGENSTRAT panels for the Greece_Minoan_POOL
# coherence test (Paper A, §2.5). NO target genome anywhere in the panel,
# the .ind, or the SNP intersection.
#
#   panel_site : Crete_Lasithi_EMBA (26) vs Crete_Chania_LBA (23)  + 6 right
#   panel_date : Crete_EARLY (>=3600 BP) vs Crete_LATE (<3600 BP)  + 6 right
#
# WHAT CHANGED FROM v1, AND WHY
#   v1 defined sub-pools by substring-matching the Locality field
#   (/lasithi|chania/). That is a heuristic standing in for a fact, and this
#   project has been burned by exactly that four separate times.
#
#   AADR ALREADY SPLITS THESE SIX SITES. The pool's 58 individuals carry six
#   distinct AADR Group IDs. The archive did the right thing; the analyst
#   overrode it with one label. So sub-pools are now defined BY THE ARCHIVE'S
#   OWN GROUP IDs -- non-circular, published, attributable to the archive.
#
#   ONE DECLARED EXCLUSION: Greece_Crete_LBA holds 23 Chania individuals AND
#   2 from Heraklion (kro008.SG, kro009.SG; Skourtanioti/Stockhammer 2023;
#   ~140 km away). It is itself a two-locality bin -- and it trips only the
#   ">=2 localities" condition, not the audit's ">=3 localities" flag, so
#   run_aadr_audit.R scores it ZERO. A second independent instance of the
#   paper's central limitation. The 2 kro are EXCLUDED from the site test so
#   that a rejection means Lasithi != Chania and nothing else.
#
# COLUMN INDICES ARE FACTS, READ FROM THE HEADER, NOT GUESSED:
#   1 = Genetic ID   11 = Date mean in BP   15 = Group ID   16 = Locality
#   The script re-derives them at runtime and DIES if they have moved.
# ============================================================================

set -euo pipefail

# ============================ CONFIG =======================================

AADR_PREFIX="/Users/ell_thales/qpadm_project2/v66.p1_1240K.aadr.patch.PUB"
ANNO="/Users/ell_thales/qpadm_project2/v66.p1_1240K.aadr.PUB.anno"
POOL_IND="/Users/ell_thales/qpadm_project2/postBA_arb.ind"   # IDs only; no genotypes
CONVERTF="/Users/ell_thales/miniconda3/pkgs/admixtools-8.0.2-ha8dfbab_0/bin/convertf"

RIGHT=(
  "Mbuti"
  "Italy_Sicily_Epigravettian"
  "Russia_Kostenki_UP"
  "Papuan"
  "Serbia_IronGates_Mesolithic"
  "Morocco_Iberomaurusian"
)

GRP_LASITHI="Greece_Crete_HgCharalambos_EMBA"
GRP_CHANIA="Greece_Crete_LBA"
EXCLUDE_FROM_CHANIA="kro008.SG kro009.SG"

# Date threshold, BP. Fixed BEFORE any f-statistic. 3600 BP ~ 1650 BCE.
# Knossos falls to the mainland ~1450 BCE (~3400 BP).
DATE_THRESHOLD=3600

OUT="clean_minoan"

# ===========================================================================

echo "=========================================================="
echo " build_minoan_panel.sh  v2"
echo " $(date)"
echo "=========================================================="

for f in "$AADR_PREFIX.geno" "$AADR_PREFIX.snp" "$AADR_PREFIX.ind" "$ANNO" "$POOL_IND"; do
  [ -f "$f" ] || { echo "STOP: missing $f"; exit 1; }
done
[ -x "$CONVERTF" ] || { echo "STOP: convertf not executable at $CONVERTF"; exit 1; }

mkdir -p "$OUT" manifest

# --- 1. the 58 --------------------------------------------------------------
awk '$3=="Greece_Minoan_POOL"{print $1}' "$POOL_IND" | sort > "$OUT/pool_ids.txt"
NPOOL=$(wc -l < "$OUT/pool_ids.txt" | tr -d ' ')
echo ""
echo "[1] Greece_Minoan_POOL members in $(basename "$POOL_IND"): $NPOOL"
[ "$NPOOL" -eq 58 ] || { echo "STOP: expected 58, got $NPOOL."; exit 1; }

# --- 2. columns, re-derived from the header ---------------------------------
COL_ID=$(awk   -F'\t' 'NR==1{for(i=1;i<=NF;i++) if ($i ~ /^Genetic ID/)     {print i; exit}}' "$ANNO")
COL_DATE=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if ($i ~ /Date mean in BP/) {print i; exit}}' "$ANNO")
COL_GRP=$(awk  -F'\t' 'NR==1{for(i=1;i<=NF;i++) if ($i ~ /^Group ID$/)      {print i; exit}}' "$ANNO")
COL_LOC=$(awk  -F'\t' 'NR==1{for(i=1;i<=NF;i++) if ($i ~ /^Locality$/)      {print i; exit}}' "$ANNO")

echo ""
echo "[2] Columns resolved from header:  ID=$COL_ID  Date=$COL_DATE  Group=$COL_GRP  Locality=$COL_LOC"
if [ "$COL_ID" != "1" ] || [ "$COL_DATE" != "11" ] || [ "$COL_GRP" != "15" ] || [ "$COL_LOC" != "16" ]; then
  echo "    STOP: columns differ from the verified 1/11/15/16 layout. Inspect the header."
  exit 1
fi
echo "    Matches the verified layout."

# --- 3. label each of the 58 ------------------------------------------------
awk -F'\t' -v ci="$COL_ID" -v cd="$COL_DATE" -v cg="$COL_GRP" -v cl="$COL_LOC" \
    -v thr="$DATE_THRESHOLD" -v gl="$GRP_LASITHI" -v gc="$GRP_CHANIA" -v ex="$EXCLUDE_FROM_CHANIA" '
  BEGIN { OFS="\t"; n=split(ex, e, " "); for(i=1;i<=n;i++) drop[e[i]]=1 }
  NR==FNR { want[$1]=1; next }
  FNR==1  { next }
  ($ci in want) {
    id=$ci; grp=$cg; loc=$cl; d=$cd+0;
    if      (grp==gl)                  site="Crete_Lasithi_EMBA";
    else if (grp==gc && !(id in drop)) site="Crete_Chania_LBA";
    else                               site="EXCLUDED";
    per = (d >= thr) ? "Crete_EARLY" : "Crete_LATE";
    print id, grp, loc, d, site, per;
  }
' "$OUT/pool_ids.txt" "$ANNO" > "$OUT/pool_labels.tsv"

NMETA=$(wc -l < "$OUT/pool_labels.tsv" | tr -d ' ')
[ "$NMETA" -eq 58 ] || { echo "STOP: matched $NMETA of 58 in the annotation."; exit 1; }

echo ""
echo "[3] AADR Group ID -> site arm  (EYEBALL EVERY ROW):"
awk -F'\t' '{c[$2" -> "$5]++} END{for(k in c) printf "    %3d  %s\n", c[k], k}' "$OUT/pool_labels.tsv" | sort -rn

echo ""
echo "    Site arm:"
awk -F'\t' '{c[$5]++} END{for(k in c) printf "      %-22s %d\n", k, c[k]}' "$OUT/pool_labels.tsv" | sort
echo "    Date arm (threshold ${DATE_THRESHOLD} BP, all 58):"
awk -F'\t' '{c[$6]++} END{for(k in c) printf "      %-22s %d\n", k, c[k]}' "$OUT/pool_labels.tsv" | sort
echo "    Date ranges:"
awk -F'\t' '{k=$6; if(!(k in lo)||$4<lo[k])lo[k]=$4; if($4>hi[k])hi[k]=$4}
  END{for(k in lo) printf "      %-22s %s - %s BP\n", k, lo[k], hi[k]}' "$OUT/pool_labels.tsv" | sort
echo "    Excluded from the site arm:"
awk -F'\t' '$5=="EXCLUDED"{printf "      %-14s %-36s %s\n", $1, $2, $3}' "$OUT/pool_labels.tsv"

# --- 4. relabelled .ind over the FULL AADR ----------------------------------
awk 'NR==FNR{ split($0,a,"\t"); if (a[5]!="EXCLUDED") lab[a[1]]=a[5]; next }
     { g = ($1 in lab) ? lab[$1] : $3; print $1"\t"$2"\t"g }' \
     "$OUT/pool_labels.tsv" "$AADR_PREFIX.ind" > "$OUT/site.ind"

awk 'NR==FNR{ split($0,a,"\t"); lab[a[1]]=a[6]; next }
     { g = ($1 in lab) ? lab[$1] : $3; print $1"\t"$2"\t"g }' \
     "$OUT/pool_labels.tsv" "$AADR_PREFIX.ind" > "$OUT/date.ind"

# --- 5. poplists ------------------------------------------------------------
printf '%s\n' "${RIGHT[@]}" > "$OUT/right.txt"
{ echo "Crete_Lasithi_EMBA"; echo "Crete_Chania_LBA"; printf '%s\n' "${RIGHT[@]}"; } > "$OUT/poplist_site.txt"
{ echo "Crete_EARLY";        echo "Crete_LATE";       printf '%s\n' "${RIGHT[@]}"; } > "$OUT/poplist_date.txt"

echo ""
echo "[5] Poplists: 2 left + 6 right per test -> dof = (2-1)*(6-1) = 5."
echo "    Nothing surplus is in either extract. Every population in an extract"
echo "    sets the basis; none of these is unused by its model."

# --- 6. convertf ------------------------------------------------------------
run_convertf () {
  local ind="$1" poplist="$2" outpre="$3"
  cat > "$OUT/par.$outpre" <<PAR
genotypename:    $AADR_PREFIX.geno
snpname:         $AADR_PREFIX.snp
indivname:       $OUT/$ind
outputformat:    EIGENSTRAT
genotypeoutname: $OUT/$outpre.geno
snpoutname:      $OUT/$outpre.snp
indivoutname:    $OUT/$outpre.ind
poplistname:     $OUT/$poplist
hashcheck:       NO
PAR
  echo ""
  echo "[6] convertf -> $outpre"
  "$CONVERTF" -p "$OUT/par.$outpre" > "$OUT/convertf_$outpre.log" 2>&1 || {
    echo "    STOP: convertf failed. See $OUT/convertf_$outpre.log"; exit 1; }
  echo "    individuals: $(wc -l < "$OUT/$outpre.ind" | tr -d ' ')   panel SNPs: $(wc -l < "$OUT/$outpre.snp" | tr -d ' ')"
  awk '{c[$3]++} END{for(k in c) printf "      %-32s %d\n", k, c[k]}' "$OUT/$outpre.ind" | sort
  local missing=0
  while read -r p; do
    [ -n "$p" ] || continue
    if ! awk -v g="$p" '$3==g{f=1} END{exit !f}' "$OUT/$outpre.ind"; then
      echo "    *** MISSING: $p ***"; missing=1
    fi
  done < "$OUT/$poplist"
  [ "$missing" -eq 0 ] || { echo "    STOP: a poplist population is absent. dof would be silently wrong."; exit 1; }
}

run_convertf site.ind poplist_site.txt panel_site
run_convertf date.ind poplist_date.txt panel_date

# --- 7. target-genome check -------------------------------------------------
echo ""
echo "[7] Target-genome check (must be CLEAN):"
if grep -qil "mydna\|FAM001\|ID001" "$OUT/panel_site.ind" "$OUT/panel_date.ind" 2>/dev/null; then
  echo "    *** FAIL: target genome present ***"; exit 1
fi
echo "    CLEAN — no target genome in either panel."

# --- 8. manifest ------------------------------------------------------------
{
  echo "script: build_minoan_panel (v2)"
  echo "date: $(date)"
  echo "purpose: Paper A 2.5 rebuilt ancients-only; sub-pools from AADR Group IDs"
  echo "aadr_prefix: $AADR_PREFIX"
  echo "anno: $ANNO"
  echo "anno_cols: id=$COL_ID date=$COL_DATE group=$COL_GRP locality=$COL_LOC (header-derived; asserted 1/11/15/16)"
  echo "pool_source: $POOL_IND (58 individual IDs only; no genotypes taken)"
  echo ""
  echo "site_arm:"
  echo "  Crete_Lasithi_EMBA = AADR $GRP_LASITHI (26)"
  echo "  Crete_Chania_LBA   = AADR $GRP_CHANIA minus $EXCLUDE_FROM_CHANIA (23)"
  echo "  exclusion_reason: Greece_Crete_LBA is a TWO-LOCALITY BIN (Chania + Heraklion,"
  echo "    ~140 km apart). kro008.SG/kro009.SG are Heraklion, Skourtanioti/Stockhammer 2023."
  echo "    Excluded so that a rejection means Lasithi != Chania and nothing else."
  echo "  audit_note: this bin trips '>=2 localities' but NOT the audit's '>=3 localities'"
  echo "    flag, so run_aadr_audit.R scores Greece_Crete_LBA ZERO. Second independent"
  echo "    instance of the paper's central limitation (see 2.4)."
  echo ""
  echo "date_arm:"
  echo "  threshold_BP: $DATE_THRESHOLD (fixed before any f-statistic)"
  echo "  all 58 individuals; no exclusions"
  echo ""
  echo "right_set:"
  printf '  - %s\n' "${RIGHT[@]}"
  echo ""
  echo "panel_site_snps: $(wc -l < "$OUT/panel_site.snp" | tr -d ' ')"
  echo "panel_date_snps: $(wc -l < "$OUT/panel_date.snp" | tr -d ' ')"
  echo "NOTE: those are CONVERTF panel sizes, NOT the maxmiss=0 basis."
  echo "  The basis is set at extract_f2. See run_minoan_split_clean.R."
  echo ""
  echo "NOTHING here is comparable to 88,776 / 121,732 / 122,704 / 122,702."
  echo "  Those bases came from a panel containing the target genome, and 88,776"
  echo "  additionally paid for Minoan_OTHER, a population its model never used."
} > manifest/build_minoan_panel.txt

echo ""
echo "=========================================================="
echo " DONE. Manifest: manifest/build_minoan_panel.txt"
echo " Next: Rscript run_minoan_split_clean.R 2>&1 | tee manifest/run_minoan_split_clean.log"
echo "=========================================================="
