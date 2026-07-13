#!/usr/bin/env bash
# ============================================================================
# build_paperB_panel.sh
#
# Ancients-only panel for the Paper B rebuild. NO target genome.
#
# WHY THIS EXISTS
#   Paper B's headline is Test 3 on the positive control
#       Mycenaean <- Minoan + Yamnaya
#   Paper A has now shown Greece_Minoan_POOL is not a population
#   (qpWave rank-0, p = 3.2e-08), and Paper A's Recommendation 2 says nothing
#   can descend from such a pool. Discovery then showed Greece_Mycenaean_POOL
#   is ALSO a bin -- five AADR Group IDs. The positive control had a bin on
#   both sides.
#
# WHAT CHANGES
#   SOURCE      Greece_Minoan_POOL (6 group IDs, 1,042 yr, straddles the
#               Mycenaean conquest)
#            -> Greece_Crete_HgCharalambos_EMBA  (Lasithi, 26, ONE AADR group,
#               3,933-4,177 BP, entirely PRE-conquest)
#
#   DESCENDANT  Greece_Mycenaean_POOL (5 group IDs)
#            -> split into two sub-pools, and TESTED before use:
#                 Myc_PELOPONNESE = Peloponnese_LBA + PalaceofNestor_BA
#                                 + Argolid_LBA + Achaea_LBA        (18)
#                 Myc_ISLANDS     = Greece_LBA (Aegina + Paros)      (8)
#               A pool is a hypothesis. This one gets tested.
#
#   Yamnaya (Russia_Samara_EBA_Yamnaya, 28) and Kotias
#   (Georgia_KotiasKlde_Mesolithic, 2) are single-group and coherent.
#   They are used under their AADR names, unaltered.
#
# Columns verified: 1=Genetic ID  11=Date mean BP  15=Group ID  16=Locality
# ============================================================================

set -euo pipefail
cd ~/qpadm_project2

AADR_PREFIX="/Users/ell_thales/qpadm_project2/v66.p1_1240K.aadr.patch.PUB"
ANNO="/Users/ell_thales/qpadm_project2/v66.p1_1240K.aadr.PUB.anno"
CONVERTF="/Users/ell_thales/miniconda3/pkgs/admixtools-8.0.2-ha8dfbab_0/bin/convertf"

RIGHT=(
  "Mbuti"
  "Italy_Sicily_Epigravettian"
  "Russia_Kostenki_UP"
  "Papuan"
  "Serbia_IronGates_Mesolithic"
  "Morocco_Iberomaurusian"
)

SOURCE_GRP="Greece_Crete_HgCharalambos_EMBA"     # -> Crete_Lasithi_EMBA (26)
MYC_MAINLAND=("Greece_Peloponnese_LBA" "Greece_PalaceofNestor_BA" \
              "Greece_Argolid_LBA" "Greece_Achaea_LBA")   # -> Myc_PELOPONNESE (18)
MYC_ISLAND=("Greece_LBA")                                  # -> Myc_ISLANDS (8)
YAMNAYA="Russia_Samara_EBA_Yamnaya"                        # 28, coherent
KOTIAS="Georgia_KotiasKlde_Mesolithic"                     # 2, one site, coherent

OUT="clean_paperB"

echo "=========================================================="
echo " build_paperB_panel.sh   $(date)"
echo "=========================================================="

for f in "$AADR_PREFIX.geno" "$AADR_PREFIX.snp" "$AADR_PREFIX.ind" "$ANNO"; do
  [ -f "$f" ] || { echo "STOP: missing $f"; exit 1; }
done
[ -x "$CONVERTF" ] || { echo "STOP: convertf not at $CONVERTF"; exit 1; }
mkdir -p "$OUT" manifest

# --- columns, re-derived ----------------------------------------------------
COL_GRP=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if ($i ~ /^Group ID$/) {print i; exit}}' "$ANNO")
[ "$COL_GRP" = "15" ] || { echo "STOP: Group ID column moved from 15 (got $COL_GRP)"; exit 1; }
echo ""
echo "[1] Group ID column = $COL_GRP (verified)"

# --- relabel the AADR .ind --------------------------------------------------
MAIN=$(IFS='|'; echo "${MYC_MAINLAND[*]}")
ISL=$(IFS='|';  echo "${MYC_ISLAND[*]}")

# Restrict to the individuals actually in postBA_arb -- i.e. the ones that
# survived this project's QC. AADR's full groups are larger (34 Lasithi, 31
# Peloponnese); the surplus were excluded upstream. Taking AADR's groups whole
# would silently discard that QC. We keep it, and DECLARE it: these populations
# are "AADR group X, restricted to the individuals retained by this project's
# QC" -- an unscripted filter, and named as such in the manifest.
awk '$3=="Greece_Minoan_POOL" || $3=="Greece_Mycenaean_POOL" {print $1}' \
  postBA_arb.ind | sort > "$OUT/qc_kept.txt"
echo "    QC-retained individuals available: $(wc -l < "$OUT/qc_kept.txt" | tr -d ' ')"

awk -F'\t' -v cg="$COL_GRP" -v src="$SOURCE_GRP" -v main="$MAIN" -v isl="$ISL" '
  BEGIN{ n=split(main,a,"|"); for(i=1;i<=n;i++) M[a[i]]=1
         m=split(isl,b,"|");  for(i=1;i<=m;i++) I[b[i]]=1 }
  NR==FNR { keep[$1]=1; next }
  FNR==1  { next }
  ($1 in keep) {
    g=$cg
    if      (g==src)  lab[$1]="Crete_Lasithi_EMBA"
    else if (g in M)  lab[$1]="Myc_PELOPONNESE"
    else if (g in I)  lab[$1]="Myc_ISLANDS"
  }
  END{ for (k in lab) print k "\t" lab[k] }
' "$OUT/qc_kept.txt" "$ANNO" > "$OUT/relabel.tsv"

echo ""
echo "[2] relabelled individuals:"
awk -F'\t' '{c[$2]++} END{for(k in c) printf "    %-22s %d\n", k, c[k]}' "$OUT/relabel.tsv" | sort

# sanity: the counts discovery reported
for pair in "Crete_Lasithi_EMBA:26" "Myc_PELOPONNESE:18" "Myc_ISLANDS:8"; do
  p="${pair%%:*}"; want="${pair##*:}"
  got=$(awk -F'\t' -v p="$p" '$2==p' "$OUT/relabel.tsv" | wc -l | tr -d ' ')
  [ "$got" -eq "$want" ] || { echo "    STOP: $p expected $want, got $got"; exit 1; }
done
echo "    counts match discovery. good."

awk 'NR==FNR{ split($0,a,"\t"); lab[a[1]]=a[2]; next }
     { g = ($1 in lab) ? lab[$1] : $3; print $1"\t"$2"\t"g }' \
     "$OUT/relabel.tsv" "$AADR_PREFIX.ind" > "$OUT/panelB.ind"

# --- poplist: ALL pops both models could need, so ONE basis serves both ------
{
  echo "Crete_Lasithi_EMBA"
  echo "Myc_PELOPONNESE"
  echo "Myc_ISLANDS"
  echo "$YAMNAYA"
  echo "$KOTIAS"
  printf '%s\n' "${RIGHT[@]}"
} > "$OUT/poplist.txt"
printf '%s\n' "${RIGHT[@]}" > "$OUT/right.txt"
echo "$YAMNAYA" > "$OUT/yamnaya.txt"
echo "$KOTIAS"  > "$OUT/kotias.txt"

echo ""
echo "[3] poplist (11 pops):"
sed 's/^/    /' "$OUT/poplist.txt"

# --- convertf ---------------------------------------------------------------
cat > "$OUT/par.panelB" <<PAR
genotypename:    $AADR_PREFIX.geno
snpname:         $AADR_PREFIX.snp
indivname:       $OUT/panelB.ind
outputformat:    EIGENSTRAT
genotypeoutname: $OUT/panelB.geno
snpoutname:      $OUT/panelB.snp
indivoutname:    $OUT/panelB.ind.out
poplistname:     $OUT/poplist.txt
hashcheck:       NO
PAR

echo ""
echo "[4] convertf (several minutes)..."
"$CONVERTF" -p "$OUT/par.panelB" > "$OUT/convertf.log" 2>&1 || {
  echo "    STOP: convertf failed. See $OUT/convertf.log"; exit 1; }
mv "$OUT/panelB.ind.out" "$OUT/panelB.ind"

echo "    individuals: $(wc -l < "$OUT/panelB.ind" | tr -d ' ')   panel SNPs: $(wc -l < "$OUT/panelB.snp" | tr -d ' ')"
awk '{c[$3]++} END{for(k in c) printf "      %-32s %d\n", k, c[k]}' "$OUT/panelB.ind" | sort

# every requested pop must be present or dof silently changes
miss=0
while read -r p; do
  [ -n "$p" ] || continue
  awk -v g="$p" '$3==g{f=1} END{exit !f}' "$OUT/panelB.ind" || { echo "    *** MISSING: $p ***"; miss=1; }
done < "$OUT/poplist.txt"
[ "$miss" -eq 0 ] || { echo "    STOP: a poplist population is absent."; exit 1; }

# --- target-genome check ----------------------------------------------------
echo ""
echo "[5] target-genome check:"
grep -qil "mydna\|FAM001\|ID001" "$OUT/panelB.ind" && { echo "    *** FAIL ***"; exit 1; }
echo "    CLEAN — no target genome."

# --- manifest ---------------------------------------------------------------
{
  echo "script: build_paperB_panel"
  echo "date: $(date)"
  echo "purpose: Paper B Test 3 rebuilt on coherent, AADR-labelled populations"
  echo ""
  echo "source:      Crete_Lasithi_EMBA = AADR $SOURCE_GRP (26)"
  echo "  replaces Greece_Minoan_POOL (6 group IDs; rejects rank-0 at p=3.2e-08)"
  echo "  Lasithi is 3,933-4,177 BP: entirely PRE-Mycenaean-conquest."
  echo ""
  echo "descendant:  Greece_Mycenaean_POOL was 5 AADR group IDs. Split:"
  echo "  Myc_PELOPONNESE (18) = ${MYC_MAINLAND[*]}"
  echo "  Myc_ISLANDS     (8)  = ${MYC_ISLAND[*]}  (Aegina, Paros)"
  echo "  These are TESTED for cladality before either is used. A pool is a hypothesis."
  echo ""
  echo "yamnaya: $YAMNAYA (28, one AADR group, Samara Oblast, 647 yr)"
  echo "kotias:  $KOTIAS (2, one site, 1 yr span -- COHERENT, so a Test-3"
  echo "  rejection of it cannot be attributed to incoherence)"
  echo ""
  echo "right_set:"
  printf '  - %s\n' "${RIGHT[@]}"
  echo ""
  echo "panel_snps: $(wc -l < "$OUT/panelB.snp" | tr -d ' ')  (convertf panel, NOT the maxmiss=0 basis)"
  echo "NOTHING here is comparable to 122,891 / 124,362 or any earlier basis."
} > manifest/build_paperB_panel.txt

echo ""
echo "=========================================================="
echo " DONE. Next: Rscript run_test3_clean.R"
echo "=========================================================="
