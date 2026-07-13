#!/usr/bin/env bash
# ============================================================================
# build_paperB_wide.sh
#
# WHY
#   Myc_PELOPONNESE <- Crete_Lasithi_EMBA + Yamnaya REJECTS at p = 8.6e-07 on
#   970,660 SNPs. Two explanations, and they have opposite consequences:
#
#     (a) POWER. At ~1M SNPs qpAdm rejects models that are approximately true.
#         If NOTHING fits Myc_PELOPONNESE at this basis, the rejection says
#         nothing about the descent relationship. Test 3 survives; it just
#         needs a stated basis regime.
#
#     (b) MISFIT. If some OTHER two-source model DOES fit while Lasithi+Yamnaya
#         does not, the descent relationship is wrong, and Test 3 has no
#         positive control. Paper B loses its exhibit.
#
#   This panel lets us ask (b). run_myc_power.R asks (a).
#
# SOURCES ADDED (AADR group IDs, coherence declared for each)
#   Crete_Chania_LBA      Greece_Crete_LBA minus kro008/kro009 (Heraklion).
#                         POST-conquest Cretans (~1350-1225 BCE). If Mycenaeans
#                         fit BETTER from these than from pre-conquest Lasithi,
#                         that is a finding, not a nuisance.
#   Turkey_N              68 in AADR; a 5-of-5-flag BIN by Paper A's own audit.
#                         Included as a competitor only, and declared as a bin.
#   Iran_GanjDareh_N      one site, coherent.
#   Khvalynsk_Eneolithic  Russia_Saratov_Eneolithic_Khvalynsk. One site cluster.
#                         (The project's Russia_Eneolithic_Steppe_QC pooled this
#                         with Ekaterinovka AND an AADR '-o' OUTLIER. Not used.)
#
#   All restricted to individuals retained by this project's QC where the pool
#   existed; taken whole from AADR otherwise. Declared either way.
# ============================================================================

set -euo pipefail
cd ~/qpadm_project2

AADR_PREFIX="/Users/ell_thales/qpadm_project2/v66.p1_1240K.aadr.patch.PUB"
ANNO="/Users/ell_thales/qpadm_project2/v66.p1_1240K.aadr.PUB.anno"
CONVERTF="/Users/ell_thales/miniconda3/pkgs/admixtools-8.0.2-ha8dfbab_0/bin/convertf"

RIGHT=("Mbuti" "Italy_Sicily_Epigravettian" "Russia_Kostenki_UP" "Papuan" \
       "Serbia_IronGates_Mesolithic" "Morocco_Iberomaurusian")

OUT="clean_paperB_wide"
mkdir -p "$OUT" manifest

echo "=========================================================="
echo " build_paperB_wide.sh   $(date)"
echo "=========================================================="

COL_GRP=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if ($i ~ /^Group ID$/) {print i; exit}}' "$ANNO")
[ "$COL_GRP" = "15" ] || { echo "STOP: Group ID column moved"; exit 1; }

# individuals retained by this project's QC (for the pools that had one)
awk '$3=="Greece_Minoan_POOL" || $3=="Greece_Mycenaean_POOL" {print $1}' \
  postBA_arb.ind | sort > "$OUT/qc_kept.txt"

awk -F'\t' -v cg="$COL_GRP" '
  NR==FNR { keep[$1]=1; next }
  FNR==1  { next }
  {
    id=$1; g=$cg
    # --- QC-restricted (these pools had an upstream QC step) ---
    if (id in keep) {
      if      (g=="Greece_Crete_HgCharalambos_EMBA")            print id "\tCrete_Lasithi_EMBA"
      else if (g=="Greece_Crete_LBA" && id!="kro008.SG" && id!="kro009.SG") print id "\tCrete_Chania_LBA"
      else if (g=="Greece_Peloponnese_LBA" || g=="Greece_PalaceofNestor_BA" \
            || g=="Greece_Argolid_LBA"     || g=="Greece_Achaea_LBA")       print id "\tMyc_PELOPONNESE"
    }
    # --- taken whole from AADR (no project pool existed, or it was a bin) ---
    if      (g=="Russia_Samara_EBA_Yamnaya")             print id "\tYamnaya_Samara"
    else if (g=="Georgia_KotiasKlde_Mesolithic")         print id "\tKotias"
    else if (g=="Turkey_N")                              print id "\tTurkey_N_BIN"
    else if (g=="Iran_GanjDareh_N")                      print id "\tIran_GanjDareh"
    else if (g=="Russia_Saratov_Eneolithic_Khvalynsk")   print id "\tKhvalynsk_Eneo"
  }
' "$OUT/qc_kept.txt" "$ANNO" | sort -u > "$OUT/relabel.tsv"

echo ""
echo "[1] populations:"
awk -F'\t' '{c[$2]++} END{for(k in c) printf "    %-24s %d\n", k, c[k]}' "$OUT/relabel.tsv" | sort

for pair in "Crete_Lasithi_EMBA:26" "Crete_Chania_LBA:23" "Myc_PELOPONNESE:18"; do
  p="${pair%%:*}"; want="${pair##*:}"
  got=$(awk -F'\t' -v p="$p" '$2==p' "$OUT/relabel.tsv" | wc -l | tr -d ' ')
  [ "$got" -eq "$want" ] || { echo "    STOP: $p expected $want, got $got"; exit 1; }
done
echo "    QC-restricted counts match. good."

awk 'NR==FNR{ split($0,a,"\t"); lab[a[1]]=a[2]; next }
     { g = ($1 in lab) ? lab[$1] : $3; print $1"\t"$2"\t"g }' \
     "$OUT/relabel.tsv" "$AADR_PREFIX.ind" > "$OUT/wide.ind"

{
  echo "Myc_PELOPONNESE"; echo "Crete_Lasithi_EMBA"; echo "Crete_Chania_LBA"
  echo "Yamnaya_Samara";  echo "Khvalynsk_Eneo";     echo "Turkey_N_BIN"
  echo "Iran_GanjDareh";  echo "Kotias"
  printf '%s\n' "${RIGHT[@]}"
} > "$OUT/poplist.txt"
printf '%s\n' "${RIGHT[@]}" > "$OUT/right.txt"

echo ""
echo "[2] poplist ($(wc -l < "$OUT/poplist.txt" | tr -d ' ') pops)"

cat > "$OUT/par.wide" <<PAR
genotypename:    $AADR_PREFIX.geno
snpname:         $AADR_PREFIX.snp
indivname:       $OUT/wide.ind
outputformat:    EIGENSTRAT
genotypeoutname: $OUT/wide.geno
snpoutname:      $OUT/wide.snp
indivoutname:    $OUT/wide.ind.out
poplistname:     $OUT/poplist.txt
hashcheck:       NO
PAR

echo ""
echo "[3] convertf (10-20 min)..."
"$CONVERTF" -p "$OUT/par.wide" > "$OUT/convertf.log" 2>&1 || {
  echo "    STOP: convertf failed. See $OUT/convertf.log"; exit 1; }
mv "$OUT/wide.ind.out" "$OUT/wide.ind"

echo "    individuals: $(wc -l < "$OUT/wide.ind" | tr -d ' ')   panel SNPs: $(wc -l < "$OUT/wide.snp" | tr -d ' ')"
awk '{c[$3]++} END{for(k in c) printf "      %-26s %d\n", k, c[k]}' "$OUT/wide.ind" | sort

# target check -- by POPULATION LABEL, not by substring.
# (The substring check false-positived on AID001.AG, an Aidonia individual.)
echo ""
echo "[4] target-genome check (by population label):"
N=$(awk '$3=="MyDNA"' "$OUT/wide.ind" | wc -l | tr -d ' ')
[ "$N" -eq 0 ] && echo "    CLEAN — MyDNA absent." || { echo "    *** FAIL: $N ***"; exit 1; }

{
  echo "script: build_paperB_wide"
  echo "date: $(date)"
  echo "purpose: can ANY two-source model fit Myc_PELOPONNESE at ~1M SNPs?"
  echo ""
  echo "QC-restricted (project QC retained these; the QC is NOT scripted):"
  echo "  Crete_Lasithi_EMBA (26 of AADR's 34)"
  echo "  Crete_Chania_LBA   (23; Greece_Crete_LBA minus kro008/kro009 = Heraklion)"
  echo "  Myc_PELOPONNESE    (18 of AADR's 31)"
  echo "Taken whole from AADR:"
  echo "  Yamnaya_Samara (Russia_Samara_EBA_Yamnaya)"
  echo "  Kotias (Georgia_KotiasKlde_Mesolithic)"
  echo "  Turkey_N_BIN (Turkey_N -- 5-of-5 coherence flags. A BIN. Competitor only.)"
  echo "  Iran_GanjDareh (Iran_GanjDareh_N)"
  echo "  Khvalynsk_Eneo (Russia_Saratov_Eneolithic_Khvalynsk)"
  echo ""
  echo "right_set:"
  printf '  - %s\n' "${RIGHT[@]}"
  echo ""
  echo "panel_snps: $(wc -l < "$OUT/wide.snp" | tr -d ' ')  (convertf panel, NOT the basis)"
} > manifest/build_paperB_wide.txt

echo ""
echo " DONE. Next: Rscript run_myc_power.R"
