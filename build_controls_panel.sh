#!/usr/bin/env bash
# ============================================================================
# build_controls_panel.sh
#
# Paper B's positive control is dead. This builds the candidates for a new one.
#
# WHAT DIED
#   Mycenaean <- Minoan + Yamnaya. Both pools are bins (rank-0 p = 3.2e-08 and
#   6.3e-05). On coherent versions the model REJECTS at every basis from 108k
#   to 970k -- while Myc <- Chania_LBA + Yamnaya FITS (p = 0.58), because
#   Chania is POST-conquest and already carries the steppe the model was
#   supposed to source. The "86.9/13.1 matches Lazaridis" was circular.
#
# THE NEW CONTROLS -- chosen from AADR by COHERENCE, not by convenience
#   Every one was selected on n / sites / span / platforms BEFORE any model ran.
#
#   1. STEPPE, two-source (canonical; Haak et al. 2015, ~75/25)
#        Germany_Esperstedt_CordedWare (13, ONE site, 247 yr)
#          <- Russia_Samara_EBA_Yamnaya + Czechia_N_GlobularAmphora
#
#   2. STEPPE, one-source (near-clade; Narasimhan et al. 2019)
#        Russia_Khakassia_Afanasievo (27, 2 sites, 400 yr)
#          <- Russia_Samara_EBA_Yamnaya
#
#   3. ANATOLIAN (European farmers descend from Anatolian farmers)
#        Slovakia_N_LBK (56, 2 sites, 164 yr)          <- Turkey_N
#        Germany_HalberstadtSonntagsfeld_EN_LBK (26, ONE site, 187 yr) <- Turkey_N
#      >>> This closes Paper B's largest admitted hole: it currently CANNOT
#      >>> screen Turkey_N, the source contributing most to the target.
#
# DELIBERATELY REJECTED, and why (all disqualified by this project's own
# Recommendation 2 -- a pool is a hypothesis, and a bin is not a population)
#   Austria_N_LBK               103 individuals, ONE site -- but 1,453 yr span.
#                               The largest sample on the list. Disqualified.
#   Czechia_EBA_CordedWare      48, but 15 sites / 706 yr. The obvious choice.
#   Poland_GlobularAmphora      33, but 7 sites / 669 yr.
#   Russia_Altai_Afanasievo     20, but 7 sites.
#
# A NOTE WE DO NOT HIDE
#   Russia_Samara_EBA_Yamnaya is itself 15 localities / 3 platforms / 647 yr.
#   It would trip our own ">= 3 localities" flag. It is the steppe source in
#   BOTH papers. We use it because there is no coherent alternative with n>=20,
#   and we say so.
#
#   IRAN: there is NO usable descendant of Iran_GanjDareh in v66.p1.
#   Iran_SehGabi_C is n=7 with an 871-yr span; Iran_TepeHissar_C is 1,523 yr
#   across 3 platforms. The Iran/CHG source CANNOT be Test-3 screened.
#   That hole is STRUCTURAL AND PERMANENT and must be stated as such.
# ============================================================================

set -euo pipefail
cd ~/qpadm_project2

AADR_PREFIX="/Users/ell_thales/qpadm_project2/v66.p1_1240K.aadr.patch.PUB"
ANNO="/Users/ell_thales/qpadm_project2/v66.p1_1240K.aadr.PUB.anno"
CONVERTF="/Users/ell_thales/miniconda3/pkgs/admixtools-8.0.2-ha8dfbab_0/bin/convertf"

RIGHT=("Mbuti" "Italy_Sicily_Epigravettian" "Russia_Kostenki_UP" "Papuan" \
       "Serbia_IronGates_Mesolithic" "Morocco_Iberomaurusian")

# targets + sources + the candidate reference, all under AADR's own names
POPS=(
  "Germany_Esperstedt_CordedWare"            # target 1   (13, 1 site)
  "Russia_Khakassia_Afanasievo"              # target 2   (27, 2 sites)
  "Slovakia_N_LBK"                           # target 3   (56, 2 sites)
  "Germany_HalberstadtSonntagsfeld_EN_LBK"   # target 3b  (26, 1 site)
  "Russia_Samara_EBA_Yamnaya"                # steppe source
  "Czechia_N_GlobularAmphora"                # farmer term for Corded Ware
  "Turkey_N"                                 # Anatolian source (a 5-flag BIN; declared)
  "Georgia_KotiasKlde_Mesolithic"            # THE CANDIDATE REFERENCE
)

# expected counts from discover_controls.sh -- the script dies if AADR disagrees
declare -a EXPECT=( 13 27 56 26 46 6 68 2 )

OUT="clean_controls"
mkdir -p "$OUT" manifest

echo "=========================================================="
echo " build_controls_panel.sh   $(date)"
echo "=========================================================="

COL_GRP=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if ($i ~ /^Group ID$/) {print i; exit}}' "$ANNO")
[ "$COL_GRP" = "15" ] || { echo "STOP: Group ID column moved"; exit 1; }
echo "[1] Group ID column = 15 (verified)"

echo ""
echo "[2] population counts (from AADR, taken WHOLE — no private QC filter):"
i=0
for P in "${POPS[@]}"; do
  n=$(awk -F'\t' -v cg="$COL_GRP" -v g="$P" 'FNR>1 && $cg==g' "$ANNO" | wc -l | tr -d ' ')
  want="${EXPECT[$i]}"
  printf "    %-42s n=%-4s (expected %s)" "$P" "$n" "$want"
  if [ "$n" -ne "$want" ]; then echo "   *** MISMATCH ***"; exit 1; fi
  echo "   ok"
  i=$((i+1))
done
echo "    All counts match discovery. No unscripted QC filter is applied here:"
echo "    these are AADR's groups, whole. That makes them PUBLICLY reproducible."

{ printf '%s\n' "${POPS[@]}"; printf '%s\n' "${RIGHT[@]}"; } > "$OUT/poplist.txt"
printf '%s\n' "${RIGHT[@]}" > "$OUT/right.txt"

echo ""
echo "[3] poplist: $(wc -l < "$OUT/poplist.txt" | tr -d ' ') populations"

cat > "$OUT/par.controls" <<PAR
genotypename:    $AADR_PREFIX.geno
snpname:         $AADR_PREFIX.snp
indivname:       $AADR_PREFIX.ind
outputformat:    EIGENSTRAT
genotypeoutname: $OUT/controls.geno
snpoutname:      $OUT/controls.snp
indivoutname:    $OUT/controls.ind
poplistname:     $OUT/poplist.txt
hashcheck:       NO
PAR

echo ""
echo "[4] convertf (10-20 min)..."
"$CONVERTF" -p "$OUT/par.controls" > "$OUT/convertf.log" 2>&1 || {
  echo "    STOP: convertf failed. See $OUT/convertf.log"; exit 1; }

echo "    individuals: $(wc -l < "$OUT/controls.ind" | tr -d ' ')   panel SNPs: $(wc -l < "$OUT/controls.snp" | tr -d ' ')"
awk '{c[$3]++} END{for(k in c) printf "      %-42s %d\n", k, c[k]}' "$OUT/controls.ind" | sort

miss=0
while read -r p; do
  [ -n "$p" ] || continue
  awk -v g="$p" '$3==g{f=1} END{exit !f}' "$OUT/controls.ind" || { echo "    *** MISSING: $p ***"; miss=1; }
done < "$OUT/poplist.txt"
[ "$miss" -eq 0 ] || { echo "    STOP: a poplist population is absent."; exit 1; }

echo ""
echo "[5] target-genome check (by population label, not substring):"
N=$(awk '$3=="MyDNA"' "$OUT/controls.ind" | wc -l | tr -d ' ')
[ "$N" -eq 0 ] && echo "    CLEAN — MyDNA absent." || { echo "    *** FAIL ***"; exit 1; }

{
  echo "script: build_controls_panel"
  echo "date: $(date)"
  echo "purpose: three literature-established positive controls for Test 3"
  echo ""
  echo "NO PRIVATE FILTER. Every population is an AADR group taken whole."
  echo "  This panel is reproducible from the public archive alone."
  echo ""
  echo "controls:"
  echo "  1. Germany_Esperstedt_CordedWare (13, 1 site, 247 yr)"
  echo "       <- Russia_Samara_EBA_Yamnaya + Czechia_N_GlobularAmphora"
  echo "       Haak et al. 2015; canonical ~75/25. Screens STEPPE."
  echo "  2. Russia_Khakassia_Afanasievo (27, 2 sites, 400 yr)"
  echo "       <- Russia_Samara_EBA_Yamnaya   [one-source; near-clade]"
  echo "       Narasimhan et al. 2019. Screens STEPPE."
  echo "  3. Slovakia_N_LBK (56, 2 sites, 164 yr) <- Turkey_N"
  echo "     Germany_HalberstadtSonntagsfeld_EN_LBK (26, 1 site) <- Turkey_N"
  echo "       Screens ANATOLIAN -- the source Paper B could not screen."
  echo ""
  echo "candidate reference under test: Georgia_KotiasKlde_Mesolithic (2, 1 site)"
  echo ""
  echo "disqualified as bins (by this project's own Recommendation 2):"
  echo "  Austria_N_LBK (103 ind, 1 site, but 1,453 yr span)"
  echo "  Czechia_EBA_CordedWare (48, 15 sites, 706 yr)"
  echo "  Poland_GlobularAmphora (33, 7 sites, 669 yr)"
  echo "  Russia_Altai_Afanasievo (20, 7 sites)"
  echo ""
  echo "DECLARED: Russia_Samara_EBA_Yamnaya is itself 15 localities / 3 platforms"
  echo "  / 647 yr and would trip our own >=3-localities flag. No coherent"
  echo "  alternative with n>=20 exists. Used, and declared."
  echo ""
  echo "DECLARED: NO usable descendant of Iran_GanjDareh exists in v66.p1."
  echo "  Iran_SehGabi_C is n=7 / 871 yr; Iran_TepeHissar_C is 1,523 yr / 3 plat."
  echo "  The Iran/CHG source CANNOT be Test-3 screened. Structural, permanent."
  echo ""
  echo "right_set:"
  printf '  - %s\n' "${RIGHT[@]}"
  echo ""
  echo "panel_snps: $(wc -l < "$OUT/controls.snp" | tr -d ' ')  (convertf panel, NOT the basis)"
} > manifest/build_controls_panel.txt

echo ""
echo " DONE. Next: Rscript run_test3_controls.R"
