#!/usr/bin/env bash
# ============================================================================
# discover_test3_pops.sh   -- INSPECT ONLY. Builds nothing, runs no model.
#
# Paper B's headline is:  Mycenaean <- Minoan + Yamnaya  (Test 3 positive control)
#
# Paper A has just shown Greece_Minoan_POOL is not a population (p = 3.2e-08),
# and Paper A's own Recommendation 2 says nothing can descend from such a pool.
# So the SOURCE is void. Before rebuilding, we must know whether the DESCENDANT
# and the other terms are populations either.
#
# This script answers, from AADR's own annotation and nothing else:
#   1. Which individuals are in each analyst-built pool?
#   2. What AADR Group IDs do they actually carry?
#   3. What sites and dates?
#
# It makes no claim. It prints what is there. Read it before we build anything.
# ============================================================================

set -uo pipefail
cd ~/qpadm_project2

ANNO="v66.p1_1240K.aadr.PUB.anno"
IND="postBA_arb.ind"            # the panel Paper B's numbers were computed on
MASTER="v66.p1_1240K.aadr.patch.PUB.ind"

# Columns, verified today: 1=Genetic ID  11=Date mean BP  15=Group ID  16=Locality

POPS=(
  "Greece_Mycenaean_POOL"
  "Greece_Minoan_POOL"
  "Russia_Yamnaya_QC"
  "Georgia_Kotias_QC"
  "Turkey_N_QC"
  "Iran_GanjDareh_N_QC"
  "Russia_Eneolithic_Steppe_QC"
)

echo "=========================================================="
echo " discover_test3_pops.sh   $(date)"
echo "=========================================================="

echo ""
echo "### What population labels actually exist in $IND?"
awk '{c[$3]++} END{for(k in c) if (c[k]>0) printf "  %-34s %d\n", k, c[k]}' "$IND" | sort

for P in "${POPS[@]}"; do
  N=$(awk -v p="$P" '$3==p' "$IND" | wc -l | tr -d ' ')
  echo ""
  echo "=========================================================="
  echo " $P   (n = $N in $IND)"
  echo "=========================================================="
  if [ "$N" -eq 0 ]; then
    echo "  NOT PRESENT under this label. Skipping."
    continue
  fi

  awk -v p="$P" '$3==p{print $1}' "$IND" | sort > /tmp/pop_ids.txt

  echo ""
  echo "  --- AADR Group IDs carried by these individuals ---"
  awk -F'\t' 'NR==FNR{w[$1]=1;next} FNR==1{next} ($1 in w){c[$15]++}
    END{for(k in c) printf "    %-42s %d\n", k, c[k]}' \
    /tmp/pop_ids.txt "$ANNO" | sort -rn -k2

  NGRP=$(awk -F'\t' 'NR==FNR{w[$1]=1;next} FNR==1{next} ($1 in w){c[$15]=1}
    END{print length(c)}' /tmp/pop_ids.txt "$ANNO")
  echo "    >>> $NGRP distinct AADR Group ID(s)"

  echo ""
  echo "  --- Locality x date range ---"
  awk -F'\t' 'NR==FNR{w[$1]=1;next} FNR==1{next} ($1 in w){
      k=$16; c[k]++
      if(!(k in lo) || $11+0 < lo[k]) lo[k]=$11+0
      if($11+0 > hi[k]) hi[k]=$11+0 }
    END{for(k in c) printf "    %3d  %6d - %6d BP   %s\n", c[k], lo[k], hi[k], k}' \
    /tmp/pop_ids.txt "$ANNO" | sort -rn

  echo ""
  echo "  --- overall date span ---"
  awk -F'\t' 'NR==FNR{w[$1]=1;next} FNR==1{next} ($1 in w){
      d=$11+0; if(mn==0||d<mn) mn=d; if(d>mx) mx=d }
    END{printf "    %d - %d BP   (span %d years)\n", mn, mx, mx-mn}' \
    /tmp/pop_ids.txt "$ANNO"
done

echo ""
echo "=========================================================="
echo " THE QUESTION THIS ANSWERS"
echo "=========================================================="
echo " If Greece_Mycenaean_POOL carries >1 AADR Group ID, then Paper B's"
echo " positive control has a bin on BOTH sides -- source and descendant --"
echo " and the 86.9/13.1 recovery it cites as proof of truth was computed"
echo " between two populations that are not populations."
echo ""
echo " Nothing is concluded here. Read the counts."
