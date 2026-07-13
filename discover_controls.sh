#!/usr/bin/env bash
# ============================================================================
# discover_controls.sh   -- INSPECT ONLY. Builds nothing, runs no model.
#
# WHY
#   Test 3 requires a descent relationship whose truth is established OUTSIDE
#   the model. Paper B's was  Mycenaean <- Minoan + Yamnaya. It is void:
#     - Greece_Minoan_POOL is not a population (rank-0 p = 3.2e-08)
#     - Greece_Mycenaean_POOL is not a population (rank-0 p = 6.3e-05)
#     - on coherent versions, Myc_PELOPONNESE <- Lasithi + Yamnaya REJECTS at
#       EVERY basis from 108k to 970k
#     - and Myc_PELOPONNESE <- Chania_LBA + Yamnaya FITS (p = 0.58), because
#       Chania is POST-conquest and already carries the steppe the model was
#       supposed to be sourcing. The old "86.9/13.1 matches Lazaridis" was
#       circular.
#
#   So Test 3 has no positive control. This script asks whether the archive
#   contains one.
#
# WHAT MAKES A USABLE CONTROL
#   1. The descent relationship is established in the literature, INDEPENDENTLY
#      of any qpAdm fit we run.
#   2. Both populations are COHERENT -- one AADR group, or a defensible subset.
#      (We are not making the Minoan mistake twice.)
#   3. The descendant is not a plausible ANCESTOR of the source, and does not
#      postdate a takeover by the very ancestry the model is meant to source.
#
# CANDIDATES (each screens ONE of the project's three sources)
#   steppe    Afanasievo <- Yamnaya                 (near-clade; Narasimhan 2019)
#             Corded Ware <- Yamnaya + CE farmer    (Haak 2015, ~75/25)
#   Anatolian Balkan/Aegean Neolithic <- Turkey_N   (European farmers descend
#                                                    from Anatolian farmers)
#   Iran/CHG  Iran Chalcolithic <- Iran_GanjDareh   (weakest; may not exist)
#
# Columns verified: 1=Genetic ID  11=Date mean BP  15=Group ID  16=Locality
# ============================================================================

set -uo pipefail
cd ~/qpadm_project2
ANNO="v66.p1_1240K.aadr.PUB.anno"

echo "=========================================================="
echo " discover_controls.sh   $(date)"
echo "=========================================================="

# --- report: n, localities, date span, publications, platforms --------------
audit () {
  local pat="$1" label="$2"
  echo ""
  echo "=========================================================="
  echo " $label   (group IDs matching: $pat)"
  echo "=========================================================="
  awk -F'\t' -v pat="$pat" '
    FNR==1 { next }
    $15 ~ pat && $15 !~ /_o$|-o$|Ignore/ {
      g=$15; n[g]++
      loc[g,$16]=1
      d=$11+0
      if (!(g in lo) || d<lo[g]) lo[g]=d
      if (d>hi[g]) hi[g]=d
      # platform from the ID suffix
      p="capture"
      if ($1 ~ /\.DG$/) p="shotgun_dip"; else if ($1 ~ /\.SG$/) p="shotgun"
      else if ($1 ~ /\.HO$/) p="HumanOrigins"
      plat[g,p]=1
    }
    END{
      for (g in n) {
        nl=0; for (k in loc) { split(k,a,SUBSEP); if (a[1]==g) nl++ }
        np=0; for (k in plat){ split(k,a,SUBSEP); if (a[1]==g) np++ }
        printf "  %-44s n=%-4d sites=%-3d plat=%-2d  %5d-%5d BP (span %d)\n",
               g, n[g], nl, np, lo[g], hi[g], hi[g]-lo[g]
      }
    }' "$ANNO" | sort -t= -k2 -rn | head -14
}

audit "^Russia_Afanasievo|^Mongolia_.*Afanasievo|Afanasievo" \
      "AFANASIEVO  -- near-clade descendant of Yamnaya (screens STEPPE)"

audit "CordedWare|Corded_Ware" \
      "CORDED WARE -- Yamnaya + CE farmer, ~75/25 (screens STEPPE)"

audit "GlobularAmphora|Globular" \
      "GLOBULAR AMPHORA -- the farmer term in the Corded Ware model"

audit "^Germany_MN|^Germany_LN|^Poland_MN" \
      "CENTRAL EUROPEAN MIDDLE NEOLITHIC -- alternative farmer term"

audit "^Greece_N$|^Greece_Peloponnese_N|^Greece_.*_N$" \
      "AEGEAN NEOLITHIC -- descendant of Turkey_N (screens ANATOLIAN)"

audit "^Bulgaria_N|^Bulgaria_.*_N$|^Serbia_N|^NorthMacedonia_N" \
      "BALKAN NEOLITHIC -- descendant of Turkey_N (screens ANATOLIAN)"

audit "^Hungary_EN|^Hungary_.*Starcevo|Koros|^Germany_LBK|LBK" \
      "CENTRAL EUROPEAN EARLY NEOLITHIC -- descendant of Turkey_N"

audit "^Iran_C$|^Iran_.*_C$|SehGabi|Hajji|^Iran_ChL" \
      "IRAN CHALCOLITHIC -- possible descendant of Iran_GanjDareh (screens IRAN)"

audit "^Russia_Yamnaya|Samara_EBA_Yamnaya|^Russia_Samara_EBA" \
      "YAMNAYA (the source itself, for reference)"

echo ""
echo "=========================================================="
echo " HOW TO READ THIS"
echo "=========================================================="
echo " A usable control needs: n >= 8, ONE locality (or few), a tight date"
echo " span, and one platform. Anything with sites >= 3 or span >= 800 is a"
echo " BIN and is disqualified -- by this project's own Recommendation 2."
echo ""
echo " If nothing usable turns up for Iran_GanjDareh, that hole is PERMANENT"
echo " and Paper B must state that its diagnostic cannot screen the source"
echo " that contributes most to the target. It already half-says this."
echo ""
echo " Nothing is concluded here. Read the counts."
