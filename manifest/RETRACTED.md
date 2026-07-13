# Retracted runs, retained deliberately

These manifests come from analyses that were **withdrawn**. They are kept
because the retraction sequence is part of the papers' argument (see the
prediction ledger), not because their numbers are usable.

- `run_test3_sourceswap_*` — a 2x2 factorial. **VOID.** All four baselines were
  inadmissible: every model of Turkey_MLBA carrying an obligatory steppe source
  returned a significantly negative weight, because Turkey_MLBA has no steppe
  ancestry. Superseded by `run_test3_controls_v3.R`.

- `run_test3_controls_*` (no v3) — same defect. Superseded by v3.

- `run_eastern_term_*` (no v2) — **BASIS-POLLUTED.** Forced all four third-source
  candidates into one extract; Turkey_Epipaleolithic (n=1) cost 17.7% of the SNP
  basis and flipped the two-source rejections from p=0.004 to p=0.11. Superseded
  by `run_eastern_term_v2.R`, which is the source of the chi-squared exhibit
  built *from* this failure.

No figure in either paper is taken from these files.

- `scripts/run_test3.R` — the FIRST version of the Test 3 script. Superseded by
  `run_test3_v2.R`. Retained for the record; no figure in either paper comes
  from it. The v2 script is canonical.

---

## 2026-07-13 — §2.5 rebuilt; two retractions

### 1. The §2.5 p-values are withdrawn.

**Retracted:** rank-0 p = 0.0074 (Lasithi vs Chania, 88,776 SNPs) and
p = 0.0014 (Early vs Late, 121,732 SNPs), from `run_anchor_integrity.R`.

**Why.** Both tests were run on `postBA_arb`, a panel into which the study's
target genome had been merged. The SNP set was therefore intersected down to a
23andMe v5 chip **before any extract ran** — so a qpWave test between two
*ancient* sub-pools, neither of which contains the target, was nonetheless
conditioned on the target's chip. Every basis produced by that script sits under
the ~123,871 SNP ceiling the chip imposes. That is not a coincidence; it is the
panel.

**A second defect, compounding.** The site arm's extract also carried
`Minoan_OTHER` — nine individuals its own model never used — at a cost of
~33,000 SNPs. Every population in an `extract_f2` call sets the basis, not only
those named in the model. This project wrote that rule and broke it here.

**A third.** Sub-pools were assigned by substring-matching the AADR Locality
field (`/lasithi|chania/`). AADR already gives these 58 individuals **six
distinct site-level Group IDs**. The archive had done the right thing and the
heuristic was unnecessary.

**Superseded by** `build_minoan_panel.sh` → `run_minoan_split_clean.R`, which
build the panels ancients-only from AADR v66.p1 with no target genome anywhere,
and define sub-pools from AADR's own Group IDs:

    Lasithi vs Chania   chisq 43.32  dof 5  p = 3.2e-08   950,383 SNPs
    Early vs Late       chisq 52.55  dof 5  p = 4.2e-10   955,096 SNPs

The direction is unchanged and the evidence is far stronger. The retracted
numbers are named in the paper rather than quietly replaced.

### 2. "The fault line is temporal" is withdrawn.

**Retracted:** the claim in §2.5 that the pool's structure is temporal rather
than spatial, inferred from the date arm rejecting harder than the site arm.

**Why.** The two arms are not independent. The date arm assigns all 49 site-arm
individuals exactly as the site arm does, and adds nine others; its higher
chi-squared reflects the extra individuals, not a sharper fault. **Site and date
are perfectly confounded in this pool**: there are no individuals between 3,525
and 3,768 BP, so Lasithi *is* the early cluster and Chania *is* the late one. No
design on this sample can separate them.

**Also withdrawn:** "the site split rejected while underpowered." It did not.
It was starved of ~860,000 SNPs by a chip-conditioned panel and a surplus
population, both of which were our doing.

**Not recomputed:** `Minoan_LATE <- Minoan_EARLY + Yamnaya` (8.6% Yamnaya,
Z = 1.54, n.s.) was run on the retracted basis and has not been redone on the
clean panel. It is marked [OWED] in the paper. The pool is heterogeneous; the
cause is not established.

### 3. `basis_census.R` did not exist.

The census establishing that this project had run across many distinct SNP bases
was executed interactively and **never written to disk**. A project whose thesis
is that unexamined provenance determines results ran its own provenance audit
off-script. `basis_census.R` now exists and is versioned. Its output is a
candidate list with false positives, not an audited count, and no figure in
either paper is taken from it.
