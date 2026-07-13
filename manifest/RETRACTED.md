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

---

## 2026-07-13 (evening) — Paper B rebuilt; three retractions

### 1. The fold is withdrawn as a reportable quantity.

**Retracted:** the claim that the degradation fold is basis-stable while the
p-value is not, and therefore the quantity to report. Reported folds:
43.8x / 41.9x / 41.9x.

**Why.** Across eleven genuinely different SNP bases the fold spans **0.8x to
80.6x**; the p-value spans 8x. At a *fixed* basis of ~225,000 SNPs, three
replicates differing only in the random SNP draw gave **4.4x / 7.7x / 80.6x**.
A ratio of two noisy p-values is noisier than either.

**Why it looked true.** The three bases cited were 122,891 / 124,362 / ~123,000
SNPs -- three near-copies of one basis. Near-identical bases give near-identical
folds. That is a tautology, and it was presented as the paper's methodological
justification.

**Superseded by** the binary admissibility flip (`run_fold_stability.R`).

### 2. The Test 3 positive control is void.

**Retracted:** `Greece_Mycenaean_POOL <- Greece_Minoan_POOL + Yamnaya`, 41.9x
fold on 122,891 SNPs, "recovers an independently established 86.9/13.1 to three
digits".

**Why.** Both pools are analyst-built bins: `Greece_Minoan_POOL` carries six AADR
group IDs and fails qpWave rank-0 at p=3.2e-08; `Greece_Mycenaean_POOL` carries
five and fails at p=6.3e-05. Rebuilt on coherent populations
(`Myc_PELOPONNESE <- Crete_Lasithi_EMBA + Yamnaya`), the model **rejects at every
basis from 108,442 to 970,005 SNPs**.

**And the relationship was never true.** `Crete_Chania_LBA` -- Late Minoan III,
c. 1350-1225 BCE, *after* the Mycenaean takeover of Knossos -- already carries
steppe ancestry. Pooling it with pre-conquest Lasithi produced a "Minoan" source
that already contained the Mycenaeans, so the model needed less Yamnaya. **The
three-digit agreement with Lazaridis et al. was the circularity, and it looked
exactly like validation.**

**Superseded by** `Germany_Esperstedt_CordedWare <- Yamnaya + Czechia_N_GlobularAmphora`
(Haak et al. 2015), on 898,057 SNPs, all populations AADR groups taken whole.

### 3. Section 3's premise is false.

**Retracted:** "chi^2 scales approximately with SNP count", and the derived rule
"compare chi^2/N".

**Why.** Under RANDOM SNP loss chi^2 does not fall. **It rises** -- 1.36x at 17.7%
loss, in 20 of 20 replicates. Fewer SNPs mean fewer jackknife blocks, a noisier
f4 covariance matrix, an inflated inverse, and an upward-biased chi^2 = E'Q^-1 E.
The chi^2/N ratio has no null and should not be used.

**The exhibit survives, restated.** Against the correct null, the observed chi^2
lies **below all 20 random-loss replicates** at the matched basis, for both thin
populations tested. The diagnostic is not "chi^2 fell faster than N" but
**"chi^2 fell when random loss says it should have risen."**

**Superseded by** `run_chisq_null.R` + `run_chisq_exhibit.R`.

### 4. What was wrong with the whole of v1.0

Every retracted result above was computed on panels that included the author's
own genome and analyst-built pools. **v2.0 is rebuilt entirely from AADR group IDs
taken whole.** The fact that v1.0 was not reproducible from public data is what
allowed most of what it got wrong.
