# The SNP basis is not a nuisance parameter

### A positive-control screen for reference-into-source gene flow, the SNP count it requires, and a χ² that rises when SNPs are removed at random

**PREPRINT DRAFT v2.0. Not peer reviewed.**

**Alexandros Karalakis** · Independent Researcher · alexandros.karalakis@me.com

---

> ## What v2.0 retracts
>
> **This version withdraws the central methodological claim of v1.0, its headline exhibit, and the premise of its second diagnostic. All three were wrong. They are named here rather than quietly replaced.**
>
> **1. The fold is not a stable statistic, and v1.0's demonstration that it was is a tautology.** v1.0 reported the degradation *fold* — the ratio p(references) / p(references + candidate) — and argued the fold rather than the p-value should be reported **because the fold was stable across SNP bases**. Measured properly, across eleven bases, **the fold spans 0.8× to 80.6× while the p-value spans 8×.** The fold is *twelve times less stable* than the quantity it was meant to replace. v1.0's "three bases" were 122,891 / 124,362 / ~123,000 SNPs — three near-copies of one basis. **Three near-identical bases produce three near-identical folds. That is not evidence of stability.**
>
> **2. The Test 3 exhibit was built on two bins and a circular relationship.** `Greece_Mycenaean_POOL ← Greece_Minoan_POOL + Yamnaya` is void: **both pools fail qpWave rank-0** (p = 6.3 × 10⁻⁵ and 3.2 × 10⁻⁸), and on coherent populations the relationship **rejects at every SNP basis from 108,442 to 970,005.** The apparent fit came from post-conquest Cretans in the "Minoan" source who already carried steppe ancestry. **The 86.9 / 13.1 that matched the literature to three digits was circularity, and it looked exactly like validation.**
>
> **3. §3's premise — "χ² scales approximately with SNP count" — is false.** Under **random** SNP loss, χ² does not fall. **It rises** (1.36× at 17.7% loss; 20 of 20 replicates). The exhibit survives, but **for the opposite reason to the one v1.0 gave.**
>
> **Every result in this version is rebuilt from AADR v66.p1 alone.** No analyst-built pool, no private genome, no unscripted filter. The panels, scripts and logs are in the repository.

---

## Abstract

qpAdm results are reported with a p-value, a set of weights, and a reference set. **They are rarely reported with the number of SNPs the statistics were computed on, and almost never with the list of populations that determined that number.** This paper is about what that omission costs.

**We propose Test 3, a positive-control screen for reference-into-source gene flow.** Harney et al. (2021) showed by simulation that a reference population which has contributed ancestry to a *source* biases qpAdm's admixture estimates; Flegontova et al. (2025) specify the violation formally (their Fig. 1a). **Neither offers a way to detect it in real data, where the truth is unknown.** Test 3 does: fit a population *independently known* to descend from a source, with and without the candidate reference. **If the candidate breaks a relationship established outside the model, the model cannot be what is wrong.**

Test 3 fires on a canonical control. `Germany_Esperstedt_CordedWare ← Yamnaya + Globular Amphora` (Haak et al. 2015) **holds** with six references (p = 0.267) and **breaks** when a CHG reference is added (p = 0.016) — while **all six legitimate references leave it intact.** The mechanism is documented: CHG is an ancestral component of Yamnaya, so the reference is an ancestor of the source.

**But Test 3 requires roughly 450,000 SNPs, and below that it does not merely lose power — it fails in two opposite directions.** At ~225,000 SNPs an *innocent* reference produces a larger effect than the guilty one does at full basis. At ~108,000 SNPs — the basis at which much 1240K analysis is actually done — **the guilty reference becomes invisible.** False positives above, false negatives below. **Running the screen below its basis requirement is worse than not running it.**

**Separately, we report a property of χ² that we have not seen stated.** Under **random** SNP loss, qpAdm's χ² **rises**: fewer SNPs mean fewer jackknife blocks, a noisier f4 covariance matrix, an inflated inverse, and an upward-biased χ² = E′Q⁻¹E. **It rose in 20 of 20 replicates at 17.7% loss.** This makes a *falling* χ² diagnostic. When a low-coverage genome is added to an `extract_f2` call at `maxmiss = 0` — **a population the model never uses** — the SNP basis collapses to that genome's coverage footprint, and χ² falls **below every one of twenty random-loss replicates at the same basis.** In one case the model's conclusion reverses: **p = 0.016 → 0.069, from a population that appears nowhere in the model.**

**We do not claim the screen is validated. One control fires.** We do not claim novelty; the citation search is owed and named as owed. What we do claim is that **the SNP basis is a free parameter that most qpAdm analyses neither report nor control, and that both of the phenomena above are invisible without it.**

---

## 1. The problem

A qpAdm model is a target, some sources, a reference set, and a p-value. Three of the four are reported. **The fourth thing — the set of SNPs on which every f4-statistic was computed — is set by a population list that need not appear anywhere in the model, and is almost never published.**

At `maxmiss = 0`, `extract_f2` retains only sites genotyped in **every population in the extract call.** One low-coverage genome in that call — not in the model, not a source, not a reference, merely present — and the basis becomes that genome's coverage footprint. **Everyone else's SNPs are deleted to match it.**

This paper reports two things that follow.

---

## 2. Test 3: a positive-control screen for reference-into-source gene flow

### 2.1 What the field already has, and what it does not

**The violation is specified. The remedy is not.**

Harney et al. (2021) showed by simulation that gene flow from a reference population into a **source** — after that source's divergence from the true admixing lineage — biases qpAdm's admixture estimates, and that the bias is worse when the gene-flow donor sits in the reference set. Flegontova et al. (2025) enumerated the assumption violations that cause a **true** model to be rejected; **their Figure 1(a) is exactly the case treated here.** The violation is not our discovery. It is named, published and formally specified.

**What neither provides is a way to detect it in real data.**

- **Harney et al. (2021) is a simulation study.** Its guidance on the reference set concerns **quantity**: qpAdm begins rejecting otherwise-plausible models once roughly thirty additional populations are added. That is a warning about how *many* references to use, not a test of whether a *particular* one is admissible.
- **Flegontova et al. (2025) recommends protocol-level hygiene**: temporal stratification, testing fewer models to raise pre-study odds, stricter conditions on estimated admixture fractions. **These reduce a screen's false-discovery rate. None interrogates a named reference population.**

### 2.2 The closest existing thing, and why it is not the same thing

**qpAdm already reports which reference is driving a misfit.** The original ADMIXTOOLS implementation prints a `dscore` block identifying which outgroup causes the model to deviate from the data when the fit is poor. Skoglund's `qpAdm_wrapper` surfaces the same information, reporting the reference with the largest |Z| for f4(Target, Fitted_target; Base, Reference) as a hint for improving the reference list. Petr's `admixr` searches reference subsets to find **redundant** references.

**These are attribution and parsimony tools. They are not validity tests, and the distinction is not cosmetic.**

`dscore` operates **inside a model the analyst is already fitting, and which is already failing.** It answers: *given that this model misfits, which reference contributes most?* It therefore **cannot separate "this reference is bad" from "this model is wrong"** — and in the situation where `dscore` is consulted, the second is usually the live hypothesis.

**Test 3 removes the confound by construction.** It runs the candidate against a relationship whose truth is established **outside the model**. **When that relationship breaks, the model cannot be the thing that is wrong.**

| | asks | confounded with model error? | needs external truth? |
|---|---|---|---|
| `dscore` / wrapper | *which reference drives this misfit?* | **yes** | no |
| `admixr` subset search | *which references are redundant?* | n/a | no |
| **Test 3** | *does this reference break a relationship known to be true?* | **no** | **yes** |

> **The contribution is not a new statistic. It is a positive control. `dscore` is a pointer; Test 3 is a test.**

### 2.3 The control

**`Germany_Esperstedt_CordedWare ← Russia_Samara_EBA_Yamnaya + Czechia_N_GlobularAmphora`**

Chosen on three properties, all fixed before any model was run:

- **Canonical.** Haak et al. (2015): Corded Ware as a ~75/25 mixture of Yamnaya-related steppe and central European farmer ancestry. Undisputed for a decade.
- **Coherent.** `Germany_Esperstedt_CordedWare`: 13 individuals, **one site**, 247-year span. Selected on the archive's own annotation, before it was fitted. *(The obvious alternative, `Czechia_EBA_CordedWare`, has 48 individuals across 15 sites and 706 years. It is a bin, and it was rejected on those grounds, not on its results.)*
- **Not ours.** All three are AADR group IDs taken whole.

**Reference set (six):** `Mbuti`, `Italy_Sicily_Epigravettian`, `Russia_Kostenki_UP`, `Papuan`, `Serbia_IronGates_Mesolithic`, `Morocco_Iberomaurusian`.

**Candidate:** `Georgia_KotiasKlde_Mesolithic` (CHG; 2 individuals, one site, **1-year date span** — coherent, so a rejection cannot be attributed to incoherence).

### 2.4 The result: an admissibility flip

`maxmiss = 0`; **898,057 SNPs**. Kotias is in the `extract_f2` call for **both** fits, so the basis is identical and **only the reference set changes.**

| Reference set | χ² | dof | p | Yamnaya | Globular Amphora | |
|---|---|---|---|---|---|---|
| **six canonical** | 5.21 | 4 | **0.267** | 64.7% (Z 15.05) | 35.3% (Z 8.20) | **HOLDS** |
| **six + `Georgia_Kotias`** | 13.93 | 5 | **0.016** | 71.6% (Z 20.08) | 28.4% (Z 7.98) | **BROKEN** |

**And the six legitimate references do not break it.** Leave-one-out and restore, same basis:

| Reference | effect | |
|---|---|---|
| `Serbia_IronGates_Mesolithic` | 0.62× | passes |
| `Papuan` | 0.62× | passes |
| `Russia_Kostenki_UP` | 0.67× | passes |
| `Italy_Sicily_Epigravettian` | 0.99× | passes |
| `Mbuti` | 2.18× | passes |
| `Morocco_Iberomaurusian` | 2.82× | passes |
| **`Georgia_Kotias`** | **16.6×** | **REJECTED** |

> **The screen discriminates. It is not merely sensitive to perturbing the reference set — it is sensitive to the *right* perturbation.**

*(The leave-one-out compares 5→6 references while the Kotias test compares 6→7. The asymmetry is unavoidable and is declared. §2.6 explains why the ratio is not the reportable quantity in any case.)*

### 2.5 The mechanism

**CHG is an ancestral component of Yamnaya** (Haak et al. 2015; Wang et al. 2019). **The reference is an ancestor of the source** — the Harney violation as specified, and Flegontova et al. (2025) Fig. 1(a), occurring in real data with a documented history.

**The weight moves in the predicted direction.** Adding Kotias pushes Yamnaya **up**, 64.7% → 71.6%. A reference sharing drift with a source distorts the estimate toward it.

**v1.0 said "the mechanism is not established." On this control it is** — named, documented and directionally consistent, though not isolated experimentally.

### 2.6 The SNP basis requirement, and two opposite failure modes

**This is the most important section in the paper, and it retracts v1.0's central methodological claim.**

Same model, same populations, same reference set. **Only the SNP basis changes** — random subsamples of the panel, three replicates per size, plus a transversion-only basis (an independent kind of thinning: it removes damage-prone transitions rather than sampling at random).

| basis | control holds? | Kotias breaks it? | do the six legit refs pass? |
|---|---|---|---|
| **898,057** | ✅ | ✅ | ✅ (max 2.8×) |
| **451,239** | ✅ | ✅ | ✅ (max 4.0×) |
| **451,387** | ✅ | ✅ | ✅ (max 2.3×) |
| **451,307** | ✅ | ✅ | ✅ (max 2.8×) |
| 239,106 *(transversions)* | ✅ | ✅ | ❌ **10.8×** |
| 225,737 | ✅ | ✅ | ❌ **16.9×** |
| 225,596 | ❌ | — | ❌ **13.7×** |
| 225,409 | ✅ | ✅ | ✅ (3.3×) |
| 108,341 | ✅ | ❌ | ✅ |
| 108,174 | ✅ | ❌ | ✅ |
| 108,219 | ✅ | ❌ | ✅ |

**Three regimes. The two bad ones fail in opposite directions.**

**≥ 450,000 SNPs — the screen works.** Four independent bases, four correct verdicts.

**~225,000–240,000 SNPs — FALSE POSITIVES.** In three of four bases, **a legitimate reference produces a larger effect than Kotias produces at full basis** (up to 16.9×). At this basis Test 3 would **convict an innocent population.** In one of the four, the control itself collapses.

**~108,000 SNPs — FALSE NEGATIVES.** In all three bases, **Kotias does not break the control at all.** The screen is blind. An analyst would clear a reference that is, in fact, an ancestor of their source.

> **Below ~450,000 SNPs, Test 3 does not lose power gracefully. It produces confident wrong answers — in one direction at intermediate basis, in the other at low basis. Running it below that threshold is worse than not running it.**

**This is not a caveat, and it does not belong in a discussion section.** A typical 1240K study at `maxmiss = 0` with several low-coverage populations works well below 450,000 SNPs. **The screen is not usable in most of the analyses it was designed to protect.**

### 2.7 The fold is retracted

**v1.0 reported a *fold* and argued that the fold, not the p-value, was the quantity to report — because the fold was stable across bases while the p-value was not.** It reported three: 43.8× / 41.9× / 41.9×.

**The fold is not stable.**

| quantity | range across 11 bases | spans |
|---|---|---|
| **fold** | 0.8× – 80.6× | **98×** |
| p (six refs) | 0.042 – 0.331 | 8× |
| p (six + Kotias) | 0.0028 – 0.289 | 104× |

**And it is unstable at a *fixed* basis size.** Three ~225,000-SNP replicates, differing only in which SNPs were drawn:

> **4.4× · 7.7× · 80.6×**

**An eighteen-fold swing from the SNP draw alone.** A ratio of two noisy p-values is noisier than either. Arithmetically obvious in hindsight; we did not check it.

**Why v1.0's claim looked true.** Its three "different" bases were **122,891 / 124,362 / ~123,000 SNPs**. **Three near-copies of one basis produce three near-identical folds. That is a tautology, and we presented it as the paper's methodological justification.**

> **Report the admissibility flip, not the fold. A threshold crossing survives noise that destroys a ratio.**

### 2.8 The exhibit we retract, and what it taught

**v1.0's Test 3 exhibit — `Greece_Mycenaean_POOL ← Greece_Minoan_POOL + Yamnaya`, 41.9× on 122,891 SNPs — is withdrawn in full.**

Both pools are bins:

| Pool | AADR group IDs | qpWave rank-0 |
|---|---|---|
| `Greece_Minoan_POOL` | **6** | **p = 3.2 × 10⁻⁸** |
| `Greece_Mycenaean_POOL` | **5** | **p = 6.3 × 10⁻⁵** |

**And the relationship was never true.** Rebuilt on coherent populations — `Myc_PELOPONNESE (18) ← Crete_Lasithi_EMBA (26) + Yamnaya` — the model **rejects at every basis from 108,442 to 970,005 SNPs** (p from 1.1 × 10⁻³ to 8.4 × 10⁻⁷). **Not a power effect. A misfit.**

**And we found what produced the apparent fit.** One basis, 956,579 SNPs:

| Model of `Myc_PELOPONNESE` | p | |
|---|---|---|
| `← Crete_Chania_LBA + Khvalynsk_Eneolithic` | **0.644** | **admissible** |
| `← Crete_Chania_LBA + Yamnaya` | **0.584** | **admissible** |
| `← Crete_Lasithi_EMBA + Yamnaya` | 2.7 × 10⁻⁶ | rejected |

**`Crete_Chania_LBA` is Late Minoan III — c. 1350–1225 BCE, a century and a half *after* the Mycenaean takeover of Knossos. Those Cretans already carry steppe ancestry.** Pooling them with pre-conquest Lasithi produced a "Minoan" source that **already contained the Mycenaeans.** The model needed less Yamnaya to reach the target — and the resulting 86.9 / 13.1 matched the published figure to three digits, **which we reported as external validation.**

> **It was not validation. It was circularity, and it looked exactly like validation.**

**Test 3 demands a truth established outside the model. We supplied one manufactured inside it** — by a pooling decision made for convenience, against the source publication's own site-level labels. *(See the companion paper, "Group IDs are not populations", §2.5.)*

### 2.9 Three limits, each found by hitting it

**(a) No purchase on a one-source clade.** `Russia_Khakassia_Afanasievo ← Russia_Samara_EBA_Yamnaya` (27 individuals, 2 sites, 400 yr; Narasimhan et al. 2019) **holds** with the six (p = 0.076, 975,248 SNPs) and **still holds** with Kotias (p = 0.095). **Test 3 does not fire — correctly.** Harney's bias acts on *admixture estimates*; where there is no mixture, there is no proportion to bias. **A two-source control is required. An analyst who tests a clade will wrongly conclude their reference set is clean.**

**(b) The Anatolian source cannot be screened, for a structural reason.** Both coherent European Neolithic descendants of `Turkey_N` reject as one-source models:

| Model | basis | p |
|---|---|---|
| `Slovakia_N_LBK ← Turkey_N` | 978,291 | **1.5 × 10⁻²⁰** |
| `Germany_HalberstadtSonntagsfeld_EN_LBK ← Turkey_N` | 975,336 | **3.6 × 10⁻¹⁵** |

**LBK carries hunter-gatherer ancestry**, so a WHG second source is needed — **and `Serbia_IronGates_Mesolithic` is in our reference set.** A WHG source on the left with a WHG-related forager on the right is *itself* a reference-into-source configuration. **The obstacle is the structure of European prehistory, not the archive.**

**(c) The Iran/CHG source cannot be screened at all.** No usable descendant of `Iran_GanjDareh` exists in AADR v66.p1: `Iran_SehGabi_C` is n = 7 across 871 years; `Iran_TepeHissar_C` spans 1,523 years and three platforms. **The hole is permanent — and it is the worst-placed one, because the reference we most want to screen against that source is Kotias, which is CHG.**

---

## 3. A χ² that rises when SNPs are removed at random

### 3.1 The claim v1.0 made, and why it was wrong

v1.0 observed that adding a low-coverage genome to an `extract_f2` call collapsed the SNP basis and reversed a model's conclusion. It then argued:

> *"χ² scales approximately with SNP count. A 17.7% loss therefore predicts χ² ≈ 12.5. Observed: 7.36 — a 51% fall. χ² falling faster than the SNP count is a diagnostic."*

**It also admitted, in its own caveats, that the inference had no null distribution:** ⟦OWED: a simulation calibrating the expected χ²/N ratio under uniform loss. Without it this is an observation, not a test with a threshold.⟧

**We built the null. The premise is false.**

### 3.2 Under random SNP loss, χ² rises

Take a **rejected** model on a clean public panel — `Germany_Esperstedt_CordedWare ← Yamnaya + Globular Amphora` with seven references (χ² = 13.93, dof 5, p = 0.016 at 898,057 SNPs) — and remove SNPs **at random**. Uniform loss by construction; this is exactly the null the diagnostic requires.

| random loss | "χ² ∝ N" predicts | **observed χ² ratio** (20 reps) |
|---|---|---|
| **17.7%** | 0.82 | **1.36** (range 1.15 – 1.52) |
| **30%** | 0.70 | **1.40** (1.04 – 1.73) |
| **50%** | 0.50 | **1.14** (0.40 – 1.65) |

**Throw away 18% of the SNPs at random and χ² goes *up* by a third. It rose in 20 of 20 replicates.**

**The mechanism is a small-sample property of the statistic, not of population history.** Fewer SNPs mean fewer jackknife blocks, which means a noisier estimate of the f4 covariance matrix Q, which means an inflated Q⁻¹ — so **χ² = E′Q⁻¹E is biased upward.**

**And the scaling law fails for accepted models too, in the same direction.** At 62% of the SNPs, the *accepted* model's χ² rose to **162%** of baseline and the *rejected* model's to **126%**. Neither tracks N.

> **v1.0 stated "χ² scales approximately with SNP count" unconditionally. It does not scale with SNP count at all. An analyst applying the χ²/N rule would be comparing a real statistic against an imaginary null.**

### 3.3 Which makes the exhibit stronger, not weaker

Same model, same reference set. **Only the extract changes** — a population is added that the **model never uses.**

| extract | N | χ² | p | |
|---|---|---|---|---|
| **core (10 pops)** | **898,057** | 13.93 | **0.016** | **REJECTS** |
| **+ `Turkey_Epipaleolithic`** | 746,953 **(−16.8%)** | **10.22** | **0.069** | **NO LONGER REJECTS** |
| **+ `Israel_Natufian`** | 639,313 **(−28.8%)** | **11.94** | 0.036 | still rejects |

**Now compare each against the random-loss null at the same basis:**

| added | naive "χ² ∝ N" prediction | **random-loss null (20 reps)** | **observed** | |
|---|---|---|---|---|
| `Turkey_Epipaleolithic` | 11.59 | median 18.34, **range 16.21 – 21.09** | **10.22** | **below all 20** |
| `Israel_Natufian` | 9.92 | median 18.70, **range 13.41 – 23.08** | **11.94** | **below all 20** |

> **χ² fell when random loss says it should have risen — below every one of twenty replicates, in both cases.**

**And the corrected null rescues a case the old rule would have missed.** `Israel_Natufian` — the genuinely low-coverage population (17 individuals, mean 107,584 SNPs) — produced a χ² of 11.94 against a naive prediction of 9.92. **It fell *less* than proportionally, so v1.0's rule would have seen nothing.** Against the real null it lies outside 20 of 20.

> **The premise was false, the diagnostic was weak, and correcting the premise makes the diagnostic strictly better.**

### 3.4 The restated diagnostic

> **Under random SNP loss, qpAdm's χ² rises. A χ² that *falls* as the basis shrinks cannot be produced by random loss. The surviving SNPs are an ascertained subset — the coverage footprint of whatever thin genome entered the extract.**
>
> **Report χ² and N together, and compare the observed χ² against a random-subsample null at the matched basis. The ratio χ²/N is meaningless.**

### 3.5 What sets the basis is not what you think

**`Turkey_Epipaleolithic` is a single individual with 891,783 of 1,233,013 SNPs covered. It is not a low-coverage genome.** It still costs **16.8% of the basis** — because at `maxmiss = 0`, every site *it* misses is deleted for *everyone*.

**And AADR contains far worse.** Reading the annotation's coverage column:

| group | n | mean SNPs on the 1240K panel |
|---|---|---|
| `Sudan_MiddleHolocene-1` | 1 | **534** |
| `Syria_EMBA_IA` | 1 | **679** |
| `Tanzania_LSA` | 2 | **817** |

**A population with 534 SNPs, dropped into an `extract_f2` call at `maxmiss = 0`, annihilates the basis for every other population in it.** Nothing in the qpAdm output says so.

### 3.6 Limits

- **Two thin populations, one model, one panel.** The upward bias of χ² under random loss is a property of the statistic and should generalise; the exhibit's magnitude need not.
- **We show the surviving SNPs are not a random subset. We do not show why.** GC content, mappability, and post-mortem damage are all plausible; none is tested.
- **This is not a threshold test.** It is a null distribution and an observation that lies outside it. Converting it into a calibrated diagnostic with a stated false-positive rate is the obvious next work, and we have not done it.

---

## 4. Recommendations

**Recommendations 3 and 5 in v1.0 were wrong, and are reversed here rather than removed.**

**For analysts**

1. **Positive-control your reference set.** Cladality and nested-fit screens are necessary and not sufficient. `Georgia_Kotias` passes both and still breaks a true relationship.
2. **Do it only at ≥ ~450,000 SNPs.** Below that the screen produces false positives at intermediate basis and false negatives at low basis. **If your basis is 120,000 SNPs, do not run Test 3 — its silence means nothing.**
3. ~~*Report the fold, not the threshold crossing.*~~ **REVERSED. Report the admissibility flip, not the fold.** The fold spans 0.8× to 80.6× across bases and swings eighteen-fold from the SNP draw alone at a fixed basis. **A ratio of two noisy p-values is noisier than either.**
4. **Publish the `extract_f2` population list, not just the model.** The basis is set by every population in the extract. A population that appears nowhere in your model can reverse your conclusion.
5. ~~*Report χ² alongside p and compare χ²/N.*~~ **REVERSED. Report χ² and N, and compare χ² against a random-subsample null at the matched basis.** χ² does not scale with N — under random loss it *rises*. **The χ²/N ratio has no null and should not be used.**
6. **Screen your populations for coverage before they enter an extract.** AADR contains groups with fewer than 1,000 SNPs on the 1240K panel. One of them in an extract at `maxmiss = 0` destroys the basis silently.
7. **Verify that a "known descendant" is what its label says.** A pool that fails qpWave rank-0 is not a population, and nothing can descend from it. **We violated this and it cost us the paper's headline.** See the companion paper.

**For archives and tool authors**

8. **Publish a per-group coverage summary** alongside locality, date, publication and platform. The field a user most needs before building an extract is the one least visible.
9. **`extract_f2` should warn when a population's coverage would cost more than a stated fraction of the basis.** The information is available at extract time. Nothing currently reports it.

---

## 5. On priority, and what is owed

**We do not claim novelty for either diagnostic.**

Six targeted searches, plus a further search of the qpAdm methods literature, returned no published equivalent of Test 3 — no protocol that validates a reference set against an **externally established** descent relationship. **That is not evidence of absence.** ⟦OWED: a citation search of the works citing Harney et al. (2021) and Flegontova et al. (2025) Fig. 1(a). It has not been done.⟧

The nearest prior art is qpAdm's own `dscore` block and its wrappers (§2.2). **They answer a different question, and we distinguish ourselves from them explicitly rather than leaving it to a referee.**

**This project has been pre-empted twice** — once on qpAdm failure modes (Harney et al. 2021), once on the cline argument (Flegontova et al. 2025). **We would rather find the third pre-emption ourselves.**

**One control fires. One.** A diagnostic demonstrated on a single relationship is a demonstration, not a validation. We have no second two-source control that both holds and fires, and we say so here rather than in a limitations paragraph.

---

## 6. Data, code and provenance

Scripts, logs and manifests: **https://doi.org/10.5281/zenodo.21333492**

**Every number in this paper is emitted by a versioned script into a log.** `manifest/` carries the `extract_f2` population list and SNP basis for every run. `manifest/RETRACTED.md` documents every analysis withdrawn during this work, and why.

**Every population in v2.0 is an AADR group ID taken whole.** No analyst-built pool, no unscripted QC filter, no private genome. **The entire paper is reproducible from AADR v66.p1 (Harvard Dataverse, doi:10.7910/DVN/FFIDCW) and the scripts in the repository.** This was not true of v1.0, and the fact that it was not true is what allowed most of what v1.0 got wrong.

**The author is not an academic and this work has not been peer reviewed.**

---

## References

- Haak W, Lazaridis I, Patterson N, et al. (2015). Massive migration from the steppe was a source for Indo-European languages in Europe. *Nature* 522: 207–211.
- Harney É, Patterson N, Reich D, Wakeley J (2021). Assessing the performance of qpAdm. *Genetics* 217(4): iyaa045.
- Flegontova O, Işıldak U, Yüncü E, Williams MP, Huber CD, Kočí J, Vyazov LA, Changmai P, Flegontov P (2025). Performance of qpAdm-based screens for genetic admixture on graph-shaped histories and stepping stone landscapes. *Genetics* 230(1): iyaf047.
- Patterson N, Moorjani P, Luo Y, et al. (2012). Ancient admixture in human history. *Genetics* 192: 1065–1093. *(ADMIXTOOLS; the `dscore` output.)*
- Skoglund P. *qpAdm_wrapper.* https://github.com/pontussk/qpAdm_wrapper
- Petr M, Vernot B, Kelso J (2019). admixr — R package for reproducible analyses using ADMIXTOOLS. *Bioinformatics* 35: 3194–3195.
- Maier R, Flegontov P, Flegontova O, et al. (2023). On the limits of fitting complex models of population history to f-statistics. *eLife* 12: e85492. *(ADMIXTOOLS 2.)*
- Narasimhan VM, Patterson N, Moorjani P, et al. (2019). The formation of human populations in South and Central Asia. *Science* 365: eaat7487.
- Wang C-C, Reinhold S, Kalmykov A, et al. (2019). Ancient human genome-wide data from a 3000-year interval in the Caucasus. *Nature Communications* 10: 590.
- Lazaridis I, Mittnik A, Patterson N, et al. (2017). Genetic origins of the Minoans and Mycenaeans. *Nature* 548: 214–218.
- Mallick S, Micco A, Mah M, et al. (2024). The Allen Ancient DNA Resource (AADR). *Scientific Data* 11: 182.
- Karalakis A (2026). *Group IDs are not populations: an audit of the Allen Ancient DNA Resource.* Companion paper, same repository.
