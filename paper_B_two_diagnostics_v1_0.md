> # ⚠ RETRACTED — superseded by v2.0
>
> **This version is retained for the record. Do not cite it. Three of its central
> claims are wrong, and they are named here rather than quietly removed.**
>
> **1. The "fold" is not a stable statistic.** This paper argued that the degradation
> fold, rather than the p-value, should be reported *because the fold was stable across
> SNP bases*. Measured across eleven genuinely different bases, **the fold spans 0.8× to
> 80.6× while the p-value spans 8×** — the fold is *twelve times less stable* than the
> quantity it was meant to replace. The three bases cited here as evidence of stability
> were 122,891 / 124,362 / ~123,000 SNPs: **three near-copies of one basis.**
>
> **2. The Test 3 exhibit was built on two bins and a circular relationship.**
> `Greece_Mycenaean_POOL` and `Greece_Minoan_POOL` each fail qpWave rank-0
> (p = 6.3 × 10⁻⁵ and 3.2 × 10⁻⁸). On coherent populations the relationship **rejects at
> every SNP basis from 108,442 to 970,005.** The apparent fit came from post-conquest
> Cretan individuals in the "Minoan" source who already carried steppe ancestry. **The
> 86.9 / 13.1 agreement with the literature was circularity, not validation.**
>
> **3. §3's premise is false.** χ² does not scale with SNP count. Under **random** SNP
> loss it **rises** (1.36× at 17.7% loss; 20 of 20 replicates). The §3 exhibit survives —
> but for the opposite reason to the one given here.
>
> **Superseded by `paper_B_v2_0.md`, which is rebuilt entirely from public AADR data.**
>
> ---
>
# Two diagnostics for qpAdm in real data

### A source-descendant test for reference-into-source gene flow, and a χ² signature of non-random SNP loss

**PREPRINT DRAFT v0.1. Not peer reviewed.**

> Both diagnostics below address gaps that the two major qpAdm performance papers **explicitly leave open**. Harney et al. (2021) and Flegontova et al. (2025) work in simulation, where the truth is known by construction. **Neither offers a way to detect their own documented failure modes in real data, where it is not.** These are two attempts.
>
> **⟦OWED⟧** marks a figure no script produces. **We claim priority for neither diagnostic without a literature search we have not completed** (§5).

---

## Abstract

Harney et al. (2021) show that gene flow from a reference population into a source, after the source's divergence from the true admixing lineage, **biases qpAdm's admixture estimates** — and that the bias is stronger when the gene-flow source sits in the reference set. Flegontova et al. (2025) list the same violation first in their exhaustive enumeration of the assumption failures that cause a true simple model to be rejected. **Both demonstrate it in simulation. Neither provides a diagnostic for real data.**

We propose **Test 3 (the source-descendant test)**: for each source *S*, identify a population *P* independently known to descend from *S*; fit *P* from *S* with and without the candidate reference *O*; **if *O* breaks a true descent relationship, *O* manufactures misfit against *S*, and *O* is rejected.**

Applied to a real analysis (AADR v66.p1 1240K, 122,891 SNPs, 708 blocks), Test 3 rejects a CHG reference (`Georgia_Kotias`) that **passes both standard screens** — qpWave cladality against every source (p from 3.9 × 10⁻⁷ to 3.5 × 10⁻¹³⁰), and a nested-fit screen in which it behaves exactly as a legitimate reference should (misfit rises 0.406 → 0.284; anchor SE **falls** 0.103 → 0.099; weights hold). Test 3 degrades a **true, admissible** descent relationship — `Mycenaean ← Minoan + Yamnaya`, which recovers an independently established 86.9 / 13.1 **to three digits** — by **41.9-fold** (p = 0.5113 → 0.0122). **All six legitimate references leave the relationship intact.** The fold is stable across three SNP bases (43.8× / 41.9× / 41.9×) while the p-value is not, and we report the fold for that reason.

Separately, we report a **χ² signature of non-random SNP loss**. Harney et al. show that under `allsnps: NO`, missing data collapses the usable SNP set — and state that the resulting behaviour is unbiased **"in cases where missing data is distributed randomly throughout the genome of all individuals."** We report the complement of that caveat. Adding **one low-coverage genome that appears in no model in the analysis** to an `extract_f2` call removed 17.7% of the SNP basis (124,362 → 102,314) and **flipped a conclusion from rejection to non-rejection** (p = 0.0044 → 0.118). Power loss predicts χ² ≈ 12.5. **Observed χ² = 7.36 — a 51% fall for an 18% loss of data**, reproduced in a second model (predicted 10.9, observed 6.17). The surviving SNPs are not merely fewer; **they carry disproportionately less of the misfit signal**, because they are the subset one low-coverage genome happens to cover.

**Test 3's limitations are severe and permanent, and we state them before its results** (§2.4).

---

## 1. The gap

Harney et al. name the violation. Flegontova et al. list it first among the assumption failures that cause rejection of true models, and show in simulation that a true clade relationship can be rejected by the rotating protocol due to left-to-right gene flows, with p-values from ~10⁻³² to ~10⁻⁵⁰.

**In simulation the analyst knows which reference is contaminated. In real data nobody does.** The standard remedies are:

- **Cladality screening** — is the reference non-cladal with every source? Necessary, and, we show, **not sufficient**.
- **Nested-fit screening** — does adding the reference behave as a legitimate reference should (misfit up, SE down, weights stable)? Necessary, and, we show, **not sufficient**.

A reference can pass both and still be wrong. Test 3 is an attempt to catch that case.

---

## 2. Test 3

### 2.1 The test

> **For each source *S*, identify a population *P* independently known to descend from *S*. Fit *P* from *S*, with and without the candidate reference *O*. If *O* breaks a true descent relationship, *O* manufactures misfit against *S*. Reject *O*.**

**The logic is direct.** If *O* received gene flow from *S* after *S* diverged from the true admixing lineage, then *O* is not symmetrically related to *S* and its descendants. Placing *O* in the reference set therefore generates f4 misfit against a relationship that is, by construction, true. **The misfit is the reference's, not the relationship's.**

**Test 3 is a positive control on the reference set.** It requires a relationship whose truth is established outside the model.

### 2.2 The baseline must be admissible — on its weights

> **No fold is computed until the baseline is proven admissible: all weights in [0,1], all |Z| > 2, model not rejected.**
>
> **A fold against a fallen model is not a weak result. It is not a result.** We learned this the expensive way: an earlier version of this analysis computed degradation ratios of 55,765× and 2,536× against baselines whose steppe weights were **significantly negative** (Z −3.57 to −6.65), accepted because their p-values (0.58, 0.41, 0.25) looked healthy. All were void. **The check is on the weights, and it comes first.**

### 2.3 Worked example

**Basis.** AADR v66.p1 1240K; `maxmiss = 0`, `blgsize = 0.05`; 16 populations, **122,891 polymorphic SNPs**, 708 jackknife blocks; `admixtools2` v2.0.10. **Every figure in this section is from that one extract.** Reference set: Mbuti, Italy_Sicily_Epigravettian, Russia_Kostenki_UP, Papuan, Serbia_IronGates_Mesolithic, Morocco_Iberomaurusian.

**The candidate reference.** `Georgia_Kotias_QC` (CHG) — an obvious remedy for a reference panel with no power on the Iran/CHG axis that separates the sources.

**Screen 1 — cladality. PASSES.**

| Kotias vs | qpWave p |
|---|---|
| Greece_Minoan_POOL | 1.07 × 10⁻¹⁸ |
| Turkey_EBA_WEST | 3.90 × 10⁻⁷ |
| Turkey_N_QC | 1.18 × 10⁻²⁷ |
| Iran_GanjDareh_N_QC | 5.77 × 10⁻²⁰ |
| Russia_Eneolithic_Steppe_QC | 3.54 × 10⁻¹³⁰ |
| Russia_Yamnaya_QC | 2.30 × 10⁻⁵⁹ |

**Screen 2 — nested fit. PASSES, textbook.**

| `MyDNA ← Minoan + Iran + Eneolithic` | six refs | + Kotias |
|---|---|---|
| model p | 0.406 | **0.284** *(misfit rises)* |
| anchor weight | 62.3% | 63.9% *(holds)* |
| anchor SE | 0.103 | **0.099** *(falls)* |
| Iran / steppe | 20.7 / 17.0 | 18.4 / 17.6 *(hold)* |

**Test 3 — REJECTS.**

| Descent relationship | six refs | + Kotias | fold |
|---|---|---|---|
| **Mycenaean ← Minoan + Yamnaya** *(independently 86.9 / 13.1)* | **0.5113**; recovers **86.9 (Z 28.21) / 13.1 (Z 4.24)** — admissible | **0.0122** | **41.9×** |

**Discrimination.** Leave-one-out over the six legitimate references — p without each: 0.366, 0.391, 0.456, 0.582, 0.370, 0.686. **All leave the true relationship intact.** With all six: **0.511**. Adding Kotias: **0.0122**.

> **Every legitimate reference leaves the relationship standing. Test 3 rejects exactly one — the one that was wrong.**

**Why report the fold rather than the threshold.** The fold is **basis-stable**: 43.8× / 41.9× / 41.9× across three SNP bases. The p-value is **not**: it moved 23% (0.0146 → 0.0113) under a change of extract that added two unrelated populations. **The stable quantity is the degradation, not the crossing.**

### 2.4 Limitations — severe, and stated before use

> **A Test-3 failure is decisive. A Test-3 pass is partial.**

- **Test 3 clears a reference only against sources for which a known descendant exists in the dataset.** In our case study, `Iran_GanjDareh` has none. The steppe control is void (Yamnaya cannot be modelled from Eneolithic steppe under any configuration: cladality 3.97 × 10⁻⁶⁵; +Turkey_N 7.29 × 10⁻¹²; +Iran 1.96 × 10⁻⁵). **Our six references are cleared on one source of three — and the unguarded source is `Turkey_N`, the largest contributor to the target.** That is not a detail: the source most able to distort the result is the one the diagnostic could not screen.
  - We tested the specific threat this creates and it did not materialise. `Serbia_IronGates_Mesolithic`, one of our six references, spans 11,465 → 7,803 BP, and farming reaches the Danube Gorges ~8,000 BP — so its late fraction could plausibly carry Anatolian-farmer ancestry, making a reference share drift with a source. It does not: f4(Mbuti, Turkey_N; IG_early, IG_late) gives Z = **−1.67** (n.s.), and `IronGates_LATE ← IronGates_EARLY + Turkey_N` returns Turkey_N at **−4.7% (Z = −2.95)** — significantly *negative*. **The late fraction is further from farmers, not closer.** The reference survives; the structural gap in the diagnostic remains.
- **It requires an externally established truth.** Where none exists, the test cannot be run — and where the "known descendant" is itself mislabelled, the test will be run on a relationship that does not exist. **This happened to us.** The population we first used as the Anatolian-source descendant was **a northern Levantine city carrying a country-based archive label**; every fold computed from it was void. (See the companion paper, *Group IDs are not populations*.)
- **It does not survive Flegontova's own α = 0.01 threshold** in our example: p = 0.0122 rejects at 0.05, not at 0.01.
- **The mechanism is not established.** We cannot say *why* Kotias breaks the relationship — whether the damage tracks the source's post-split CHG intake, as the Harney violation predicts. The controls that would isolate it were all inadmissible. **Test 3's result does not depend on the mechanism; the diagnostic works without it. But we cannot offer the mechanism, and we do not.**
- **Demonstrated on one dataset, by the analyst who devised it, on a case where the answer was already known.** **Its value depends on it firing on someone else's published analysis. That has not been done, and it is the highest-value next step.**

### 2.5 It applies to the rotating design

Rotation **deliberately places candidate sources into the reference set** — precisely where a post-split gene-flow relationship does its damage. **A rotating set must itself be Test-3 screened**, which to our knowledge is not current practice.

---

## 3. A χ² signature of non-random SNP loss

### 3.1 Harney's caveat, and its complement

Harney et al. examine the `allsnps` option directly, and quantify the collapse: at 80%, 85% and 90% missing data, `allsnps: NO` leaves an expected 181,988 / 37,304 / 1,610 SNPs, with standard errors rising accordingly. They recommend `allsnps: YES`, on the grounds that it improves model discrimination **without creating biases *in cases where missing data is distributed randomly throughout the genome of all individuals*.**

**That conditional is the whole of this section.** When missingness is driven not by uniform degradation but by **one low-coverage genome**, the surviving SNPs are **that genome's coverage footprint** — an ascertained subset, not a random one. Harney et al. bracket this case out. We report it.

### 3.2 The exhibit

`extract_f2` at `maxmiss = 0` fixes the SNP basis using **every population in the extract, not only those in the model.** Against a 10-population core extract of **124,362** SNPs:

| Added to the core | SNPs | cost |
|---|---|---|
| Georgia_Kotias_QC | 124,377 | +15 |
| Turkey_PPN_QC | 124,372 | +10 |
| Iran_GanjDareh_N_QC | 122,862 | −1,500 (−1.2%) |
| **Turkey_Epipaleolithic** (n = 1) | **102,314** | **−22,048 (−17.7%)** |
| Israel_Natufian_QC (n = 1) | 38,041 | −86,321 (−69.4%) |

**Now fix the model. Fix the populations *in* the model. Change only the extract** — by adding a population **the model never uses:**

| Extract | SNPs | Model | **χ²** (dof 4) | p | |
|---|---|---|---|---|---|
| core, 10 pops | **124,362** | Minoan + Eneolithic | **15.14** | **0.0044** | **rejects** |
| + Turkey_Epipaleolithic | 102,314 | *(identical)* | **7.36** | 0.118 | **does not reject** |
| core, 10 pops | **124,362** | Minoan + Yamnaya | **13.23** | **0.0102** | **rejects** |
| + Turkey_Epipaleolithic | 102,314 | *(identical)* | **6.17** | 0.187 | **does not reject** |

**A conclusion reverses because of a population that appears nowhere in the model.**

### 3.3 The signature

χ² scales approximately with SNP count. A 17.7% loss therefore predicts **χ² ≈ 12.5**. **Observed: 7.36 — a 51% fall.** The second model reproduces it: predicted 10.9, observed 6.17.

> **χ² falling faster than the SNP count is a diagnostic.** It distinguishes *"there is no signal"* from *"the signal-bearing sites were dropped."* p alone cannot tell them apart, and the difference decides whether a non-rejection means anything.

**Practical rule.** When a model's p-value changes after an extract change, **report the χ² and the SNP count, and compare their ratios.** If χ²/N held constant, the change is power. If χ² fell faster, the surviving SNP set is compositionally different, and the non-rejection is an artefact of ascertainment.

### 3.4 Caveats

- **χ² ∝ N is an approximation**, and the f4 covariance matrix is itself noisier at 102k SNPs. The effect is large enough to survive that, but we claim the observation is **consistent with** compositional bias, not that it proves it.
- ⟦**OWED**: a simulation calibrating the expected χ²/N ratio under uniform vs. footprint-driven SNP loss. **Without it this is an observation, not a test with a threshold.** It is the obvious next piece of work and we have not done it.⟧
- **We ran the entire case study at `maxmiss = 0`**, i.e. in the `allsnps: NO` regime Harney et al. advise against — while Flegontova et al. advise the opposite remedy (use large SNP sets with low missing-data rates, i.e. **exclude thin populations** rather than patch around them). **We did not choose between these prescriptions, and we should have.** Our single catastrophic basis collapse came from admitting precisely the low-coverage genome both papers warn about.

---

## 4. Recommendations

1. **Test-3 screen every reference**, and every rotating set. Cladality and nested-fit screens are necessary, not sufficient.
2. **Prove a baseline admissible on its weights before computing any degradation ratio.**
3. **Report the fold, not the threshold crossing**, when the p-value is basis-fragile.
4. **Report the `extract_f2` population list.** The basis is set by every population in the extract, not only those in the model.
5. **Report χ² alongside p** whenever the SNP basis changes, and compare χ²/N.
6. **Verify that a "known descendant" is what its label says** before using it in Test 3, or in any control. See the companion paper.

---

## 5. On priority

**We have not completed a systematic literature search, and we do not claim priority for either diagnostic.** Six targeted searches did not return a published equivalent of Test 3, and Harney et al.'s randomness caveat appears not to have been revisited. **Neither is evidence of novelty.** Before either claim is made, a citation search of Harney et al. (2021) Fig. 5 and Flegontova et al. (2025) Fig. 1(a) is required.

**This project has already been pre-empted twice** — once on qpAdm failure modes (Harney), once on the cline argument (Flegontova). **We would rather find the third pre-emption ourselves.**

---

## Data and code availability

**Scripts, logs and manifests: https://doi.org/10.5281/zenodo.21333492** (archived from https://github.com/alexandroskaralakis/aadr-group-ids).

Every number in this paper is emitted by a versioned script in `scripts/` into a log in `logs/`. `manifest/` carries the `extract_f2` population list and SNP basis for each run — the quantity §3 shows to be decisive. `manifest/RETRACTED.md` documents the withdrawn analyses, including the basis-polluted run from which §3's exhibit was built.

**Genotype data are public** (AADR v66.p1, Harvard Dataverse, doi:10.7910/DVN/FFIDCW).

**The target genome is the author's own.** It is not deposited publicly, because it carries information about living relatives who have not consented to its release. It is available to researchers on request.

**A reproducibility gap.** The working panels were assembled interactively with EIGENSOFT before this repository existed; their construction is not scripted. All downstream analysis is reproducible given the panels; the panels themselves are not.

## References

- Harney É, Patterson N, Reich D, Wakeley J (2021). Assessing the performance of qpAdm. *Genetics* 217(4): iyaa045.
- Flegontova O, Işıldak U, Yüncü E, Williams MP, Huber CD, Vyazov LA, Changmai P, Flegontov P (2025). Performance of qpAdm-based screens for genetic admixture on admixture-graph-shaped histories and stepping-stone landscapes. *Genetics* 230(1): iyaf047.
- Maier R, Flegontov P, Flegontova O, Işıldak U, Changmai P, Reich D (2023). On the limits of fitting complex models of population history to f-statistics. *eLife* 12: e85492.
- Lazaridis I et al. (2017). Genetic origins of the Minoans and Mycenaeans. *Nature* 548: 214–218.
- Mallick S et al. (2024). The Allen Ancient DNA Resource (AADR). *Scientific Data* 11: 182.
- [Companion paper] *Group IDs are not populations: an audit of the Allen Ancient DNA Resource.* Same repository: https://doi.org/10.5281/zenodo.21333492
