# Group IDs are not populations

### An audit of the Allen Ancient DNA Resource, and a worked case in which an archive label determined a qpAdm result end to end

**PREPRINT DRAFT v1.2. Not peer reviewed.**

> **v0.2 has results.** The audit ran on AADR v66.p1: **23,089 individuals, 3,897 group IDs, 2,090 with n ≥ 2.**
>
> **Two of the seven original flags are withdrawn.** One was circular; one was measuring our own chronology table rather than the archive. Both are reported and neither is counted. **This is stated here rather than quietly fixed, because a paper about unexamined metadata cannot ship unexamined metrics.**
>
> **A second worked case is added** (§2.5), and it is not the archive's fault. `Greece_Minoan_POOL` — a pool the *analyst* built — is a bin of six sites spanning a millennium and straddling the Mycenaean conquest of Crete. **AADR gives its 58 individuals six separate site-level group IDs; we overwrote all six with one.** It fails qpWave rank-0 at **p = 3.2 × 10⁻⁸**.
>
> **v1.1 reframed the paper around Eisenmann et al. (2018).** That paper did not merely warn about country labels; it proposed a **coherence criterion**, and asked archives and studies to publish the tables that would let anyone check it. **This audit is the first measurement of that criterion.** Two of Eisenmann et al.'s eleven authors maintain the AADR.
>
> **v1.2 rebuilds §2.5 on an ancients-only basis, and retracts its previous p-values.** The old figures (p = 0.0074 and p = 0.0014) were computed on a panel into which the target genome had been merged, so their SNP set was capped by a consumer genotyping chip — even though neither test contains the target. **Rebuilt from AADR alone at ~950,000 SNPs, both splits reject far harder: p = 3.2 × 10⁻⁸ and p = 4.2 × 10⁻¹⁰.** The conclusion is unchanged and the evidence is stronger, but **the retracted numbers are named rather than quietly replaced.** §2.5 is now reproducible from the public archive; **the scripts are in the repository and touch no private data.**
>
> **v1.2 also withdraws the claim that "the fault line is temporal."** Site and date are perfectly confounded in this pool — there are no individuals in a 243-year window between the two clusters — and no design on this sample can separate them.
>
> **⟦OWED⟧** marks two figures no script yet produces: a regionalised chronology table, and a systematic count of within-group duplicate genotype calls. Both are named where they arise.

---

## Abstract

A qpAdm source population is a string in an `.ind` file. In archaeogenetic practice those strings come, very largely, from the group IDs distributed with the Allen Ancient DNA Resource (AADR). **Group IDs are curatorial conveniences. They are not population labels, and they are not audited as such.**

**The field already knows this, and said so.** Eisenmann et al. (2018) proposed five criteria for any archaeogenetic nomenclature — brevity, **coherence**, accessibility, flexibility, stability — where coherence means that <q>individuals from genetically distinguishable groups should not be given the same name</q>. They warned specifically that modern nation-state names are politically charged and should not be assumed to track genetic clusters, and they asked that every study publish a table giving each sample's site, date and cluster assignment. **Two of that paper's eleven authors are David Reich and Iosif Lazaridis.** **AADR's group IDs are nonetheless overwhelmingly of the form `Country_Period`, and no coherence table is published with any release.**

**This paper is the first measurement of that coherence criterion.** We audit all 2,090 AADR v66.p1 group IDs with n ≥ 2 against the archive's own annotation. **One group in five mixes genotyping platforms. One in six spans three or more localities. One in eleven spans eight centuries or more. One in twelve spans 300 km or more. One in seven trips two or more independent coherence flags.**

We then report what one such label cost. **`Turkey_MLBA` is Alalakh (Tell Atchana), Hatay — a northern Levantine city 25 km from the Syrian border and ~900 km from the Aegean.** The publication that generated most of these genomes classifies it as *Northern Levant in its title*. AADR names it `Turkey_` because Hatay is inside the modern Republic of Turkey. In a personal-ancestry qpAdm analysis carried through nineteen SNP bases and six drafts, this one label supplied the ground-truth control, the f4 negative control, and the study's central negative result — *"a design that rejects a true descent relationship cannot be used to adjudicate a disputed one."* **The relationship was not true. qpAdm was right to reject it.** Every fit, standard error and p-value behaved normally throughout.

**And the audit cannot see it.** `Turkey_MLBA` is one site, one platform, tightly dated, geographically pointlike: **it collects zero flags.** The automated audit catches *incoherent* groups; it is blind to *misnamed coherent* ones — which are the more dangerous kind, because they behave perfectly and lie only in the label. **All prevalence figures here are lower bounds, and we say so in the abstract rather than the limitations.**

Finally, we turn the audit on the analyst. **`Greece_Minoan_POOL`, a pool we built ourselves, is also a bin** — six sites, 1,042 years, straddling the fall of Knossos — and it fails qpWave rank-0 tests at **p = 3.2 × 10⁻⁸** by site and **p = 4.2 × 10⁻¹⁰** by date. **AADR gives those 58 individuals six separate site-level group IDs. Lazaridis et al. (2017) had deliberately kept the Cretan sites apart, and Eisenmann et al. (2018) held that separation up as an example of good practice. We overwrote six correct published labels with one wrong label of our own, because it was convenient.** **The archive supplied one bad population. We supplied the other — against explicit published advice, and against the archive's own correct labels, in the very field we were auditing.**

---

## 1. The problem is old, the criterion exists, and nobody measured it

Eisenmann et al. (2018) — writing from the Max Planck Institute for the Science of Human History and Harvard, and **including David Reich and Iosif Lazaridis among its eleven authors** — set out what an archaeogenetic nomenclature has to do. They proposed five criteria: **brevity, coherence, accessibility, flexibility, stability.**

**Coherence is the one this paper measures.** Their definition: <q>individuals from genetically distinguishable groups should not be given the same name</q> — and, conversely, that genetically indistinguishable individuals should share one. They warned that present-day country names are political constructs, arose largely after the eighteenth century, are ethno-nationally charged, and must not be assumed a priori to correspond to genetic clusters. They recommended that a cluster be named for a region only when its samples actually span that region, and that **every study publish a table listing each individual's site, relative and absolute date, and cluster assignment.**

**None of this was implemented in the archive that the field then built its analyses on.** AADR's group IDs are overwhelmingly `Country_Period`. No coherence table accompanies any release. The AADR descriptor states that the project errs on the side of inclusivity, admitting data even where meta-information is incomplete, and relies on **user-community feedback** to identify erroneous meta-information. **The archive is explicit that its metadata is community-corrected rather than systematically audited, and it invites the correction.**

**In seven years nobody has run the audit.** This paper is it — plus the case that forced it, and the two bins we found in our own work.

**This is not a charge of carelessness.** AADR's genotype curation is careful and documented, and **the disconfirming evidence is distributed with the archive** — in the Locality field, one column from the one we used. **The failure is an interface failure: analysts build populations from Group IDs and never open Locality.** We did exactly that, for six drafts. *(The corrections in §2 have been reported to the AADR maintainers.)*

---

## 2. Results

### 2.1 The audit

**Each flag below is an operationalisation of Eisenmann et al.'s coherence criterion**: a group that spans many sites, many centuries, many publications or many platforms is a group in which genetically distinguishable individuals may be sharing one name.

AADR v66.p1, 1240K. **23,089 individuals; 3,897 group IDs; 2,090 with n ≥ 2** (`run_aadr_audit.R`).

| Flag | groups | prevalence |
|---|---|---|
| **≥ 2 genotyping platforms** | 416 | **19.9%** |
| **≥ 3 distinct localities** | 348 | **16.7%** |
| **temporal span ≥ 800 years** | 192 | **9.2%** |
| **geographic span ≥ 300 km** | 178 | **8.5%** |
| **≥ 3 publications** | 86 | **4.1%** |

| **≥ 2 of the 5 flags** | 290 | **13.9%** |
| **≥ 3 of the 5 flags** | 116 | **5.6%** |

> **One AADR group in seven trips two or more independent coherence flags. One in eighteen trips three or more. Eighteen groups trip all five.**

**The worst offenders are not obscure.** All five flags, with the two most extreme columns:

| Group | n | localities | pubs | platforms | span km | **span yr** |
|---|---|---|---|---|---|---|
| **Spain_C** | 112 | **23** | **10** | 3 | 724 | 1,854 |
| **France_Mesolithic** | 18 | 9 | 4 | 2 | 724 | **5,840** |
| **Turkey_N** | 68 | 6 | **7** | 3 | 743 | **3,046** |
| **Turkey_PPN** | 32 | 4 | 6 | 2 | 482 | **3,247** |
| **Scotland_N** | 74 | **20** | 6 | 3 | 458 | 1,898 |
| **Spain_EBA** | 72 | **19** | 4 | 2 | 714 | 883 |
| **Kazakhstan_IA** | 14 | 5 | 3 | 2 | **2,296** | 1,201 |
| **Turkey_EBA** | 30 | 8 | 4 | 3 | 1,293 | 1,350 |
| Ukraine_Mesolithic | 19 | 4 | 4 | 2 | 199 | **5,411** |

**`France_Mesolithic` spans 5,840 years. `Ukraine_Mesolithic`, 5,411. These are not populations; they are epochs.** `Spain_C` pools 112 individuals from 23 sites across 10 publications.

> **And `Turkey_N` — 68 individuals, 6 localities, 7 publications, 3 platforms, 3,046 years — is the primary Anatolian source of the case study in §2.3–2.5.** The analyst used it for six drafts. It trips all five flags.

### 2.2 Two flags withdrawn, and why

**`F_multipoly` (≥ 2 political entities) is circular.** Group IDs are `Country_Period`; the Political Entity field *is* the country. The flag tests only whether AADR agrees with itself. It fired on **1 group in 2,090** — and that one, `BantuSA`, is modern. **Useless by construction.**

**`F_period` was measuring our chronology table, not the archive.** It fired on **74 of 74** `Scotland_N` individuals, **46 of 46** `England_N`, **11 of 11** `Wales_N`, and on 2,210 individuals overall. **A flag that fires on 100% of a group is not detecting an error; it is announcing that the window is wrong.** Our `_N` window (7,000–11,000 BP) is Near-Eastern; the northwest European Neolithic runs ~7,500–4,000 BP, so `Belgium_N` at 4,725 BP is *correct* and our table was *wrong* — as is `Bulgaria_Dzhulyunitsa_EBA` at 5,153 BP, correct for the Balkan Early Bronze Age. The v0.1 figure of 13.3% is **withdrawn**.

> **But the screen was not worthless, and we say so rather than bury it.** Among the false positives it produced a real one: **`Albania_Cinamak_IA` contains an individual dated 525 BP — approximately 1425 CE, medieval — in a group labelled Iron Age.** ⟦OWED: a regionalised chronology table, after which this check becomes usable. Until then it measures the analyst, and its prevalence figure is not reported.⟧

> **We report this rather than fix it quietly, because a paper arguing that unexamined metadata determines results cannot ship unexamined metrics.** The period check remains a good idea; it needs region-specific windows, and until it has them it measures the analyst.

### 2.3 The exhibit: `Turkey_MLBA` is Alalakh

All individuals carry the prefix `ALA`. From AADR's own annotation:

| Field | Value |
|---|---|
| **Group ID** | `Turkey_MLBA` |
| **Locality** | **"Atchana-Alalakh (Hatay, Reyhanlı)"** (28) / **"Tell-Atchana (Hatay, Reyhanlı)"** (10) |
| **Political Entity** | Turkey |
| Dates | 3337–3889 BP (≈ 1940–1390 BCE) |
| Publications | Skourtanioti/Krause *Cell* 2020; Ingman/Stockhammer *PLoS ONE* 2021 |

**Tell Atchana — ancient Alalakh — lies in the Amuq valley, Hatay, 25 km from the Syrian border and ~900 km from the Aegean.** The publication reporting most of these genomes is titled *"Genomic History of Neolithic to Bronze Age Anatolia, **Northern Levant**, and Southern Caucasus."*

> **The source publication classifies it as Northern Levant. The archive labels it Turkey. Both are correct. Only one is a population.**

*(A minor defect worth its own line: the one site carries two Locality strings, so any script grouping by Locality splits a single population in two.)*

**A third defect class, found incidentally.** The period screen surfaced **duplicate genotype calls of the same individual inside one group ID** — `I17183.AG` and `I17183.SG` both in `Armenia_Keti_MBA`; `I2520.AG`, `I2520.DG` and `I2520.SG` all in `Bulgaria_Dzhulyunitsa_EBA`. **An analyst pooling by Group ID pseudo-replicates that individual two- or three-fold.** ⟦OWED: a systematic count of within-group duplicate calls across v66.p1. It is one `sed` and one `uniq -d`, and we have not run it.⟧

### 2.4 What the audit cannot see — the paper's central limitation

| | n | localities | pubs | platforms | span km | span yr | **flags** |
|---|---|---|---|---|---|---|---|
| **Turkey_EBA** | 30 | **8** | 4 | 3 | **1,293** | 1,350 | **5 of 5** |
| **`Turkey_MLBA`** | 38 | 2 | 2 | 1 | **0** | 552 | **0 of 5** |

**Zero flags.** By every automated metric, `Turkey_MLBA` is an *exemplary* population: one site, one platform, tightly dated, geographically pointlike. **It is the group that destroyed six drafts of analysis, and the audit is blind to it.**

> **The audit finds incoherent groups. It cannot find misnamed coherent ones — and those are the more dangerous kind, because they behave perfectly in every statistic and lie only in the label.** Detecting `Turkey_MLBA` required reading the Locality field and knowing where Hatay is. **Every prevalence figure in §2.1 is therefore a lower bound.**

**And `Turkey_EBA` is an archive bin** — eight localities across ~1,293 km, four publications, three platforms — **including `ulu117.SG`, whose own contextual date (4000–3000 BCE; 5450 ± 289 BP) places it in the Late Chalcolithic**, before the Anatolian Early Bronze Age begins. The remaining individuals fall between ~4,439 and 5,166 BP.

**A second blind spot, of the opposite kind.** `Greece_Crete_LBA` holds 23 individuals from **Chania** and 2 from **Heraklion** (`kro008.SG`, `kro009.SG`; Skourtanioti/Stockhammer 2023) — two localities ~140 km apart under one group ID.

**It trips the "≥ 2 localities" condition. Our audit flags at "≥ 3."** `Greece_Crete_LBA` therefore scores **zero flags**, exactly as `Turkey_MLBA` does, for the opposite reason: it is a *small* bin, and small bins fall below the threshold.

**And it cannot be tested.** A qpWave of Chania against Heraklion has n = 2 on one side. **It would not reject, and a non-rejection would mean nothing.** The remedy this paper recommends — *a pool is a hypothesis, test it* — has no power here.

> **A two-locality bin with a two-individual minority is invisible to the metric and untestable by the method. We report it because we happened to read the Locality column — which is this paper's entire recommendation, and its entire limitation.**

### 2.5 And now the analyst: `Greece_Minoan_POOL`

**The pool we built ourselves is a bin — and this time the archive is blameless.**

`Greece_Minoan_POOL` is 58 individuals. **AADR gives those 58 individuals six different group IDs:**

| AADR Group ID | n | Locality | Date range (BP) |
|---|---|---|---|
| `Greece_Crete_HgCharalambos_EMBA` | **26** | Lasithi (Ierapetra, Schinokapsala, Vornospilia) | 3,933 – 4,177 |
| `Greece_Crete_LBA` | **25** | **Chania (23) + Heraklion (2)** | 3,135 – 3,300 |
| `Greece_Crete_MoniOdigitria_BA` | 4 | Moni Odigitria (Heraklion) | 3,900 |
| `Greece_Crete_Aposelemis_LBA` | 1 | Heraklion | 3,525 |
| `Greece_Crete_Zakros_Karaviadaina_BA` | 1 | Karaviadaina (Zakros) | 3,768 |
| `Greece_Crete_NecropolisofArmenoi_BA` | 1 | Necropolis of Armenoi (NW coast) | 3,300 |

**Hagios Charalambos is the cave on the Lasithi plateau.** The archive named it by site and period. It named all six by site and period. **Four of them are single individuals, kept as singletons because they are singletons.**

**We overwrote six published labels with one label of our own, and then built six drafts on the result.**

**And the field had already told us not to.** Lazaridis et al. (2017) deliberately did not pool the Cretan Bronze Age samples; they named them by site. Eisenmann et al. (2018) cite that decision **approvingly, by name**, as an illustration of naming a cluster only at the geographic scale its sampling supports. **Both the source publication and the archive kept these sites apart. We merged them because one string is easier to type into a `left` set than six.**

**The test.** Sub-pools are defined by **AADR's own Group IDs** — never by Locality string-matching, never by genetic similarity. Two arms:

- **By site.** `Crete_Lasithi_EMBA` = `Greece_Crete_HgCharalambos_EMBA` (26). `Crete_Chania_LBA` = `Greece_Crete_LBA` (23), **excluding `kro008.SG` and `kro009.SG`**, which are Heraklion, ~140 km away (§2.4). The exclusion is declared so that a rejection means *Lasithi ≠ Chania* and nothing else.
- **By date.** All 58, split at **3,600 BP** — a threshold fixed before any f-statistic was computed. Knossos falls to the Mycenaean mainland ~1450 BCE (~3,400 BP).

Right set in both arms: `Mbuti`, `Italy_Sicily_Epigravettian`, `Russia_Kostenki_UP`, `Papuan`, `Serbia_IronGates_Mesolithic`, `Morocco_Iberomaurusian`. Two left, six right → dof 5.

**qpWave rank-0** (`build_minoan_panel.sh` → `run_minoan_split_clean.R`; panels built ancients-only from AADR v66.p1 via `convertf`, `maxmiss = 0`, a fresh extract per arm, **no target genome in either panel**):

| Split | n vs n | χ² (dof 5) | **rank-0 p** | basis | |
|---|---|---|---|---|---|
| **Lasithi vs Chania** | 26 vs 23 | 43.32 | **3.2 × 10⁻⁸** | 950,383 | **REJECTED** |
| **Early vs Late** (3,600 BP) | 31 vs 27 | 52.55 | **4.2 × 10⁻¹⁰** | 955,096 | **REJECTED** |

> **Our own project record said "Minoan pool homogeneous, qpWave rank 0, p = 0.09."** It is not homogeneous. **It rejects at p = 3 × 10⁻⁸.**

*(The site arm was extracted twice, in independently created directories, and returned χ² = 43.32425 both times.)*

**The threshold cannot be tuned.** **There are no individuals between 3,525 and 3,768 BP.** The early cluster runs 3,768–4,177; the late cluster runs 3,135–3,525. **The pool contains a 243-year void, and our pre-registered threshold of 3,600 BP falls inside it. Every threshold from 3,526 to 3,767 BP produces the identical partition.** The date arm is not a cut we placed; it is a gap the data already had.

**What we cannot show, and will not claim.** **Site and date are perfectly confounded in this pool.** Lasithi *is* the early cluster; Chania *is* the late one. The date arm assigns all 49 site-arm individuals exactly as the site arm does, and adds nine others. The two tests are not independent, and **the higher χ² of the date arm reflects nine additional individuals, not a sharper fault line.**

> **We therefore claim only this: `Greece_Minoan_POOL` is not one population. We do not claim to know whether the structure is spatial, temporal, or — as is most likely on a small island over a millennium — both. The design cannot separate them, and nothing in these data can.**

⟦OWED⟧ **Whether the late fraction carries steppe ancestry from the Mycenaean takeover.** On the superseded basis, `Minoan_LATE ← Minoan_EARLY + Yamnaya` returned 8.6% Yamnaya at Z = 1.54 (n.s.; model rejected, p = 0.003). **That result is withdrawn together with its basis and has not been recomputed on the clean panel.** The pool is heterogeneous; the cause is not established.

**The retraction.** **Earlier drafts of this section reported p = 0.0074 (site, 88,776 SNPs) and p = 0.0014 (date, 121,732 SNPs). Both are withdrawn.** Neither was wrong in direction, but both were computed on a panel into which the study's target genome had been merged — so the SNP set had been intersected down to a consumer genotyping chip *before any extract ran*. **A test between two ancient sub-pools, neither of which contains the target, was nonetheless conditioned on the target's chip.** The site arm was further degraded by carrying `Minoan_OTHER` — nine individuals its own model never used — into the extract, at a cost of ~33,000 SNPs.

**Every population in an `extract_f2` call sets the basis, not only those named in the model.** We wrote that rule, enforced it on others, and broke it here.

> **This is the paper's fairest exhibit and its most uncomfortable. `Turkey_MLBA` was the archive's error. `Greece_Minoan_POOL` is entirely ours — committed against explicit, published, seven-year-old advice, and against the archive's own correct labels, by an analyst who was at the time auditing SNP bases to the last digit and pre-registering every prediction.**
>
> **The failure mode is active, not passive.** It is not that archives lump. **It is that analysts re-lump what archives and publications carefully separated**, because one label is easier to type than six. Every safeguard downstream of that keystroke is irrelevant.

### 2.6 Two negative controls, both reported

Not every suspicion survived contact with the data. Both of the following were pre-registered predictions, and both **failed** (`run_anchor_integrity.R`):

- **`Turkey_N_QC` contains `Kumtepe` (5,615 BP), ~2,400 years younger than the rest of the pool and Late Chalcolithic — in a group named Neolithic.** We predicted it would contaminate the analysis's headline model. **It does not.** With Kumtepe in (34 individuals, 122,704 SNPs): **60.1 / 26.2 / 13.7**, p = 0.775. With it out (33, 122,702): **60.0 / 26.2 / 13.7**, p = 0.767. **Shift: 0.1 points.** *(Note both panels differ; the two rows are compared to each other and to nothing else.)*
- **`Serbia_IronGates_Mesolithic` spans 11,465 → 7,803 BP, and farming reaches the Danube Gorges ~8,000 BP**, so the late fraction could carry farmer ancestry — which would make a *reference* share drift with a *source*. **It does not.** f4(Mbuti, Turkey_N; IG_early, IG_late) Z = **−1.67** (n.s.); and `IronGates_LATE ← IronGates_EARLY + Turkey_N` returns Turkey_N at **−4.7%, Z = −2.95** — *significantly negative*. **The late fraction is further from farmers, not closer.**

**A metadata flag is a hypothesis, not a verdict.** Two of the four bins we found turned out to be harmless. **The audit generates candidates; it does not adjudicate, and a paper that reported only the damaging cases would be doing the thing it criticises.**

---

## 3. Recommendations

**For analysts**

1. **Before using any archive population, print its Locality, Political Entity, date range, publication and platform — per individual.** One `cut -f`. In the project described here that command would have saved six drafts, and it was never run.
2. **A pool is a hypothesis. Test it.** qpWave rank on sub-pools defined by **published site and period** — never by genetic similarity, which is circular. A pool needing rank ≥ 1 is not a population, and nothing can descend from it.
3. **Name pools by what they are, not by what you selected for.** We narrowed `Turkey_EBA` by a Levantine-ancestry screen and named the result `_WEST` — a claim about geography from a screen on ancestry. **The name then did the reasoning.**
4. **A period token in a group name is not a date.** Check the dates.
5. **Verify every ground-truth control as a control before trusting it as one.**

**For archives**

6. **Publish a group-ID coherence table with each release** — per group: locality count, geographic and temporal span, publications, platforms. **This is not a new request.** Eisenmann et al. (2018) asked that every *study* publish a per-sample table of site, date and cluster assignment. We are asking that the *archive* publish the aggregate. **AADR already holds every field required. It is a `group_by`.**
7. **Flag group IDs whose modern political label misleads about the ancient region.** `Turkey_MLBA` → `Levant_Alalakh_MLBA`, or minimally an alias.
8. **Normalise Locality strings.** One site, one string.
9. **Adopt a nomenclature whose primary key is not a modern nation state.** Eisenmann et al. (2018) proposed three, and their geographic-temporal system was favoured by the archaeologists among their own authors. **The proposal is seven years old and two of its authors maintain the AADR.**

---

## 4. Limitations

- **The audit measures incoherence, not misnaming, and misnaming is the worse failure** (§2.4). **All prevalence figures are lower bounds.**
- **The withdrawn period check needs region-specific chronologies** before it can be used. As it stands it measures the analyst.
- **Distance is a poor proxy for population discontinuity.** Two sites 300 km apart may be one population; two 50 km apart across a mountain range may not be. **The metrics surface candidates for human reading; they do not adjudicate** — as §2.6 demonstrates.
- **The audit's thresholds create a second blind spot.** Flags fire at ≥ 3 localities, ≥ 3 publications, ≥ 800 years. **A two-site, two-publication, 400-year bin scores zero.** `Greece_Crete_LBA` (§2.4) is exactly that. Lowering the thresholds would raise prevalence and lower precision, and we have not characterised the trade-off. **Prevalence is a lower bound for this reason as well as for the misnaming reason.**
- **§2.5 establishes that `Greece_Minoan_POOL` is not one population. It does not establish why.** Site and date are perfectly confounded in these data, and no design on this sample can separate them.
- **One archive, one release.** We do not know whether other compendia do better.
- **We have not quantified how often mislabelled groups appear in the published qpAdm literature.** That is the follow-up, and a larger paper.
- The analyst is not an academic and this work has not been peer reviewed.

---

## 5. Conclusion

The problem was named in 2018 — by a paper whose authors include the maintainers of the archive, the authors of the Minoan study, and the institute that produced the genomes we mislabelled. **They gave the field a criterion (coherence), a warning (country names are not populations), and a remedy (publish the sample table). Nobody measured whether any of it was being followed.**

The disconfirming evidence was distributed with the archive, one column from the one we used. Every statistical safeguard in a careful pipeline — pre-registration, basis auditing, weight and Z checks, reference screens — **sits downstream of a string that nobody reads.**

And when we finally built the instrument to catch it, **it could not see the case that motivated it**, and it turned up a bin of our own making that was worse than the archive's.

> **Before the panel, before the outgroups, before the SNP basis, before the pre-registration: read the label.**

---

## Data and code availability

**Scripts, logs and manifests: https://doi.org/10.5281/zenodo.21333493** (archived from https://github.com/alexandroskaralakis/aadr-group-ids).

Every number in this paper is emitted by a versioned script in `scripts/` into a log in `logs/`. `manifest/` carries the `extract_f2` population list and SNP basis for each run. `manifest/RETRACTED.md` documents the analyses that were withdrawn during this work, and why — including two whose folds were computed against inadmissible baselines, and one that polluted its own SNP basis.

**Genotype data are public** (AADR v66.p1, Harvard Dataverse, doi:10.7910/DVN/FFIDCW) and are not redistributed here.

**The target genome in the case study is the author's own.** It is not deposited publicly, because it carries information about living relatives who have not consented to its release. It is available to researchers on request.

**§2.5 is fully reproducible from public data.** Its panels are built from the AADR v66.p1 release alone by `build_minoan_panel.sh` (via EIGENSOFT `convertf`) and tested by `run_minoan_split_clean.R`. **Neither script touches the target genome**, and both are scripted end to end.

**The reproducibility gap is now confined to §2.6.** The Kumtepe control targets the author's genome and cannot be reproduced without it. The Iron Gates control runs on the older working panels (`postBA_west`, `postBA_arb`), which were assembled interactively with `convertf`/`mergeit` before this repository existed; **their construction is not scripted.** Both are negative controls that *failed*, so no conclusion in this paper rests on them — but in a paper about unexamined provenance, the asymmetry is stated here rather than left for a reader to find.

## References

- Eisenmann S, Bánffy E, van Dommelen P, Hofmann KP, Maran J, Lazaridis I, Mittnik A, McCormick M, Krause J, Reich D, Stockhammer PW (2018). Reconciling material cultures in archaeology with genetic data: the nomenclature of clusters emerging from archaeogenomic analysis. *Scientific Reports* 8: 13003. doi:10.1038/s41598-018-31123-z
- Mallick S, Micco A, Mah M, Ringbauer H, Lazaridis I, Olalde I, Patterson N, Reich D (2024). The Allen Ancient DNA Resource (AADR): a curated compendium of ancient human genomes. *Scientific Data* 11: 182.
- Skourtanioti E, Krause J et al. (2020). Genomic history of Neolithic to Bronze Age Anatolia, Northern Levant, and Southern Caucasus. *Cell* 181: 1158–1175.
- Ingman T, Stockhammer PW et al. (2021). Human mobility at Tell Atchana (Alalakh). *PLoS ONE* 16: e0241883.
- Koptekin D, Somel M et al. (2023). Spatial and temporal heterogeneity in human mobility in Anatolia. *Current Biology* 33: 41–57.
- Lazaridis I et al. (2022). The genetic history of the Southern Arc. *Science* 377: eabm4247.
- Lazaridis I et al. (2017). Genetic origins of the Minoans and Mycenaeans. *Nature* 548: 214–218.
- Harney É, Patterson N, Reich D, Wakeley J (2021). Assessing the performance of qpAdm. *Genetics* 217(4): iyaa045.
- Flegontova O, Işıldak U, Yüncü E, Williams MP, Huber CD, Vyazov LA, Changmai P, Flegontov P (2025). Performance of qpAdm-based screens for genetic admixture on admixture-graph-shaped histories and stepping-stone landscapes. *Genetics* 230(1): iyaf047.
