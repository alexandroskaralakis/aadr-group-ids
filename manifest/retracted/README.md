# Retracted manifests

These files are the output of analyses that were **withdrawn**. They are kept
because the retraction sequence is part of the papers' argument, not because
their numbers are usable. **No figure in either paper is taken from them.**

- `run_chisq_calibration_*` — the superseded §3 calibration. Computed on a panel
  into which the author's genome had been merged (so the SNP basis was capped by
  a consumer genotyping chip), and using `Greece_Minoan_POOL` as a source — a
  pool that fails qpWave rank-0 at p = 3.2e-08 and is therefore not a population.
  Its premise, that chi^2 scales with SNP count, is **false**: under random SNP
  loss chi^2 *rises*. See `manifest/RETRACTED.md`.

  Superseded by `run_chisq_null.R` + `run_chisq_exhibit.R`, which are built
  entirely from public AADR data.
