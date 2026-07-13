# Group IDs are not populations

Scripts, logs and manifests for two preprints:

- **Paper A** — an audit of AADR group IDs, and a case in which an archive
  label determined a qpAdm result end to end.
- **Paper B** — Test 3 (a source-descendant diagnostic for reference-into-source
  gene flow) and a chi-squared signature of non-random SNP loss.

## Data

Not included. Download AADR v66.p1 (1240K) from Harvard Dataverse:
https://doi.org/10.7910/DVN/FFIDCW
Place the .geno / .snp / .ind / .anno files in the project root.
Cite: Mallick et al. 2024, Sci Data 11:182.

## Reproduce

Run in this order — later scripts read manifests written by earlier ones:

    Rscript scripts/run_uniform_3source.R    > logs/run_uniform_3source.log    2>&1
    Rscript scripts/run_eastern_term_v2.R    > logs/run_eastern_term_v2.log    2>&1
    Rscript scripts/run_test3_v2.R           > logs/run_test3_v2.log           2>&1
    Rscript scripts/run_test3_controls_v3.R  > logs/run_test3_controls_v3.log  2>&1
    Rscript scripts/run_anchor_integrity.R   > logs/run_anchor_integrity.log   2>&1
    Rscript scripts/run_aadr_audit.R         > logs/run_aadr_audit.log         2>&1

Every number in both papers is emitted by one of these into `logs/`.
`manifest/` carries the population list and SNP basis for each script.

## Known gap

The working panels (`postBA_west`, `postBA_arb`) were built interactively with
EIGENSOFT `convertf`/`mergeit` before this repository existed. That build is
NOT scripted. All downstream analysis is reproducible given the panels; the
panels themselves are not. This is stated in both papers.

## Environment

admixtools2 v2.0.10 | R 4.6.0 | EIGENSOFT 8600 | PLINK v1.90 | macOS
