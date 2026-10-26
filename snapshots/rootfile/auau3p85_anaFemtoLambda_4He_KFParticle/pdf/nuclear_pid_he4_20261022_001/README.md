# He4-only PID / regional purity: 2026-10-22

- [Illustrated Japanese report, 10 pages](he4_pid_purity_report_20261022.pdf)
- [All 199 He4 fit intervals](fit_all_regions.pdf)
- [Purity, full and low-R views](purity_vs_p_over_q.pdf)
- [Low R Gaussian+linear, 19 pages](fit_low_linear.pdf)
- [Middle R Gaussian nuclei only, 120 pages](fit_low_none.pdf)
- [High R He4 Gaussian+linear only, 60 pages](fit_high_linear.pdf)
- [Candidate PID header](candidate_pid/NuclearIdDeDxVsMom.h), [coefficient diff](candidate_pid/he4_coefficients.patch), [selection YAML](pid_selection.yaml)

All d/t/He3 PID coefficients and validity domains are exactly Git 8f7d539. Only He4 polynomial coefficients change. The R>0.61 cut is separate from the candidate header. No production files are changed. High-R exclusion of other nuclei is the specified model assumption; original four-species TPC competition remains.

76/88 aggregate intervals pass historical support (4/4, 24/24, 48/60). The first 0.61–0.65 bin has P=63.65813 +/- 0.06628 percent, conditional yield uncertainty only. Unsupported intervals remain in the atlas; no interpolation, veto or new single upper cutoff.

Reproducible tools: tools/nuclear_pid_he4_20261022/README.md in the repository. Full validation is in completed.json: 398 baseline/candidate fits, 556 independent raw counts, 417 He3 half-height nodes, unchanged 102 old inputs and all 415 PDF pages read back. Dataset is an existing 2021 P24iy partial merged sample. See input_provenance.json, coefficients.json, native.tsv and summary.json.
