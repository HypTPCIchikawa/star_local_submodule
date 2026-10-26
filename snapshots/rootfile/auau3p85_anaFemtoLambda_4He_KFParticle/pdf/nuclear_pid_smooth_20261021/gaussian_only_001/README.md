# Gaussian-only nuclear PID fit atlas

All four nuclear Gaussian components; no smooth background. ALL below R=2, existing TOF m² gates above. He4 uses the d TOF gate and strict R>0.61. He4 means at 0.61<R<0.70 are fixed to the preceding Gaussian+linear fit; widths and yields remain free.

- [All 1,579 intervals](gaussian_all_intervals.pdf)
- [All 199 He4 intervals](fit_gaussian_he4_all.pdf)
- [Updated purity](purity_vs_p_over_q.pdf)

| Species | ALL PDF | TOF PDF | All-volume pages (ALL / TOF) |
|---|---|---|---|
| d | [160 pages](fit_low_d.pdf) | [300 pages](fit_high_d.pdf) | 1–160 / 620–919 |
| t | [160 pages](fit_low_t.pdf) | [300 pages](fit_high_t.pdf) | 161–320 / 920–1219 |
| He3 | [160 pages](fit_low_he3.pdf) | [300 pages](fit_high_he3.pdf) | 321–480 / 1220–1519 |
| He4 | [139 pages](fit_low_he4.pdf) | [60 pages](fit_high_he4.pdf) | 481–619 / 1520–1579 |

One R bin per page, linear/log y scales. `valid` is numerical validity, `quality` is historical support, `mu fixed` identifies conditional fixed means. Unassessed intervals remain in the atlas. 175/364 aggregate bins pass support; near-unity purity under this background-free assumption is not evidence that the removed background is absent. The first He4 aggregate bin fails the existing width-ratio condition despite its high computed purity.

See settings.yaml, summary.json, native.tsv, fixed_center_validation.json, warm_seed_validation.json, halfheight_validation.json and completed.json. All 11 PDFs (3,359 pages including merged copies) passed readback, all ROOT spectra have zero smooth background, identical raw data and total equal to the nuclear sum. Source snapshots and commands are saved. merge_visual_validation.json records 16 source/merged boundary-page comparisons. Historical inputs and production PID are preserved.
