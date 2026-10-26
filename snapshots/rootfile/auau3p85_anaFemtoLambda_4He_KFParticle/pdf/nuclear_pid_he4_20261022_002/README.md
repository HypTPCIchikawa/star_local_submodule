# He4 low-R Gaussian-only / fixed-center revision (002)

- [Illustrated revised report, 11 pages](he4_pid_purity_report_20261022.pdf)
- [New low-R fixed-center fits, 19 pages](fit_low_fixed.pdf)
- [All 199 He4 intervals](fit_all_regions.pdf)
- [Combined purity, 2 pages](purity_vs_p_over_q.pdf)
- [Fixed-center validation](fixed_center_validation.json), [old/new low-R comparison](low_comparison.json), [final manifest](completed.json)

Only R<0.8 is refitted, with He4 mean fixed to preceding Gaussian+linear fit values and no smooth background. Width and nuclear yields are free. All four nuclear Gaussians remain in low-R contamination. R>=0.8 (180 native intervals, 84 aggregates) and PID coefficients are exactly preserved. The mean is removed from the optimization vector; zero mean error is conditional, not absence of calibration uncertainty.

The first aggregate bin has calculated P=99.6185% but fails historical width-ratio support. It is not plotted as a supported purity point, and its fit remains in the atlas. 75/88 aggregates are supported. Near-unity purity after removing background is not proof of physical purity or fit improvement.

All PDFs (416 pages including combined copies) passed readback. Source, settings, commands, raw/count/half-height/parameter validation and unchanged-input evidence are saved. Original 001 outputs and production settings are preserved. Reproduce with tools/nuclear_pid_he4_fixedmean_20261022/README.md in the repository.

## TOFマップ追加版

reportは12頁。6頁に既存TAG_D（dと共通m²窓）のR–dE/dxと確定PID中心・±2σを追加。全R表示とR≥2拡大を掲載。Purity・係数・fit結果は不変更。旧11頁reportは解析rootfile配下tests/report_tof_map_20261022_001/へ保持。更新記録：report_tof_update.json。
