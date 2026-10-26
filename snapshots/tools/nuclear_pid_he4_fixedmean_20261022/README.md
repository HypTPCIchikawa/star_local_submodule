# He4 low-R fixed-center Gaussian-only revision

This revision changes only ALL R<0.8. The He4 Gaussian center is fixed per native R bin to the previous 2026-10-22 Gaussian+linear fit. Width and nuclear yields remain free. Smooth background columns are absent. The same four nuclear Gaussians are retained for contamination. Previous PID coefficients and R>=0.8 outputs are copied exactly; coefficient calibration is not rerun.

From the repository root, run in order:

```
python3 tools/nuclear_pid_he4_fixedmean_20261022/run.py prepare
python3 tools/nuclear_pid_he4_fixedmean_20261022/run.py fit
python3 tools/nuclear_pid_he4_fixedmean_20261022/analyze.py
python3 tools/nuclear_pid_he4_fixedmean_20261022/validate.py
python3 tools/nuclear_pid_he4_fixedmean_20261022/render.py
python3 tools/nuclear_pid_he4_fixedmean_20261022/build_report.py
python3 tools/nuclear_pid_he4_fixedmean_20261022/finish.py
```

Settings select previous output, fixed centers and model boundaries explicitly. The baseline free-fit identity reference is preserved when setting starts. Mean is removed from the optimizer vector, its conditional error is zero, width stays free and NDF counts the remaining free parameters. This excludes uncertainty in the mean estimated from the same data, not a claim of exact physical calibration. Historical support criteria are retained.

The 19 new fits are combined with 180 unchanged native intervals. Validation checks exact fixed centers, the entire unchanged region JSONs/ROOT products, all PID coefficients, independent raw counts, He3 half-height coverage, zero background below R=2 and no other high-R nuclei. ROOT5 runs via the configured STAR wrapper; the illustrated Japanese report uses XeLaTeX. Original outputs are preserved in 001 and the revision uses 002. No event processing, production changes, farm submission, hadd, commit or push.
