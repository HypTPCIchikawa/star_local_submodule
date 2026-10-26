# Gaussian-only nuclear PID purity study

Run from the star-analyzer repository root, in order:

```
python3 tools/nuclear_pid_smooth_20261021/gaussian_only/run.py prepare
python3 tools/nuclear_pid_smooth_20261021/gaussian_only/run.py fit
python3 tools/nuclear_pid_smooth_20261021/gaussian_only/analyze.py
python3 tools/nuclear_pid_smooth_20261021/gaussian_only/render.py
python3 tools/nuclear_pid_smooth_20261021/gaussian_only/finish.py
```

Settings are explicit in settings.yaml. Preparation refuses an existing output directory. Fit resumes only completed runs and refuses incomplete directories. The engines are frozen historical study snapshots; production sources are untouched. Background template columns are absent, not merely assigned a small normalization. All four nuclear Gaussians remain in the contamination calculation.

Linear-fit centers, widths and nuclear yields initialize each fit. The original identity reference is retained. The He4 target mean is fixed to the linear-fit mean only in the configured low-R ALL interval; sigma and nuclear yields remain free. A fixed mean has zero conditional fit error, excluding calibration uncertainty. Missing linear yields use the explicit documented count/4 fallback; invalid historical responses remain unassessed.

The renderer uses the configured mainconf with script/singularity_root_batch.sh. Every native-R bin, including unassessed fits, appears in an atlas. Purity plots use the historical support criteria for this Gaussian-only model. Finish verifies original-input hashes, exact seeds and fixed means, Gaussian-only component sums, He3 half-height coverage, PDF page counts and strict full readback. It merges all-species and He4-only volumes with ROOT CropBox preserved and generates representative PNGs.

See Todo/2026/10/21/he4_r061_20261021.md for results and output links. No event processing, farm submission or production PID change.
