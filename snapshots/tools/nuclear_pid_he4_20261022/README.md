# He4-only PID and regional purity study, 2026-10-22

Explicit settings are in settings.yaml. From the repository root, run these Python3 commands in order:

```
python3 tools/nuclear_pid_he4_20261022/run.py prepare
python3 tools/nuclear_pid_he4_20261022/run.py baseline
python3 tools/nuclear_pid_he4_20261022/calibrate.py
python3 tools/nuclear_pid_he4_20261022/run.py candidate
python3 tools/nuclear_pid_he4_20261022/analyze.py
python3 tools/nuclear_pid_he4_20261022/validate.py
python3 tools/nuclear_pid_he4_20261022/render.py
python3 tools/nuclear_pid_he4_20261022/build_report.py
python3 tools/nuclear_pid_he4_20261022/finish.py
```

Preparation refuses an existing output directory. Fit resumes only completed regions and refuses incomplete output. Mathematical sources are copied from frozen historical study snapshots, not current production sources. A reorthogonalized QR solve calibrates only He4 coefficients without a host BLAS dependency. The original Git header is preserved and a complete candidate header changes only the He4 mean/sigma rows. The strict He4 lower R limit is specified separately in pid_selection.yaml; no production header is changed.

Low ALL bins [0.61,0.80) use four nuclear Gaussians plus linear, [0.80,2.00) use four nuclear Gaussians without smooth background. TOF bins [2.00,5.00) use He4 Gaussian plus linear only. Disabled high-R nuclear slots have exactly zero amplitudes and predictions and are not fitted degrees of freedom. The original four-species TPC best-species gate is preserved, with d/t/He3 exactly Git. He4 response center/width are free in every region.

Validation checks 556 fit-independent raw gate counts, 417 He3 half-height nodes, high-R excluded nuclei, middle-R zero background, old-input hashes and shape-likelihood invariance to gate changes. ROOT5 macros run through the configured mainconf and STAR wrapper; PDF is built with XeLaTeX. Finish checks all PDFs and exports a final manifest. The input is an existing 2021 P24iy partial merged sample, not a new full-statistics event run.

Results and limitations: Todo/2026/10/22/report_he4_20261022.md. No job submission, hadd, production PID adoption, commit or push.
