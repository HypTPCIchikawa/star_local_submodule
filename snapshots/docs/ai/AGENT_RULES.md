# Agent rules

This document is the source of truth for project-specific agent workflow rules.
Framework-level design principles remain in `../../PHILOSOPHY.md`.

## Editing and implementation rules

- **Approval infrastructure failures (user decision, 2026-10-19)**: For auto-review capacity/timeouts or sandbox namespace failures, follow [approval-recovery](skills/approval-recovery.md). Promptly report unexecuted work and offer the manual-approval workflow instead of repeatedly stalling. Do not change your own approval policy or bypass a denial. Continue unaffected authorized work.

- **File Modification Approval**: The agent MUST NOT modify, create, or delete any source, configuration, or documentation files in the repository without explicit prior approval from the user.
- **Macros**: Do not implement custom event loops. Build `StChain`, run `Init()`, `Make(i)`, and `Finish()`. The Maker-using macro must be ACLiC-compiled and linked through `AddLinkedLibs`.
- **Makers**: Implement `Init()`, `Make()`, `Clear()`, and `Finish()`. Keep histogram lifecycle in `DeclareHistograms()` and `WriteHistograms()` / `Finish()`.
- **Config**: Prefer YAML updates over recompilation when changing cuts and histogram values. Main entry is `config/mainconf/main_<anaName>.yaml`.
- **No hardcoded analysis parameters**: Put thresholds and cuts in YAML configs.
- **Comments explain the current code to a first-time reader.** Concise, careful English: what the code does, why,
  units, assumptions and edge cases. No dates, no "decided by", no change logs ("was X until ..."), no meeting notes and
  no pointers to personal or volatile files (`local/`, `scratch/`, analysis notes) in code, configs or READMEs. *Why:*
  such text ages immediately and hides what the code does now; the history is in the git log (commit messages) and in
  the analysis notes, where it belongs.
- **Each document has one job; write it in the present tense.** `PHILOSOPHY.md`: framework principles with their reasons.
  `docs/ai/AGENT_RULES.md`: rules for this repository with their reasons. `docs/ai/skills/*.md`: step-by-step procedures.
  A component's `README.md`: what it does, how to use it, its configuration and what has been verified, completely and
  as it is now. `docs/REFERENCE.md`: commands and entry points. Decisions, studies and their history: analysis notes and
  commit messages, not these files.
- **analysis_info sync**: Keep `anaName`, macro base names, output naming, and `starTag` consistent.
- **Todo updates (user decision, 2026-10-13)**: When a task has a record under `Todo/`, keep its Markdown plan/report and directory index current as analysis results, code, or YAML settings change. Record what is implemented versus only proposed, links to current plots and settings, validation results, limitations, and remaining work. Preserve historical results and mark superseded conclusions rather than silently erasing them; update the task record during the work, not only at the final handoff.
- **Build environment for debug/repro**: For debug work and any build used for farm submission, use a **batch-matched STAR toolchain** before submit. **Either** (a) interactive SL7: `sl7`, then `source ./script/setup.sh <mainconf>` and `make`, **or** (b) from any host with Singularity: **`./script/singularity_make.sh <mainconf>`** — option (b) runs `make` inside `star-bnl/star-sw:latest` with the same `sl73_*` / `sl74_*` `STAR_HOST_SYS` resolution as batch and **counts as the SL7-class build** for agents and reproducibility. Do not rely on host-only `make` for farm-bound libraries or heap/exit-crash diagnosis. For `root4star` / QA on fragile login nodes, use the other **`script/singularity_*`** wrappers in `../../docs/REFERENCE.md`.

## Batch safety: how values reach a job, and what happens when they do not

These rules exist because three defects of the same shape reached a 241-job farm pilot on
2026-09-15, each invisible to every local test
(`analysisnote/auau3p85fxt_anaFemtoPhi/root-tree-study/results/pilot-farm-20260915.md`).

- **A value the analysis needs arrives as an argument, not through the environment.** The
  environment may *override* for local convenience; it may not *supply*. Locally the wrapper
  scripts export `STAR_ANA_*`; the batch joblist calls `root4star` directly and exports nothing, so
  anything read only from the environment is empty on the farm — and empty was silently accepted
  for the job id, for the mainconf path, and for every tree schema key in turn.
- **Never continue with a default when a configuration cannot be read.** `if (missing) return
  kTRUE;` and "fall back to `.current_mainconf`" both turn a broken job into a job that produces
  plausible, wrong output. A job that cannot find its configuration must fail.
- **No built-in default mainconf.** A macro that is given no configuration stops.
- **Resolve relative paths against `gSystem->WorkingDirectory()`, never `gSystem->Getenv("PWD")`.**
  A SUMS job is a `csh` script that `cd`s into its runtime bundle, and csh does not update `PWD`.
- **Print the resolved configuration path and where it came from, unconditionally**, so that a job
  running with the wrong configuration is visible in the first lines of its log rather than in the
  output months later.
- **Required keys are checked, not defaulted.** A YAML read of the form "use the key if present,
  else keep the default" makes a typo silent. Keys that decide what an output *contains* (schema
  version, which species are stored, the storage envelope) must be present or the maker fails.
- **Provenance that was not supplied is recorded as `unrecorded`, never as an empty string,** and a
  warning is printed. An empty `gitRevision` reads like "no revision", which is what the pilot
  trees carried.

## Adding analyses and scripts

- **New analysis** requires: Maker code under `StMaker/StXXXMaker/` (auto-built; no Makefile edit), `run_anaXxx.C` that Loads/links `libStCommon` before the Maker `.so`, `anaXxx.C`, `script/run_anaXxx.sh`, mainconf, referenced YAMLs, and analysis_info.
- **Joblist generation**: Use `script/generate_joblist.sh <mainconf>` and produce `job/joblist/joblist_<anaName>.xml` (`analysis.anaName`); the batch runtime mainconf is the one passed to `generate_joblist.sh`. For input files you already have (a PicoDst list of your own, the candidate stores of a pair stage) use `script/generate_joblist_from_list.py` (same template and analysis info; input list or directory, files per subjob, output directories and a label are arguments; writes `joblist_<anaName>_<label>.xml`).
- **Skill sync after docs edits**: If a task adds/updates `docs/ai/skills/*.md`, run `script/sync_and_check_skills.sh` before finishing.
- **Scripts and workflow docs**: When script or run flow changes, update:
  - `../../docs/REFERENCE.md`
  - `../../INSTALL.md` if onboarding changed
  - `../../README.md` when index/table changes
  - `../../job/run/README.md` for submit/cleanup/archive details

## StRoot vendoring

- **When to vendor:** only when `$STAR/lib` is insufficient and the code cannot live in `StMaker/`; see [`docs/ai/skills/reuse-star-stroot.md`](skills/reuse-star-stroot.md).
- **Canonical example:** [`StRoot/StRefMultCorr/`](../../StRoot/StRefMultCorr/) with `PROVENANCE.md`, `Makefile` target, and repo-local `.so` load in `run_ana*.C` (not `$STAR/lib` duplicate).
- **Agent procedure:** follow `reuse-star-stroot.md` when searching `$STAR/StRoot` or past trees and copying into `StRoot/<Package>/`.

## 3.85 GeV dataset policy (user decision, 2026-10-01 JST)

- For this user's future 3.85 GeV Lambda / Lambda-nucleus analyses and studies, use
  **2021 P24iy only** (`production_3p85GeV_fixedTarget_2021`). Do not substitute
  2019 P24iy or another year/production when access fails; report the problem.
- Preserve historical ROOT outputs, submitted XML/config snapshots and study records
  unchanged. They are not approved inputs or production templates for future work.
  Never merge the 2019 results into the 2021 sample.
- Check the actual catalog URL and explicit local input list, not just `nFiles: all`
  or the analysis name. The old 2019 validation list is historical only.
- Dataset selection is not evidence of complete statistics: distinguish catalog event
  counts from measured PicoDst Entries, and verify coverage before claiming a full run.

## Centrality work

- **Bin conventions:** [`StRoot/StRefMultCorr/README.md`](../../StRoot/StRefMultCorr/README.md) is canonical (cent9: larger index = more central; 0–60% → cent9 2–8).
- **Agent procedure:** follow [`docs/ai/skills/centrality-strefmultcorr.md`](skills/centrality-strefmultcorr.md) when changing centrality YAML, `CentralityHelper`, cent axes, or `cfCent9Min`/`Max`.

## Femto architecture (reconstruction parts vs femto pairing)

- **Reconstruction is a part, femto is the pairing layer.** A particle reconstruction (KFParticle,
  helix/V0, track PID, ...) is a part that produces candidates; one femto Maker combines them into
  SE/ME pairs. The femto Maker must not depend on how a candidate was reconstructed.
- **KFParticle is not special in the femto layer.** No femto Maker links a KFParticle Maker or library,
  and there is no KFParticle-specific femto build rule. KFParticle code lives in the parts
  (`StMaker/kfparticle/`, `*KFParticleMaker` reconstruction Makers), which the runner loads.
- **Two femto Makers exist.** `StFemtoMaker` (shared; φ–p, φ–d, Λ–nucleus; builders and a Λ provider)
  and `StCandidateFemtoMaker` (parts given by the macro through `include/CandidateSource.h`).
  Extend one of them with a new species/channel (YAML + a source or builder); do not add a new
  dedicated pairing Maker per channel. Do not change `StFemtoMaker` for a channel of the other Maker.
- **Several particles from one KFParticle pass** come from one part (`KfParticleCandidateSource`): add
  the particle there, not a second Finder pass.
- **CF, sidebands, purity corrections and BDT thresholds are chosen outside the Maker**
  (`common/macro/candidate_femto/`), from the pair THnSparse; the Maker fills, it does not decide.
- **Read the PicoDst once (`StCandidateFemtoMaker`).** Production runs in two stages with the same macro: the
  reconstruction stage (`writeCandidateStore: true`) runs the event selection and the parts and stores every candidate
  with its event record (`CandidateStore.h`); the pair stage gets the store files instead of PicoDsts and does BDT
  scoring, pairs, mixing and training trees. A part's candidates must be fully described by what it gives through
  `CandidateSource` (no PicoDst access in the pair stage); the pair stage must reproduce a one-pass run bit by bit (check it
  when a part or the store changes). Mixing pools live inside one job: more events per pair-stage job means fuller pools.
- **Special methods stay out of the common pairing code.** A study that needs another mixing rule, extra axes or a skim
  is prototyped separately and compared with the standard result; it enters the Maker only if it changes the answer.
- Details: [`StMaker/StCandidateFemtoMaker/README.md`](../../StMaker/StCandidateFemtoMaker/README.md),
  [`StMaker/StFemtoMaker/README.md`](../../StMaker/StFemtoMaker/README.md).

## Candidate BDT (machine-learned candidate selections)

These rules keep a BDT reproducible and keep training and application from drifting apart. The
contract below is enforced in code (`StMaker/StCandidateFemtoMaker/CandidateBdtScorer.h`, `CandidateTreeEnsemble.h`,
`common/macro/datadriven_bdt/`); do not bypass a check to make a model load. Procedure:
[`docs/ai/skills/candidate-bdt.md`](skills/candidate-bdt.md). Design and how to add a BDT for a new particle (what is generic,
what the analyser decides): [`docs/ai/candidate-bdt-guide.md`](candidate-bdt-guide.md).

- **What a BDT is for.** Candidates that form a mass peak (decays reconstructed with KFParticle). The labels come from
  the data (mass window vs side bands; no Monte Carlo), so the selection is defined by the analysis sample itself and needs
  no simulation tuning. Identified tracks (K⁺) are selected with PID cuts; a BDT for them only if their purity is not
  good enough.
- **A score belongs to a candidate, a threshold to the analysis.** One model per species (not per
  channel): the same candidate has the same score in every channel. The femto Maker computes scores
  and stores them on the score axes of the pair THnSparse; thresholds are chosen in the macros
  (`species_<s>_scoreMin`), never in the Maker. Score axis edges (`<s>_scoreEdges`) must contain every
  threshold you intend to study, so changing a threshold needs no rerun.
- **Training is outside the chain, with one tool.** Train with `common/macro/datadriven_bdt/` (Python + xgboost, one
  YAML recipe per model; data only: mass window vs side bands) on the candidate trees the Maker writes
  (`writeCandidateTrees: true`). Do not train inside a Maker and do not write a per-channel training macro; extend the
  toolkit instead (a new way to define signal and background is a new `labels.method`, not a new script).
- **Training/application contract:**
  - *Variables:* model expressions may use only the branches of the candidate trees (`runId`, `eventId`,
    `cent9`, `mass`, `px`, `py`, `pz` and the part's `VariableNames()`); derived quantities are
    expressions (`sqrt(px*px+py*py)`). Train on trees written by the same part version the model is
    applied with; a change of a part that alters a variable means retraining.
  - *Inputs:* the input expressions of a model live only in its tree-ensemble files (`input` records); the model file
    does not repeat them. `inputs.use: all` takes every variable of the part; inputs may not use `mass`, `runId`,
    `eventId` or `score`, and pT stays out of the inputs (the side bands are reweighted in it instead).
  - *Folds:* `runParity` only. Fold f is trained on runs with `runId % 2 == f` and scores the other
    parity, so no candidate is scored by a model that saw its run and no event is split. The Maker refuses
    any other fold scheme.
  - *Domain:* a model carries the event selection and preselection of its training sample; candidates
    outside get `outsideTrainingScore` (below -1, in the first score bin of its own). The domain must not use the
    flag branches `sharedTrackRemoved_<channel>`.
  - *Trees and flags:* `writeCandidateTrees: true` requires `candidateTreeSpecies` (the species that get a tree). Trees keep
    every usable candidate; each tree has one flag per channel of its species, `sharedTrackRemoved_<channel>` (1 = removed
    from that channel by the shared-track rule; each channel decides independently). Flags are never inputs.
  - *Training selection:* `data.trainingSelection` (required, may be `""`) restricts the training samples only. Use it to
    train without the candidates removed in the channel(s) the model serves (for a species in several channels: all of
    them); do not put flags into `eventSelection`/`preselection`, which define the scoring domain. A condition the
    trees cannot express (e.g. "the event contains a partner") cannot be part of a Maker model: state
    it in the model's `particle` text; pairs impose it implicitly, single-species QA does not.
  - *Expressions:* one grammar for both sides (a TTreeFormula subset): `+ - * /`, comparisons,
    `&& || !`, `log10 log exp sqrt abs max min`; parenthesise mixed `&&`/`||`; no `#` (YAML comment). It is
    implemented twice, `CandidateExpression` (C++) and `candidate_expression.py` (line-by-line port): change both
    together and rerun `check`.
  - *Score:* margin = sum of the tree leaves + base margin, score = 2 / (1 + exp(-margin)) - 1, defined by the
    tree-ensemble file (`tree_ensemble.py` is the reference; `CandidateTreeEnsemble.h` must agree to 1e-6).
- **A model is a file.** `<name>.model.yaml` (written by `datadriven_bdt` action `model`) plus its two
  tree-ensemble files `<name>_f{0,1}.trees.txt`. Every key is required and records provenance (training data,
  candidate counts, recipe, git revision, learner and its parameters, validation). Never hand-edit a model file; retrain and rewrite it.
  Production models live under `config/bdt/<anaName>/` so batch jobs find them (`config/` is git-ignored). Keep them
  local for day-to-day work; force-add only when publishing that analysis through the shared repository. Study models
  stay with the study outputs (not tracked). A model used by a batch job must be under `config/`: the joblist ships
  `analysis config lib include StMaker StRoot`, nothing else (a test model gets its own clearly named
  `config/bdt/<anaName>_<test>/` directory with a README and is not committed). Keep one catalog of all models per analysis
  (e.g. `bdt_models.md` in the version's `share/figure` directory).
- **Thresholds are chosen by a scan of the stored scores, and the scan is kept.** Scan the thresholds on one production
  output with the CF macro options (`scoreMinA` / `scoreMinB`), compare SE and ME purities, failed purity fits, statistical
  error and the stability of the CF, and choose the loosest thresholds with high purity and a stable CF (procedure:
  `common/macro/candidate_femto/README.md`, "Choosing BDT thresholds"). *Why:* tight cuts reduce the reliance on the purity
  correction but lose pairs and destabilise the per-k* purity fits; loose cuts do the opposite. The chosen values live in a
  `_bdtSelection_<thresholds>` YAML; the scan stays as a study for the systematic uncertainty.
- **Purities are measured where the pairs are.** Purities for the CF come from the pair THnSparse, SE and ME separately
  (ME events are not SE events, so their purities differ); per-species histograms are QA only.
- **Mainconf wiring.** Optional key `candidateBdt` → YAML with `speciesKeys`, `outsideTrainingScore`,
  `<s>_model` (relative to `config/`), `<s>_scoreEdges`. No key = no BDT (one score bin). The Maker
  writes each model's description and weight checksums into `CandidateFemtoConfiguration`.
- **Before a model is used for physics:** (1) `datadriven_bdt check` on trees written with the model
  (Maker score = offline score for every candidate, domain identical); (2) the sculpting checks of the
  `datadriven_bdt` report (ratio at equal purity, ε(m), background-model variation, mass-correlation warnings);
  (3) the overtraining, SHAP and training-size pages of the report (enough training statistics, no input that
  dominates for a non-physical reason); (4) purities quoted from the same sample as the figure. Put the report path
  into the model's `validation` key.
- **Batch environment.** Scoring runs in ROOT 5.34 / C++98 jobs, with no Python at job time: the Maker evaluates the
  tree-ensemble file itself (`CandidateTreeEnsemble.h`). Another model type needs another C++98 evaluator behind
  `CandidateBdtScorer` and a reference implementation in the toolkit, not a change of the pairing code.
- **Extending.** New species: a model plus a `speciesKeys` entry, no code change. New input variable:
  add it to the part's `VariableNames()` (unique names; it then appears in the trees and in
  expressions). `StFemtoMaker` has no BDT; if one is needed there, move `CandidateBdtScorer.h` to a shared
  place and reuse it rather than writing a second scorer.
- **Farm training trees.** The joblist template copies only `<stem>_$JOBID.root`; for a job that writes candidate trees,
  generate the joblist with `script/generate_joblist_from_list.py --tree-output-dir`, which adds the `<stem>_$JOBID.tree.root` output.

## Femto φ–p checkHist / CF

- **CF and sideband logic live in** `common/macro/checkHistAnaFemtoPhiProton.C` (not in Maker). Skill: [`docs/ai/skills/femto-species-naming.md`](skills/femto-species-naming.md).
- **Dual PDF:** `singularity_checkHistAnaFemtoPhiProton.sh` writes QA PDF + `{anaName}_checkHistAnaFemtoPhiProton_CF_{jobid}.pdf` under `share/figure/<anaName>/`.
- **YAML (maker):** `cfCentSlices` (defaults in `FemtoConfig::SetDefaults`), `cfCentSlicesQaPdfInclude`, `cfPdfExcludeQaSlices`, `sidebandSubtractAlpha`, `sidebandAlphaMode`, `negativeBinPolicy`.
- **YAML (mixing):** `mixingMode`, `maxMixedPairsPerEvent`, `mixBothDirections`, `bufferSize` for ME statistics (`StFemtoMaker::FillMixedEventPairs`). Mixing bin = same Vz × **same cent9** × EP; ME is cross-event only.
- **SB-LR:** no new maker channel; checkHist `combineSidebandLR` on projected leftSB + rightSB histograms.
- **Sub CF keys (checkHist cache):** `CF_sig_raw`, `CF_sig_sub_SBL`, `CF_sig_sub_SBR`, `CF_sig_sub_SBLR` per slice id.
- **Plan / ops record:** `analysisnote/YYYYMMDD/femto_cent_sb_cf_plan.md` (when present).

## Naming and terminology

- **anaName**: `{system}_{anaId}[_condition]` and use it consistently across configs, outputs, and job artifacts.
- **jobid**: 32-char hex from star-submit, used in configlog/output naming.
- **mainconf**: Canonical configuration entry point under `config/mainconf/`.
- **baseRunMacro/baseAnaMacro**: Names without `.C`; joblist and ACLiC output names derive from them.

## Directory and path conventions

- `analysis/`: paired `run_anaXxx.C` and `anaXxx.C`
- `StMaker/`: compiled Maker source to `lib/libStXXXMaker.so`
- `config/`: concern-separated YAMLs; mainconf references others relative to `config/`
- `script/`: setup, local run, joblist generation helpers
- `job/run/`: submit and run-artifact management tools
- `lib/`: generated binaries, not committed

## Output conventions

- ROOT output: `rootfile/<anaName>/`
- After verified hadd, retain merged ROOTs at the output directory's top level and
  relocate only the verified subjob inputs to `subjobs/<jobid>/`, with a README and
  relocation manifest. Do not delete ROOT data or move failed/unmerged inputs.
  Use `watch_job_merge_notify.py` (automatic) or `archive_merged_subjobs.py --runmeta ...`
  (already-validated runs); preserve historical provenance and never mix merged and
  archived inputs in another hadd. This is the user's standing output policy (2026-10-01).
- Submitted config snapshot: `job/run/configlog/config_<anaName>_<jobid>.txt`
- Submitted joblist snapshot: `job/run/joblistlog/joblist_<anaName>_<jobid>.xml`
- Submit-time runmeta manifest: `job/run/runmeta/runmeta_<anaName>_<jobid>.json`
- Submit-time runmeta sidecars: `job/run/runmeta/gitstatus_<anaName>_<jobid>.txt`, `gitdiff_<anaName>_<jobid>.patch`, `gitsubmodules_<anaName>_<jobid>.txt`, `runtime_bundle_<anaName>_<jobid>.tar.gz`, `sums_artifacts_<anaName>_<jobid>.tar.gz`, `submit_stdout_<anaName>_<jobid>.txt`
- Watch-merge log: `job/run/watchmerge/watchmerge_<anaName>_<jobid>.pid` / `.log`; completion metadata in runmeta `postProcess.watchMerge`
- QA PDF: `share/figure/<anaName>/<anaName>_checkHistAnaPhi[_<jobid>].pdf`

## Submit reproducibility artifacts

- Treat `jobid` as the key used to join config, submitted XML, runmeta, output ROOT, and QA.
- `submit.sh` should save per-`jobid` reproducibility artifacts immediately after jobid extraction, not as a later manual step.
- `runmeta` is the machine-readable index for submit-time provenance; use it before scanning loose `job/run/<anaName><jobid>*` files.
