# STAR analysis agent entrypoint

This file is the entry point for agent tools that follow the `AGENTS.md` convention (Codex and
others). `CLAUDE.md` is the same entry point for Claude Code, and both point at the same
source-of-truth documents — the rules live in `docs/ai/`, not here.

## Required reading order

1. `PHILOSOPHY.md`
2. `docs/ai/AGENT_RULES.md`
3. `docs/ai/skills/` — choose the procedure file that matches the task

## Intent

- Keep the StChain / StMaker architecture and the YAML-driven reproducibility constraints intact.
- Use `mainconf` as the operational entry point for setup, run and job tooling.
- Treat `docs/ai/*` and `PHILOSOPHY.md` as source-of-truth. `.cursor/*` is a wrapper for Cursor
  integration and is generated from `docs/ai/skills/`; do not edit it by hand.

## Rules that are easy to get wrong

These are in `docs/ai/AGENT_RULES.md` in full. They are repeated here because each one was written
after a defect reached the farm.

- **Do not modify, create or delete repository files without explicit prior approval.**
- **A value the analysis needs arrives as an argument, not through the environment.** The
  environment may override for local convenience; it may not supply. A batch job exports nothing.
- **Never continue with a default when a configuration cannot be read.** A job that cannot find its
  configuration must fail, loudly. A silent default produces plausible, wrong output.
- **No hardcoded analysis parameters.** Thresholds and cuts live in YAML, and anything that changes
  a result has to be visible in the configuration and in the provenance.
- **Build with the batch-matched toolchain** before anything farm-bound:
  `./script/singularity_make.sh <mainconf>`. A host-only `make` does not count.
- **Rebuild the libraries after touching a config struct.** Adding a member to `FemtoConfig` and
  running against the old `libStarAnaConfig.so` is a header/library layout mismatch, and it
  segfaults somewhere unrelated.

## Skill sync completion rule

For approval-review capacity/timeouts or sandbox startup failures, follow
`docs/ai/skills/approval-recovery.md`. Promptly notify the user when manual approval is needed;
do not silently keep retrying or report unexecuted edits as saved.

- Source-of-truth for skills is `docs/ai/skills/*.md`.
- After adding or editing any file under `docs/ai/skills/`, always run:
  - `script/sync_and_check_skills.sh`
- The task is not complete until that command passes.

## Human-facing references

- `README.md`
- `INSTALL.md`
- `docs/REFERENCE.md`
- `job/run/README.md`

## A note on `.agents`

There used to be a `.agents` symlink pointing at `.cursor`. Nothing in the repository read it, it
was gitignored, and the Codex sandbox rejected it. It was removed on 2026-09-17 and replaced by
this file, which is what the convention actually asks for. Do not recreate the symlink; add the
entry point a tool needs as a real file instead.
