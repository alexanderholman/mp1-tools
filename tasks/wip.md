# WIP

## Active
- Phase 3 planning active: migrate `artemis_utilities` into `mlp_tools` package with safer execution controls.

## Next
- Add `mlp-tools` Artemis CLI wrapper with path overrides and `--dry-run`.
- Resolve template rendering mismatch in Artemis SBATCH setup as a separate tested change.
- Add docs/contracts for VASP and Artemis command behavior in `mp1-tools` docs.

## Done
- Phase 1 complete: package scaffold, `mp1-id`, `mp1-energies`, tests passing.
- Phase 2 complete: `mp1-vasp` added with explicit POTCAR directory configuration and loud failure mode.
- Namespace standardized to `mlp_tools` (from `mp1_tools`), tests passing (`5 passed`).
