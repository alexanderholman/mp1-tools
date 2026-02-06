<!-- filename: Coder.md -->

## [SPEC] Specialism Addendum — Coder (v1.0)

### Purpose
Defines standards for writing usable code artifacts.

### Operating Rules
- Functions named by intent (verbs).
- Inputs/outputs must be explicit.
- Include a minimal "How to run" section.
- Prefer fewer dependencies.
- Handle basic failure modes or document them.

### Outputs (typical)
- scripts (python/bash)
- config files
- parsers/validators
- small CLIs
- notebooks (only if requested)

### Quality Gates (Coder)
- Runs end-to-end from clean checkout (or run steps provided)
- Dependencies listed
- Basic error handling or documented failure modes
- No hidden environment assumptions

### [TEST] Acceptance Checks
- A user can copy/paste the file and execute it with stated prerequisites
- At least one example invocation exists
