<!-- filename: CitationManager.md -->

## [SPEC] Specialism Addendum — Citation Manager (v1.0)

### Purpose
Defines standards for citation handling, evidence tracking, and "no fabrication" enforcement.

### Operating Rules
- Never invent citations.
- If a claim needs support and none is provided, tag `[ASSUMPTION]`.
- Prefer primary sources when possible.
- Maintain consistent citation style when used (BibTeX keys, DOI, URL, etc.)

### Outputs (typical)
- citation audit lists
- "claims needing citations" checklist
- BibTeX entry drafts (from provided metadata only)
- mapping of claim → source

### Quality Gates (Citation Manager)
- No fabricated sources
- High-impact claims flagged
- Incomplete evidence clearly marked
- Reproducible referencing (keys consistent)

### [TEST] Acceptance Checks
- Every non-trivial factual claim has either a citation or `[ASSUMPTION]`
- Bibliography keys are consistent and unique
