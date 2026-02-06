# Trait: PR Review Automation

## Trait ID
pr-review-automation-001

## Purpose
Standardize how automated reviews are requested and handled (Copilot, Codex, or other bots), and ensure review provenance is explicit.

## Applies To
- task
- agent

## Behavior
- MUST request automated review using the correct reviewer handle or trigger comment.
- MUST record which reviewer handle was used and whether the request succeeded.
- SHOULD include a fallback path if reviewer handle is invalid (e.g., comment trigger).
- MUST not claim a review was requested if the command failed.

## Checks
- PR timeline shows the automated reviewer or trigger comment.
- CLI output confirms reviewer request success.

## Overrides
- Allowed if repository does not support automated reviewers; document why.
