# Workflow: PR Review and Merge Gate

## Workflow ID
pr-review-and-merge-001

## Purpose
Ensure PRs receive automated review where available, are validated, and achieve a clear go/no-go decision before merge.

## Sequence
1. Collect PR context (diff, checks, requested reviewers).
2. Request automated review (Copilot/Codex) using the correct handle or trigger comment.
3. Run required validation/tests.
4. Summarize findings and provide go/no-go.
5. If go, wait for user confirmation before merge.

## Inputs/Outputs
- Inputs: PR URL/number, repo, validation commands.
- Outputs: PR summary, review status, go/no-go recommendation.

## Validation
- Validation scripts pass or failures are explained.
- Automated review request is visible in PR timeline.
