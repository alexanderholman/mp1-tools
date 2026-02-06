# Trait: PR Authorship Clarity

## Trait ID
pr-authorship-clarity-001

## Purpose
Ensure PR comments, reviews, and actions clearly reflect their true author (human vs automation) and avoid misattributed authorship.

## Applies To
- task
- agent

## Behavior
- MUST state when actions are performed using a user-authenticated token.
- MUST avoid implying automated/bot authorship unless a bot identity is actually used.
- SHOULD prefer bot/app credentials for automated reviews or comments.
- MUST document the chosen identity in the PR summary (e.g., “commented as alexanderholman via user token”).

## Checks
- PR comments indicate author identity and method.
- If a bot identity is desired, `gh auth status` shows bot/app credentials.

## Overrides
- Allowed when user explicitly requests comments as their account.
