# Traits

Traits are reusable behavioral modules that can be attached to agents, specialisms, or tasks. They define constraints, expectations, and operational guardrails that shape how work is executed.

## Trait Format
Each trait MUST include:

- **Trait ID** (unique)
- **Purpose** (what it enforces)
- **Applies To** (agents, specialisms, tasks)
- **Behavior** (required behaviors)
- **Checks** (how to confirm adherence)
- **Overrides** (how to disable or override)

Use this template:

```markdown
# Trait: <Name>

## Trait ID
<trait-id>

## Purpose
<what this trait enforces>

## Applies To
- agent
- specialism
- task

## Behavior
- MUST ...
- SHOULD ...
- MUST NOT ...

## Checks
- <check 1>
- <check 2>

## Overrides
- <how to override, if allowed>
```
