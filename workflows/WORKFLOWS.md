# Workflows

Workflows are task-level execution patterns that define sequencing, agent handoffs, and required validation steps. Attach a workflow to a task when it must follow a specific pipeline.

## Workflow Format
Each workflow MUST include:

- **Workflow ID** (unique)
- **Purpose** (what the workflow achieves)
- **Sequence** (ordered steps and agent roles)
- **Inputs/Outputs** (artifacts and checkpoints)
- **Validation** (required checks)

Use this template:

```markdown
# Workflow: <Name>

## Workflow ID
<workflow-id>

## Purpose
<what the workflow achieves>

## Sequence
1. <step>
2. <step>

## Inputs/Outputs
- Inputs: <...>
- Outputs: <...>

## Validation
- <check 1>
- <check 2>
```
