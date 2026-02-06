# GitHub Copilot Workspace Instructions

## Repository Overview
This template ships with AgentFactory-compatible agent definitions, validation tooling, and append-only logs. It provides a flexible, organized structure compatible with major AI platforms, and is designed for plug-and-play use with OpenCode and GitHub Copilot.

## Key Architecture Principles

### 1. Flexible Directory Structure (SHOULD)
- Agent files SHOULD be organized in the `agents/` directory
- Nested subdirectories are ALLOWED for better organization (e.g., `agents/testing/`, `agents/security/`)
- Configuration and documentation files are typically stored at the repository root
- Organization SHOULD follow conventions compatible with GitHub Copilot, OpenAI ChatGPT, Google Gemini, agent-based IDEs (OpenCode.ai), and Google Colab
- Directory structure SHOULD be logical and self-documenting

### 2. Required Agent File Format (MUST)
Every agent file MUST include these headings in this exact order:
1. `## Purpose` - What the agent does and why it exists
2. `## Inputs` - Required and optional inputs the agent needs
3. `## Outputs` - What the agent produces
4. `## Behavior` - How the agent processes inputs to produce outputs
5. `## Constraints` - Limitations, boundaries, and operational constraints

### 3. Append-Only Files (MUST)
The following files are append-only and MUST NOT have content removed or modified:
- `specs.md` - Technical specifications and requirements
- `agent_runs.md` - Log of agent execution runs
- `decisions.md` - Architectural and design decisions

**Only add new entries at the end of these files. Never modify or delete existing entries.**

### 4. No Fabrication Policy (MUST)
- MUST NOT fabricate citations, references, or test results
- All claims MUST be verifiable and traceable to sources
- Use `[ASSUMPTION]` tags for speculative content
- Provide verifiable sources for all external information

### 5. Tagging Conventions
Use these tags consistently throughout documentation:
- `[SPEC]` - Specification entries
- `[ASSUMPTION]` - Assumptions made
- `[RISK]` - Identified risks
- `[TODO]` - Pending tasks
- `[DECISION]` - Design decisions
- `[TEST]` - Test cases
- `[DONE]` - Completed items

## Core Files and Their Purpose

### Configuration Files
- **agents.yaml** - Central agent registry and configuration schema
  - Defines all agents with metadata (id, name, description, tags, status)
  - Specifies validation rules (MUST/SHOULD/MAY requirements)
  - Lists allowed tags for categorization

- **agents.md** - Agent documentation rules and guidelines
  - Template for creating new agents
  - MUST/SHOULD/MAY requirements explained
  - Validation procedures

### Agent Files
Located in `agents/` directory:
- **Architect.md** - Spec Author + System Designer
- **Builder.md** - Implementer / Artifact Producer
- **Skeptic.md** - Adversarial Reviewer / Breaker
- **Editor.md** - Clarity + Structure Editor
- **ProjectManager.md** - Packaging + Orchestration
- **CitationOfficer.md** - Evidence Auditor + Claim Tracker
- **ChatGPT.md** - Generalist Execution Agent
- **OpenCodeManager.md** - Session orchestration for OpenCode

### Documentation Files
- **README.md** - Quick start and overview
- **specs.md** - Technical specifications (append-only)
- **decisions.md** - Design decisions (append-only)
- **agent_runs.md** - Execution log (append-only)
- **INSTALL.md** - Installation instructions
- **CHANGELOG.md** - Version history

## Development Workflow

### Adding a New Agent
1. Define the agent in `agents.yaml` with all required fields:
   - Unique ID
   - Name and description
   - At least one tag from the allowed list
   - Status (active/inactive)
   - File path in agents/ directory
   
2. Create the agent markdown file in `agents/` directory following the template
   
3. Ensure all five required headings are present in correct order
   
4. Run `./validate_agents.sh` to verify compliance
   
5. Document the addition in `decisions.md` (append-only)

### Modifying Existing Agents
1. Update the agent definition in `agents.yaml` if metadata changes
2. Update the markdown file with changes
3. Increment the version number
4. Document changes in agent's version history
5. Record the decision in `decisions.md` (append-only)
6. Run validation: `./validate_agents.sh`

### Running Validation
```bash
./validate_agents.sh
```

This script validates:
- Required headings in agent files
- Unique agent IDs
- Proper tag usage
- File existence
- Markdown format compliance

## Code Style and Standards

### Markdown
- Use GitHub Flavored Markdown
- Use code blocks with language specification
- Use H2 (`##`) for required agent headings
- Keep lines reasonably short for readability

### YAML
- Use 2-space indentation
- Follow schema defined in `agents.yaml`
- Validate YAML syntax before committing

### Versioning
- Follow Semantic Versioning (SemVer 2.0.0)
- Document version changes in agent files
- Update CHANGELOG.md for releases

## Important Constraints

### What You MUST Do
- Validate all changes with `./validate_agents.sh` before committing
- Include all five required headings in agent files
- Use unique IDs for all agents
- Tag all agents with at least one allowed tag
- Append to append-only files (never modify existing content)
- Mark assumptions clearly with `[ASSUMPTION]` tags

### What You MUST NOT Do
- Modify or delete content in append-only files (specs.md, decisions.md, agent_runs.md)
- Fabricate citations or test results
- Remove or override working code without explicit reason
- Create agents without required headings
- Use duplicate agent IDs

## Copilot Agent Integration

This repository is designed to work seamlessly with GitHub Copilot Agents. The agent definitions in the `agents/` directory can be:

1. **Referenced by Copilot** to understand system architecture and roles
2. **Used as templates** for creating new AI agent definitions
3. **Validated automatically** via the validation script in CI/CD

When working with GitHub Copilot Workspace:
- Copilot can read agent definitions to understand the agent factory pattern
- Copilot can help create new agent definitions following the established format
- Copilot can validate changes against the specifications
- Copilot respects the flat-file structure and append-only constraints

## Testing and Quality

### Validation Tests
Run the validation suite:
```bash
./validate_agents.sh
```

Expected output: All tests should pass with green checkmarks.

### Manual Review Requirements
- Citation verification (TEST-005-1) requires human review
- Ensure all external references are accurate
- Verify test results match actual execution

## Questions and Support

For questions or issues:
1. Check the [SUPPORT.md](/.github/SUPPORT.md) guide
2. Review [agents.md](/agents.md) for agent documentation rules
3. Check [specs.md](/specs.md) for technical specifications
4. Refer to [decisions.md](/decisions.md) for design rationale

## License
See [LICENSE](/LICENSE) for usage rights and restrictions.
