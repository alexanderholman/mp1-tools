# AGENTS.md

This file provides practical guidance for AI agents working on repositories created from this template.

## Project Overview

This template is preloaded with AgentFactory-style agent definitions, validation, and append-only logs so OpenCode and Copilot can work plug-and-play.

## Key Files and Structure

- **agents/** - Agent definition files (markdown format)
  - Can be flat or nested in subdirectories
  - Each agent has 5 required headings: Purpose, Inputs, Outputs, Behavior, Constraints
  
- **specialisms/** - Domain-specific standards and guidelines
  - Testing.md, Security.md, Coder.md, etc.
  
- **agents.yaml** - Central registry
  - All agents MUST be registered here with unique IDs
  - Defines tags, file paths, and metadata
  
- **Append-only files** (DO NOT modify existing content, only append):
  - `specs.md` - Technical specifications
  - `decisions.md` - Architectural decisions
  - `agent_runs.md` - Execution logs

- **traits/** - Reusable behavioral modules for agents/specialisms/tasks
- **workflows/** - Task-level execution patterns and validation gates
- **tasks/** - Recurring and stateful task tracking (cadence, wip, triage)

## Development Workflow

### Adding a New Agent

1. Create the agent markdown file in `agents/` (or subdirectory)
2. Include all 5 required headings in order:
   - `## Purpose`
   - `## Inputs`
   - `## Outputs`
   - `## Behavior`
   - `## Constraints`
3. Add entry to `agents.yaml` with:
   - Unique ID (e.g., "my-agent-001")
   - Name, description, version
   - At least one tag from allowed_tags list
   - File path relative to repository root
   - Status (active/inactive)
4. Run validation: `./validate_agents.sh`
5. If making architectural decisions, append to `decisions.md`

### Validation

Always run before committing:
```bash
./validate_agents.sh
```

This validates:
- Required headings present and in correct order
- Unique agent IDs
- Tags from allowed list
- File paths exist
- Markdown format

### Python Validation (Alternative)

For more detailed validation:
```python
import yaml

with open('agents.yaml', 'r') as f:
    data = yaml.safe_load(f)
    
# Check unique IDs
ids = [agent['id'] for agent in data['agents']]
assert len(ids) == len(set(ids)), "Duplicate agent IDs found"
```

## Coding Conventions

### Agent Files (Markdown)
- Use GitHub Flavored Markdown
- Use H2 (`##`) for required headings
- Use H3 (`###`) for subsections
- Include code blocks with language specification
- Keep lines reasonably short for readability

### YAML Files
- Use 2-space indentation
- Quote strings containing special characters
- Validate YAML syntax before committing

### Append-Only Files
- **NEVER** delete or modify existing content
- **ALWAYS** add new entries at the end
- Include date stamps in entry headers
- Use horizontal rules (`---`) to separate entries

## Testing

Run all validation tests:
```bash
./validate_agents.sh
```

Check specific aspects:
```bash
# Check for required headings in a specific file
grep -E "^## (Purpose|Inputs|Outputs|Behavior|Constraints)" agents/MyAgent.md

# Validate YAML syntax
python3 -c "import yaml; yaml.safe_load(open('agents.yaml'))"

# Check for nested files (if needed)
find agents/ -name "*.md" -type f
```

## Common Tasks

### View Current Agents
```bash
# List all agent files
find agents/ -name "*.md" -type f

# View agents in registry
grep "  - id:" agents.yaml
```

### Check Append-Only File Status
```bash
# View recent additions to specs.md
tail -50 specs.md

# View recent decisions
tail -50 decisions.md

# View recent agent runs
tail -50 agent_runs.md
```

### Update Agent Metadata
1. Edit `agents.yaml` to update metadata
2. Increment version number
3. Update `last_modified` date
4. Run `./validate_agents.sh`
5. Append decision to `decisions.md` explaining why

## Important Rules

### MUST DO
- Include all 5 required headings in agent files
- Register agents in agents.yaml with unique IDs
- Use tags from the allowed_tags list
- Run validation before committing
- Append to append-only files (never modify existing content)
- Mark assumptions with `[ASSUMPTION]` tags

### MUST NOT DO
- Modify or delete content in append-only files (specs.md, decisions.md, agent_runs.md)
- Fabricate citations or test results
- Use duplicate agent IDs
- Create agents without required headings
- Skip validation before committing

## Directory Structure

Both flat and nested structures are supported:

```
# Flat structure (original, still valid)
agents/
├── Architect.md
├── Builder.md
├── Tester.md
└── SecurityReviewer.md

# Nested structure (now supported)
agents/
├── core/
│   ├── Architect.md
│   └── Builder.md
├── quality/
│   ├── Tester.md
│   └── SecurityReviewer.md
└── documentation/
    └── Editor.md
```

The agents.yaml file paths should match the actual file locations.

## Troubleshooting

### Validation Fails
- Check that all 5 headings are present and spelled correctly
- Verify agent ID is unique in agents.yaml
- Ensure tags are from the allowed_tags list
- Confirm file path in agents.yaml matches actual file location

## Further Documentation

- `.github/copilot-instructions.md` - Detailed instructions for GitHub Copilot
- `agents.md` - Agent documentation guidelines and template
- `specs.md` - Technical specifications (append-only)
- `decisions.md` - Architectural decisions (append-only)
- `README.md` - Project overview and quick start

## Platform Compatibility

This repository structure is compatible with:
- GitHub Copilot Workspace
- OpenAI ChatGPT (this AGENTS.md file)
- Google Gemini
- Google Colab
- Agent-based IDEs (OpenCode.ai)
