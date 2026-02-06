# Example Documentation Agent

## Purpose
This is an example agent that demonstrates the required structure for agent documentation in the Agent Factory. It serves as a template and reference for creating new agents.

## Inputs
- **Source Documents**: List of markdown files to process
- **Format Specification**: Desired output format (markdown, HTML, etc.)
- **Validation Rules**: Set of rules to validate documentation against

## Outputs
- **Processed Documentation**: Formatted and validated documentation files
- **Validation Report**: Summary of validation results
- **Error Log**: List of any errors or warnings encountered

## Behavior
The Example Documentation Agent processes documentation files through the following steps:

1. **Intake**: Receives source documents and configuration
2. **Parsing**: Parses markdown files to extract structure and content
3. **Validation**: Checks documents against specified rules
4. **Processing**: Applies formatting and transformations
5. **Output**: Generates final documentation and reports

The agent validates the following:
- Required headings are present
- Links are valid and not broken
- Code blocks have language specifications
- Images have alt text

## Constraints
- **File Size**: Maximum file size of 10MB per document
- **Format**: Only supports markdown (.md) files as input
- **Performance**: Processes up to 100 files per run
- **Dependencies**: Requires markdown parser library
- **Execution Time**: Maximum 5 minutes per run

## Tags
Tags defined in agents.yaml: documentation, automation, utility

## Version History
- v1.0.0 (2026-01-28): Initial version for demonstration purposes
