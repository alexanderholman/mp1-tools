#!/bin/bash
# validate_agents.sh
# Validation script for Agent Factory
# This script runs all MUST requirement tests

# Note: Do not use 'set -e' here; tests are expected to fail without aborting the script.

echo "======================================"
echo "Agent Factory Validation Suite"
echo "======================================"
echo ""

PASSED=0
FAILED=0
TOTAL=0

# Color codes
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Check if agents/ directory exists when agents are defined in agents.yaml
if [ -f "agents.yaml" ]; then
    agent_count=$(grep -c '  - id:' agents.yaml 2>/dev/null || echo "0")
    if [ "$agent_count" -gt 0 ] && [ ! -d "agents" ]; then
        echo -e "${YELLOW}Warning: agents.yaml defines $agent_count agent(s) but agents/ directory does not exist.${NC}"
        echo "Create the directory with: mkdir -p agents"
        echo ""
    fi
fi

# Test function
run_test() {
    local test_id=$1
    local test_name=$2
    local test_command=$3
    
    TOTAL=$((TOTAL + 1))
    echo "Running ${test_id}: ${test_name}"
    
    if eval "$test_command"; then
        echo -e "${GREEN}✓ PASS${NC}: ${test_id}"
        PASSED=$((PASSED + 1))
    else
        echo -e "${RED}✗ FAIL${NC}: ${test_id}"
        FAILED=$((FAILED + 1))
    fi
    echo ""
}

# Helper: verify all agent files referenced in agents.yaml exist
check_agent_files_exist() {
    if [ ! -f "agents.yaml" ]; then
        return 0  # No agents.yaml, nothing to check
    fi
    
    # Extract file_path entries from agents.yaml and check they exist
    local missing_files=0
    while IFS= read -r filepath; do
        if [ -n "$filepath" ] && [ ! -f "$filepath" ]; then
            echo "Missing file: $filepath"
            missing_files=$((missing_files + 1))
        fi
    done < <(grep 'file_path:' agents.yaml | awk '{print $2}' | tr -d '"')
    
    [ $missing_files -eq 0 ]
}

# TEST-007-1: Verify agent files exist (replaces flat structure test)
echo "=== SPEC-007: Flexible Directory Structure ==="
if [ -d "agents" ]; then
    run_test "TEST-007-1" \
        "Verify all agent files referenced in agents.yaml exist" \
        "check_agent_files_exist"
else
    echo "Note: agents/ directory does not exist yet - skipping TEST-007-1"
    echo ""
fi

# TEST-002-1: Verify required headings in agent files (presence and order)
echo "=== SPEC-002: Agent File Format ==="
if [ -d "agents" ] && [ "$(find agents/ -name '*.md' -type f 2>/dev/null | wc -l)" -gt 0 ]; then
    # Find all .md files in agents/ directory (including nested)
    while IFS= read -r file; do
        [ -f "$file" ] || continue
        # Check both presence and order of headings
        run_test "TEST-002-1" \
            "Verify required headings in $(basename $file)" \
            "grep -n '^## Purpose' '$file' > /dev/null && \
             grep -n '^## Inputs' '$file' > /dev/null && \
             grep -n '^## Outputs' '$file' > /dev/null && \
             grep -n '^## Behavior' '$file' > /dev/null && \
             grep -n '^## Constraints' '$file' > /dev/null && \
             [ \$(grep -n '^## Purpose' '$file' | cut -d: -f1) -lt \$(grep -n '^## Inputs' '$file' | cut -d: -f1) ] && \
             [ \$(grep -n '^## Inputs' '$file' | cut -d: -f1) -lt \$(grep -n '^## Outputs' '$file' | cut -d: -f1) ] && \
             [ \$(grep -n '^## Outputs' '$file' | cut -d: -f1) -lt \$(grep -n '^## Behavior' '$file' | cut -d: -f1) ] && \
             [ \$(grep -n '^## Behavior' '$file' | cut -d: -f1) -lt \$(grep -n '^## Constraints' '$file' | cut -d: -f1) ]"
    done < <(find agents/ -name '*.md' -type f)
    
    # TEST-002-2: Verify heading format follows markdown H2 convention
    while IFS= read -r file; do
        [ -f "$file" ] || continue
        run_test "TEST-002-2" \
            "Verify H2 heading format in $(basename $file)" \
            "grep -E '^## (Purpose|Inputs|Outputs|Behavior|Constraints)$' '$file' > /dev/null"
    done < <(find agents/ -name '*.md' -type f)
else
    echo "Note: No agent markdown files found yet - skipping TEST-002-1 and TEST-002-2"
    echo ""
fi

# TEST-003-1: Verify each agent has at least one tag
# TEST-003-2: Verify agent IDs are unique
# TEST-003-3: Verify tags are from allowed list
echo "=== SPEC-003: Tags and Metadata ==="
if [ -f "agents.yaml" ]; then
    # TEST-003-1: Check each agent has at least one tag
    run_test "TEST-003-1" \
        "Verify each agent has at least one tag" \
        "[ -z \"\$(awk '/  - id:/{flag=1; next} flag && /tags:/{flag=2; next} flag==2 && /^[[:space:]]*-/{flag=0; next} flag==2 && /^  - id:/{print \"no_tags\"; flag=0}' agents.yaml | grep no_tags)\" ]"
    
    # TEST-003-2: Check agent IDs are unique
    run_test "TEST-003-2" \
        "Verify agent IDs are unique in agents.yaml" \
        "! grep 'id:' agents.yaml | awk '{print \$3}' | sed 's/\"//g' | sort | uniq -d | grep -q ."
    
    # TEST-003-3: Verify tags are from allowed list
    if grep -q 'allowed_tags:' agents.yaml; then
        # Extract allowed tags and agent tags, then compare
        allowed_tags=$(awk '/^allowed_tags:/,/^# / {if (/^  - /) print $2}' agents.yaml | tr '\n' '|' | sed 's/|$//')
        if [ -n "$allowed_tags" ]; then
            run_test "TEST-003-3" \
                "Verify all tags are from allowed list" \
                "! awk '/^agents:/,/^# Validation/ {if (/^    tags:/) {flag=1; next} if (flag && /^      - /) {print \$2} if (flag && /^    [a-z_]/) flag=0}' agents.yaml | grep -vE \"^($allowed_tags)\$\" | grep -q ."
        else
            echo "Note: No allowed_tags found - skipping TEST-003-3"
        fi
    else
        echo "Note: No allowed_tags defined in agents.yaml - skipping TEST-003-3"
    fi
else
    echo "Error: agents.yaml not found"
    FAILED=$((FAILED + 3))
    TOTAL=$((TOTAL + 3))
    echo ""
fi

# TEST-004-1: Verify append-only files are only appended to
echo "=== SPEC-004: Append-Only Files ==="
echo "Note: TEST-004-1 requires git history - checking files exist"
run_test "TEST-004-1" \
    "Verify append-only files exist" \
    "[ -f specs.md ] && [ -f agent_runs.md ] && [ -f decisions.md ]"

# TEST-005-1: Citations verification (manual review required)
echo "=== SPEC-005: No Fabrication ==="
echo "Note: TEST-005-1 (Citation verification) requires manual review"
echo ""

# TEST-006-1: Verify documentation uses .md extension
echo "=== SPEC-006: Markdown Output ==="
run_test "TEST-006-1" \
    "Verify key documentation files use .md extension" \
    "[ -f agents.md ] && [ -f specs.md ] && [ -f agent_runs.md ] && [ -f decisions.md ]"

# Summary
echo "======================================"
echo "Validation Summary"
echo "======================================"

if [ $TOTAL -eq 0 ]; then
    echo -e "${YELLOW}Warning: No tests were executed!${NC}"
    exit 1
fi

echo "Total tests run: $TOTAL"
echo -e "${GREEN}Passed: $PASSED${NC}"
if [ $FAILED -gt 0 ]; then
    echo -e "${RED}Failed: $FAILED${NC}"
    echo ""
    echo "Please fix the failing tests before proceeding."
    exit 1
else
    echo -e "${GREEN}All tests passed!${NC}"
    exit 0
fi
