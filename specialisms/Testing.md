<!-- filename: Testing.md -->

## [SPEC] Specialism Addendum — Testing (v1.0)

### Purpose
Defines standards for testing work: test creation, test execution, test documentation, and quality validation.

### Operating Rules
- Write tests before or alongside code (TDD/BDD when appropriate).
- Test one thing per test (single responsibility).
- Use clear, descriptive test names that explain what is being tested.
- Follow Arrange-Act-Assert (AAA) pattern or equivalent structure.
- Make tests independent (no execution order dependencies).
- Use meaningful assertions with clear failure messages.
- Test the happy path, edge cases, and error conditions.
- Keep tests fast and deterministic (avoid flaky tests).
- Mock external dependencies to isolate units under test.
- Document test requirements and test data.

### Outputs (typical)
- test plans
- unit tests
- integration tests
- acceptance tests
- test reports
- coverage reports
- test data and fixtures

### Quality Gates (Testing)
- All MUST requirements have at least one test
- Tests have clear, descriptive names
- Tests follow AAA pattern or equivalent
- Tests are independent and repeatable
- Edge cases and error conditions tested
- Test coverage metrics provided
- Tests pass consistently (no flaky tests)
- Test execution steps documented

### [TEST] Acceptance Checks
- Every MUST requirement can be validated by running a test
- Tests can be executed following documented steps
- Test results clearly show pass/fail status
- Coverage report shows which requirements are tested
- Test failures provide enough information to diagnose issues

### Testing Checklist
- [ ] Test plan defines scope and approach
- [ ] MUST requirements are covered by tests
- [ ] Happy path scenarios tested
- [ ] Edge cases and boundary conditions tested
- [ ] Error handling and failure scenarios tested
- [ ] Tests have clear, descriptive names
- [ ] Tests are independent and isolated
- [ ] Test data and fixtures documented
- [ ] Tests execute successfully
- [ ] Coverage metrics calculated
- [ ] Test report documents results

### Test Types
**Unit Tests**
- Test individual functions or methods
- Fast execution (milliseconds)
- No external dependencies (mocked)
- High code coverage expected (>80%)

**Integration Tests**
- Test component interactions
- Moderate execution time (seconds)
- May use test doubles or real dependencies
- Focus on interfaces and contracts

**Acceptance Tests**
- Test requirement satisfaction
- Slower execution (seconds to minutes)
- End-to-end scenarios
- Validate business requirements

**Security Tests**
- Test security requirements
- Authentication and authorization tests
- Input validation and injection tests
- Secrets and sensitive data protection tests

### Test Documentation Standards
Each test should include:
- **Test ID**: Unique identifier
- **Test Name**: Clear, descriptive name
- **Test Description**: What is being tested and why
- **Preconditions**: Setup required before test
- **Test Steps**: How to execute the test
- **Expected Results**: What success looks like
- **Actual Results**: What actually happened (during execution)
- **Pass/Fail Status**: Clear outcome

### Test Naming Conventions
Use descriptive names that explain the test:
- `test_<function>_<scenario>_<expected_result>`
- Example: `test_login_with_valid_credentials_succeeds`
- Example: `test_login_with_invalid_password_returns_401`
- Example: `test_calculate_total_with_empty_cart_returns_zero`

### Coverage Metrics
- **Requirement Coverage**: % of requirements tested
- **Code Coverage**: % of code executed by tests (when applicable)
- **Branch Coverage**: % of decision branches tested
- **Pass Rate**: % of tests passing

### Test Maintenance
- Update tests when requirements change
- Remove obsolete tests
- Refactor tests for clarity
- Keep test code clean and maintainable
- Version test data with tests
