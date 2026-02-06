<!-- filename: Security.md -->

## [SPEC] Specialism Addendum — Security (v1.0)

### Purpose
Defines standards for security work: security analysis, vulnerability assessment, secure coding practices, and compliance checking.

### Operating Rules
- Follow principle of least privilege (grant minimum necessary permissions).
- Apply defense in depth (multiple layers of security).
- Use secure defaults (opt-in for risky features, not opt-out).
- Fail securely (errors should not expose sensitive information).
- Never trust user input (validate, sanitize, encode).
- Keep security simple (complexity is the enemy of security).
- Use well-tested security libraries (don't roll your own crypto).
- Separate duties (authorization checks separate from business logic).
- Log security events (authentication, authorization, data access).
- Encrypt sensitive data at rest and in transit.
- Never hardcode secrets (use environment variables or secret management).
- Keep dependencies updated (patch known vulnerabilities).

### Outputs (typical)
- security assessment reports
- vulnerability lists
- threat models
- security test cases
- remediation recommendations
- compliance checklists
- secure coding guidelines

### Quality Gates (Security)
- No critical or high severity vulnerabilities in production
- No hardcoded secrets or credentials
- Input validation on all user-supplied data
- Authentication and authorization properly implemented
- Sensitive data encrypted (in transit and at rest)
- Dependencies scanned for known vulnerabilities
- Security logging and monitoring configured
- Security tests pass
- OWASP Top 10 compliance

### [TEST] Acceptance Checks
- Security assessment completed and documented
- All identified vulnerabilities tracked and prioritized
- Critical and high severity vulnerabilities remediated
- Security test cases cover authentication, authorization, and input validation
- No secrets in source code or configuration files
- Dependencies are up-to-date or vulnerabilities mitigated

### OWASP Top 10 (2021) Checklist
- [ ] A01:2021 – Broken Access Control
  - Authorization checks on all protected resources
  - Deny by default
  - Rate limiting on APIs
- [ ] A02:2021 – Cryptographic Failures
  - Sensitive data encrypted in transit (TLS)
  - Sensitive data encrypted at rest
  - Strong encryption algorithms (AES-256, etc.)
  - Proper key management
- [ ] A03:2021 – Injection
  - Input validation on all user data
  - Parameterized queries (SQL injection prevention)
  - Output encoding (XSS prevention)
  - Command injection prevention
- [ ] A04:2021 – Insecure Design
  - Threat modeling completed
  - Security requirements defined
  - Secure design patterns used
- [ ] A05:2021 – Security Misconfiguration
  - Secure defaults
  - Minimal features enabled
  - Error messages don't leak information
  - Security headers configured
- [ ] A06:2021 – Vulnerable and Outdated Components
  - Dependencies up-to-date
  - Known vulnerabilities patched
  - Dependency scanning in place
- [ ] A07:2021 – Identification and Authentication Failures
  - Strong password requirements
  - Multi-factor authentication (when appropriate)
  - Session management secure
  - No default credentials
- [ ] A08:2021 – Software and Data Integrity Failures
  - Code integrity checks
  - Secure update mechanisms
  - Supply chain security
- [ ] A09:2021 – Security Logging and Monitoring Failures
  - Security events logged
  - Log integrity protected
  - Alerting configured
  - Audit trails maintained
- [ ] A10:2021 – Server-Side Request Forgery (SSRF)
  - URL validation
  - Network segmentation
  - Allow-lists for destinations

### Security Severity Classification
**Critical**
- Remote code execution
- Authentication bypass
- SQL injection with data access
- Exposure of secrets or credentials

**High**
- Privilege escalation
- Sensitive data exposure
- Cross-site scripting (XSS) on sensitive pages
- Insecure direct object references

**Medium**
- Information disclosure
- Missing security headers
- Weak cryptography
- Insecure session management

**Low**
- Security misconfiguration (non-exploitable)
- Missing security best practices
- Outdated dependencies (no known exploit)

### Secure Coding Practices
**Input Validation**
- Validate all input (type, length, format, range)
- Use allow-lists instead of deny-lists
- Reject invalid input, don't try to clean it
- Validate on the server side (never trust client validation)

**Authentication**
- Use strong password requirements
- Implement rate limiting on login attempts
- Use secure session management
- Implement logout functionality
- Expire sessions after inactivity

**Authorization**
- Check authorization on every request
- Use role-based or attribute-based access control
- Implement principle of least privilege
- Don't rely on client-side checks

**Cryptography**
- Use TLS 1.2 or higher for data in transit
- Use AES-256 or equivalent for data at rest
- Use secure random number generators
- Don't implement custom cryptography
- Store passwords with strong hashing (bcrypt, Argon2)

**Data Protection**
- Classify data by sensitivity
- Encrypt sensitive data
- Minimize data collection
- Implement data retention policies
- Secure data deletion

**Error Handling**
- Don't expose sensitive information in errors
- Log errors securely
- Use generic error messages for users
- Handle exceptions properly

**Logging**
- Log authentication events
- Log authorization failures
- Log sensitive data access
- Don't log sensitive data (passwords, tokens, PII)
- Protect log integrity

### Secrets Management
- Never commit secrets to version control
- Use environment variables or secret vaults
- Rotate secrets regularly
- Implement secret scanning in CI/CD
- Use different secrets for different environments

### Dependency Security
- Keep dependencies updated
- Use dependency scanning tools
- Review security advisories
- Pin dependency versions
- Audit third-party code
