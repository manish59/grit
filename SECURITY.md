# Security Policy

## Supported Versions

| Version | Supported          |
| ------- | ------------------ |
| 0.1.x   | :white_check_mark: |

## Reporting a Vulnerability

If you discover a security vulnerability in GRIT, please report it responsibly:

1. **Do not** open a public GitHub issue for security vulnerabilities
2. Email the maintainer directly with details of the vulnerability
3. Include steps to reproduce the issue if possible
4. Allow reasonable time for a fix before public disclosure

## Security Considerations

GRIT is a command-line tool for processing genomic interval data. Security considerations include:

### Input Validation
- GRIT validates BED file format but does not sanitize file paths
- Users should ensure input files come from trusted sources

### File System Access
- GRIT reads and writes files as specified by command-line arguments
- No network access or remote file capabilities

### Memory Safety
- Written in Rust, which provides memory safety guarantees
- No unsafe code blocks in core algorithms

## Response Timeline

- **Acknowledgment**: Within 48 hours
- **Initial Assessment**: Within 1 week
- **Fix Timeline**: Depends on severity, typically 1-4 weeks
