# Contributing to PEN-ASSEMBLE

Thank you for your interest in contributing. PEN-ASSEMBLE is a research pipeline for IS110-family bridge recombinase design. Contributions that improve reproducibility, extend the scoring framework, or add new design strategies are especially welcome.

---

## Development Setup

```bash
git clone https://github.com/ahmedanees-m/pen-assemble.git
cd pen-assemble
pip install -e ".[dev,docs]"
```

This installs the package in editable mode with all development and documentation dependencies.

---

## Running Tests

```bash
pytest tests/ -v
```

All 63 tests must pass before submitting a pull request. The CI matrix runs Python 3.10, 3.11, and 3.12.

---

## Code Style

This project uses [Ruff](https://github.com/astral-sh/ruff) for linting and formatting, and [mypy](https://mypy-lang.org/) for type checking.

```bash
# Check and fix style
ruff check pen_assemble/ tests/ --fix
ruff format pen_assemble/ tests/

# Type check
mypy pen_assemble/ --ignore-missing-imports
```

The CI lint job runs both checks on every push. Pull requests with lint failures will not be merged.

---

## Submitting a Pull Request

1. Fork the repository and create a feature branch from `main`.
2. Make your changes with appropriate tests.
3. Ensure `pytest`, `ruff check`, `ruff format --check`, and `mypy` all pass locally.
4. Open a pull request with a clear description of what changed and why.

For significant changes (new scoring axes, new strategies, algorithm modifications), please open an issue first to discuss the approach.

---

## Reporting Bugs

Please open a [GitHub Issue](https://github.com/ahmedanees-m/pen-assemble/issues) with:
- A minimal reproducible example
- The Python version and OS
- The full error traceback

---

## Scientific Contributions

If you are extending the PenScore framework or adding a new design strategy:
- Document the scientific rationale clearly in the PR description
- Add a corresponding entry to `DESIGN_PROVENANCE.md` if any change to the documented methodology is introduced
- Include unit tests that cover the new scoring logic

---

## License

By contributing, you agree that your contributions will be licensed under the [MIT License](LICENSE).
