"""CLI for pen-assemble."""

from __future__ import annotations

import click


@click.group()
@click.version_option()
def main() -> None:
    """PEN-ASSEMBLE - computational design of programmable genome-writing editors."""


if __name__ == "__main__":
    main()
