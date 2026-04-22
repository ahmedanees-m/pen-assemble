""CLI for pen-assemble.""
from __future__ import annotations

import click

@click.group()
@click.version_option()
def main() -> int:
    ""$repo CLI.""
    return 0

if __name__ == '__main__':
    main()

