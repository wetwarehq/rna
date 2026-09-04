"""Load packaged auth.seq for a claimed cassette id."""

from __future__ import annotations

import re
from importlib.resources import files

from ..seq import normalize

SLOT_NAMES = ("5_utr", "cds", "3_utr", "polya")
_SAFE = re.compile(r"^[A-Za-z0-9._-]+$")


def parse_fa(text: str) -> dict[str, str]:
    slots: dict[str, str] = {}
    name: str | None = None
    buf: list[str] = []
    for line in text.splitlines():
        if line.startswith(">"):
            if name is not None:
                rna, _, _ = normalize("".join(buf))
                slots[name] = rna
            name = line[1:].strip().split()[0]
            buf = []
        else:
            buf.append(line.strip())
    if name is not None:
        rna, _, _ = normalize("".join(buf))
        slots[name] = rna
    return slots


def load_auth(claimed_id: str | None) -> dict[str, str] | None:
    """Return slot → RNA for a packaged cassette, or None if unset/unknown."""
    if not claimed_id:
        return None
    if not _SAFE.match(claimed_id):
        return None
    root = files("rna.auth")
    target = f"{claimed_id}.fa"
    if not root.joinpath(target).is_file():
        return None
    return parse_fa(root.joinpath(target).read_text(encoding="utf-8"))
