"""Sequence normalisation. Never used as an annotator."""

from __future__ import annotations


def normalize(text: str) -> tuple[str, list[str], str]:
    """Strip FASTA, whitespace, and digits; convert T→U. Returns (rna, notes, raw_compact)."""
    notes: list[str] = []
    s = (text or "").strip()
    if s.startswith(">"):
        lines = s.splitlines()
        body = [ln.strip() for ln in lines[1:] if not ln.startswith(">")]
        s = "".join(body)
        notes.append("Parsed a FASTA header and used the first record.")
    raw_chars = []
    for ch in s:
        if ch.isspace() or ch.isdigit():
            continue
        raw_chars.append(ch)
    raw = "".join(raw_chars)
    upper = raw.upper()
    t_count = upper.count("T")
    u_count = upper.count("U")
    if t_count and u_count:
        notes.append(f"Mixed T ({t_count}) and U ({u_count}); converting T to U.")
    elif t_count:
        notes.append(f"Converted {t_count} T residue(s) to U (DNA template → RNA).")
    rna = upper.replace("T", "U")
    return rna, notes, raw.upper()
