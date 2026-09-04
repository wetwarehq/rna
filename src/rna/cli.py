"""Command line: rna verify [card.json|-]"""

from __future__ import annotations

import argparse
import json
import sys

from .verify import VERSION, verify


def _load(path: str | None, sequence: str | None, rna_class: str) -> dict:
    if sequence:
        return {"agent": {"sequence": sequence, "rna_class": rna_class}}
    if path is None or path == "-":
        text = sys.stdin.read()
    else:
        with open(path, encoding="utf-8") as handle:
            text = handle.read()
    text = text.strip()
    if not text:
        raise SystemExit("empty input")
    if text[0] in "{[":
        payload = json.loads(text)
        if isinstance(payload, dict) and "agent" not in payload and "sequence" not in payload:
            payload = {"agent": payload}
        return payload
    return {"agent": {"sequence": text, "rna_class": rna_class}}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="rna",
        description="Verify a synthetic RNA Sequence Card.",
    )
    parser.add_argument("command", nargs="?", default="verify", help="verify (default)")
    parser.add_argument("path", nargs="?", help="JSON card, FASTA, or raw sequence. '-' reads stdin.")
    parser.add_argument("--sequence", "-s", help="Raw sequence (skips file input).")
    parser.add_argument("--class", dest="rna_class", default="mRNA", help="RNA class (default mRNA).")
    parser.add_argument("--version", action="version", version=f"rna {VERSION}")
    args = parser.parse_args(argv)
    if args.command not in {"verify", "card"}:
        # allow `rna card.json` without the verify verb
        args.path = args.command if args.path is None else args.path
        args.command = "verify"
    payload = _load(args.path, args.sequence, args.rna_class)
    card = verify(payload)
    json.dump(card, sys.stdout, indent=2)
    sys.stdout.write("\n")
    return 0 if card["verifier"]["verdict"] != "fail" else 2


if __name__ == "__main__":
    raise SystemExit(main())
