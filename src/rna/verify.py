"""Sequence Card pipeline: agent fill in, verifier fill out."""

from __future__ import annotations

from datetime import datetime, timezone

from .genetic import longest_orf, translate
from .rules import (
    alphabet_checks,
    coding_checks,
    is_coding,
    list_checks,
    polya_len,
    record_composition,
    record_structure,
    run_lists,
)

VERSION = "0.1.0"
SCHEMA = "rna.card.v1"

CODING_CLASSES = ("mRNA", "cRNA", "coding")


def utcnow() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def normalize(text: str) -> tuple[str, list[str], str]:
    """Strip FASTA, whitespace, and digits; convert T→U. Returns (rna, notes, raw_compact)."""
    notes: list[str] = []
    s = text.strip()
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


def _agent(payload: dict) -> dict:
    src = payload.get("agent") if isinstance(payload.get("agent"), dict) else payload
    seq = src.get("sequence") or src.get("seq") or ""
    rna_class = src.get("rna_class") or src.get("class") or "mRNA"
    whitelist = src.get("whitelist") or []
    blacklist = src.get("blacklist") or []
    if isinstance(whitelist, str):
        whitelist = [ln.strip() for ln in whitelist.splitlines() if ln.strip()]
    if isinstance(blacklist, str):
        blacklist = [ln.strip() for ln in blacklist.splitlines() if ln.strip()]
    mods = src.get("modifications") or []
    if isinstance(mods, str):
        mods = [m.strip() for m in mods.replace(",", " ").split() if m.strip()]
    cds_start = src.get("cds_start")
    cds_end = src.get("cds_end")
    if cds_start is not None:
        cds_start = int(cds_start)
    if cds_end is not None:
        cds_end = int(cds_end)
    default_bl = src.get("use_default_blacklist")
    if default_bl is None:
        default_bl = True
    return {
        "name": src.get("name") or src.get("agent") or "agent",
        "model": src.get("model") or "",
        "sequence": seq,
        "rna_class": rna_class,
        "topology": src.get("topology") or "linear",
        "host": src.get("host") or src.get("organism") or "Homo sapiens",
        "product": src.get("product") or "",
        "cds_start": cds_start,
        "cds_end": cds_end,
        "modifications": mods,
        "purpose": src.get("purpose") or "",
        "notes": src.get("notes") or "",
        "whitelist": whitelist,
        "blacklist": blacklist,
        "use_default_blacklist": bool(default_bl),
    }


def _verdict(checkers: list[dict]) -> str:
    statuses = {c["status"] for c in checkers}
    if "fail" in statuses:
        return "fail"
    if "warn" in statuses:
        return "warn"
    return "pass"


def verify(payload: dict | str) -> dict:
    """Fill a Sequence Card. Accepts a card dict, an agent dict, or a raw sequence string."""
    if isinstance(payload, str):
        payload = {"agent": {"sequence": payload, "rna_class": "mRNA"}}
    agent = _agent(payload)
    rna, notes, raw = normalize(agent["sequence"])
    coding = is_coding(agent["rna_class"])

    cds_start = agent["cds_start"]
    cds_end = agent["cds_end"]
    cds_source = "agent"
    if coding and (cds_start is None or cds_end is None):
        orf = longest_orf(rna)
        if orf:
            cds_start, cds_end = orf
            cds_source = "auto"
            notes.append(f"CDS auto-detected as longest AUG ORF [{cds_start}, {cds_end}).")
        else:
            cds_source = "missing"
            notes.append("No AUG ORF found; CDS checkers will fail.")

    checkers = alphabet_checks(rna, notes, coding)
    translation = ""
    if coding:
        checkers.extend(coding_checks(rna, cds_start, cds_end, agent["modifications"]))
        if cds_start is not None and cds_end is not None and 0 <= cds_start < cds_end <= len(rna):
            translation = translate(rna[cds_start:cds_end])
    else:
        checkers.append(
            {
                "id": "coding_rules",
                "title": "Coding-RNA checkers",
                "status": "skip",
                "detail": f"rna_class={agent['rna_class']}. Molecular-biology ORF checkers apply to mRNA/cRNA only.",
                "evidence": "",
            }
        )

    lists = run_lists(rna, raw, agent["whitelist"], agent["blacklist"], agent["use_default_blacklist"])
    checkers.extend(list_checks(lists))

    composition = record_composition(rna)
    structure = record_structure(rna)
    tail = polya_len(rna)
    verdict = _verdict(checkers)
    counts = {
        "pass": sum(1 for c in checkers if c["status"] == "pass"),
        "warn": sum(1 for c in checkers if c["status"] == "warn"),
        "fail": sum(1 for c in checkers if c["status"] == "fail"),
        "skip": sum(1 for c in checkers if c["status"] == "skip"),
        "info": sum(1 for c in checkers if c["status"] == "info"),
    }

    regions = []
    if cds_start is not None and cds_end is not None and 0 <= cds_start < cds_end <= len(rna):
        if cds_start:
            regions.append({"name": "5utr", "start": 0, "end": cds_start})
        regions.append({"name": "cds", "start": cds_start, "end": cds_end})
        utr3_end = len(rna) - tail if tail >= 15 else len(rna)
        if cds_end < utr3_end:
            regions.append({"name": "3utr", "start": cds_end, "end": utr3_end})
        if tail >= 15:
            regions.append({"name": "polya", "start": len(rna) - tail, "end": len(rna)})
    elif tail >= 15:
        regions.append({"name": "body", "start": 0, "end": len(rna) - tail})
        regions.append({"name": "polya", "start": len(rna) - tail, "end": len(rna)})
    else:
        regions.append({"name": "body", "start": 0, "end": len(rna)})

    verifier = {
        "name": "rna",
        "version": VERSION,
        "timestamp": utcnow(),
        "sequence_normalized": rna,
        "length": len(rna),
        "rna_class": agent["rna_class"],
        "coding": coding,
        "cds_start": cds_start,
        "cds_end": cds_end,
        "cds_source": cds_source,
        "polya_length": tail,
        "regions": regions,
        "translation": translation,
        "verdict": verdict,
        "counts": counts,
        "checkers": checkers,
        "composition": composition,
        "structure": structure,
        "lists": lists,
        "notes": notes,
    }
    return {"schema": SCHEMA, "agent": agent, "verifier": verifier}
