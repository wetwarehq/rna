"""Card pipeline: agent fills slots.seq; stamps are written beside them."""

from __future__ import annotations

from datetime import datetime, timezone

from .auth.load import SLOT_NAMES, load_auth
from .rules import (
    alphabet_checks,
    coding_checks,
    is_coding,
    list_checks,
    record_composition,
    record_structure,
    run_lists,
)
from .seq import normalize
from .stamps import cds_translation, polymer_from_slots, stamp_slot

VERSION = "0.2.0"
SCHEMA = "rna.card.v1"

STRUCTURE_CHECKER_IDS = frozenset({"inverted_repeat", "five_prime_structure"})


def utcnow() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _as_list(value) -> list:
    if not value:
        return []
    if isinstance(value, list):
        return value
    if isinstance(value, str):
        return [ln.strip() for ln in value.replace(",", "\n").splitlines() if ln.strip()]
    return []


def _agent_meta(src: dict) -> dict:
    mods = src.get("modifications") or []
    if isinstance(mods, str):
        mods = [m.strip() for m in mods.replace(",", " ").split() if m.strip()]
    return {
        "name": src.get("name") or src.get("agent") or "agent",
        "model": src.get("model") or "",
        "rna_class": src.get("rna_class") or src.get("class") or "mRNA",
        "topology": src.get("topology") or "linear",
        "host": src.get("host") or src.get("organism") or "Homo sapiens",
        "product": src.get("product") or "",
        "modifications": mods,
        "purpose": src.get("purpose") or "",
        "notes": src.get("notes") or "",
        "whitelist": _as_list(src.get("whitelist")),
        "blacklist": src.get("blacklist") or [],
        "use_default_blacklist": True if src.get("use_default_blacklist") is None else bool(src.get("use_default_blacklist")),
    }


def _read_slots(payload: dict, src: dict) -> tuple[dict, list[str], str]:
    """Agent-owned slots. seq is stored as submitted. Never overwritten by stamps."""
    notes: list[str] = []
    raw_slots = payload.get("slots") if isinstance(payload.get("slots"), dict) else src.get("slots")
    slots = {name: {"seq": ""} for name in SLOT_NAMES}
    filled = False
    if isinstance(raw_slots, dict):
        for name in SLOT_NAMES:
            item = raw_slots.get(name)
            if isinstance(item, str):
                slots[name] = {"seq": item}
                filled = filled or bool(item)
            elif isinstance(item, dict):
                seq = item.get("seq") or item.get("sequence") or ""
                slots[name] = {"seq": seq}
                filled = filled or bool(seq)
    if filled:
        return slots, notes, polymer_from_slots(slots)

    blob = src.get("sequence") or src.get("seq") or payload.get("sequence") or ""
    if not blob:
        return slots, notes, ""
    notes.append("Legacy blob fallback. Slots were not filled; CDS is not inferred.")
    rna, n2, _ = normalize(blob)
    notes.extend(n2)
    cds_start = src.get("cds_start")
    cds_end = src.get("cds_end")
    if cds_start is not None and cds_end is not None:
        a, b = int(cds_start), int(cds_end)
        if 0 <= a < b <= len(rna):
            tail = 0
            for ch in reversed(rna):
                if ch != "A":
                    break
                tail += 1
            body_end = len(rna) - tail if tail >= 15 else len(rna)
            slots["5_utr"] = {"seq": rna[:a]}
            slots["cds"] = {"seq": rna[a:b]}
            slots["3_utr"] = {"seq": rna[b:body_end]}
            slots["polya"] = {"seq": rna[body_end:]}
            notes.append("Blob sliced with agent-supplied CDS bounds. No ORF search.")
            return slots, notes, polymer_from_slots(slots)
    return slots, notes, rna


def _card_reduction(slot_stamps: dict, extra_fail: bool, claimed: bool) -> str:
    if extra_fail or any(slot_stamps[n]["reduction"] == "fail" for n in SLOT_NAMES):
        return "fail"
    if claimed and all(slot_stamps[n]["identity"] for n in SLOT_NAMES):
        return "match"
    return "cleared"


def verify(payload: dict | str) -> dict:
    """Stamp a card. Agent seq is not rewritten."""
    if isinstance(payload, str):
        payload = {"sequence": payload, "rna_class": "mRNA"}
    src = payload.get("agent") if isinstance(payload.get("agent"), dict) else payload
    if not isinstance(src, dict):
        src = {}
    agent = _agent_meta(src)
    claimed_id = payload.get("claimed_id") or src.get("claimed_id") or None
    if claimed_id == "":
        claimed_id = None
    slots, notes, polymer = _read_slots(payload, src)
    auth = load_auth(claimed_id)
    claimed = bool(claimed_id)

    slot_stamps = {name: stamp_slot(name, slots[name]["seq"], auth, claimed) for name in SLOT_NAMES}
    coding = is_coding(agent["rna_class"])

    extra_fail = False
    list_rows = []
    if polymer:
        lists = run_lists(
            polymer,
            polymer,
            agent["whitelist"],
            agent["blacklist"],
            agent["use_default_blacklist"],
        )
        list_rows = list_checks(lists)
        extra_fail = extra_fail or any(r["status"] == "fail" for r in list_rows)
    else:
        lists = {"whitelist": [], "blacklist": [], "default_blacklist": agent["use_default_blacklist"]}

    heuristics = []
    cds_rna, _, _ = normalize(slots["cds"]["seq"])
    if coding and cds_rna:
        u5, _, _ = normalize(slots["5_utr"]["seq"])
        start = len(u5)
        end = start + len(cds_rna)
        raw_checks = coding_checks(polymer, start, end, agent["modifications"])
        heuristics = [
            c
            for c in raw_checks
            if c["id"] not in STRUCTURE_CHECKER_IDS and c["status"] in {"warn", "info"}
        ]
    elif not coding:
        heuristics.append(
            {
                "id": "coding_rules",
                "title": "Coding-RNA checkers",
                "status": "skip",
                "detail": f"rna_class={agent['rna_class']}. ORF grammar applies to mRNA/cRNA CDS slots.",
                "evidence": "",
            }
        )

    alpha = alphabet_checks(polymer, notes, coding) if polymer else []
    heuristics.extend(c for c in alpha if c["status"] in {"warn", "info"})

    reduction = _card_reduction(slot_stamps, extra_fail, claimed)
    stamps = {
        **slot_stamps,
        "reduction": reduction,
        "composition": record_composition(polymer),
        "structure": record_structure(polymer),
        "lists": lists,
        "list_checkers": list_rows,
        "heuristics": heuristics,
        "translation": cds_translation(slots),
        "notes": notes,
        "name": "rna",
        "version": VERSION,
        "timestamp": utcnow(),
    }
    return {
        "schema": SCHEMA,
        "claimed_id": claimed_id,
        "slots": slots,
        "stamps": stamps,
        "agent": agent,
    }
