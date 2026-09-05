"""Slot stamps. Identity against packaged auth.seq. Structure is not consulted."""

from __future__ import annotations

from .genetic import STOPS, translate
from .seq import normalize

SLOT_NAMES = ("5_utr", "cds", "3_utr", "polya")


def check(id: str, title: str, status: str, detail: str, evidence: str = "") -> dict:
    return {"id": id, "title": title, "status": status, "detail": detail, "evidence": evidence}


def cryptic_orfs(seq: str) -> list[tuple[int, int]]:
    """AUG-to-in-frame-stop spans of at least one sense codon plus stop."""
    hits: list[tuple[int, int]] = []
    n = len(seq)
    for i in range(0, n - 5):
        if seq[i : i + 3] != "AUG":
            continue
        for j in range(i + 3, n - 2, 3):
            if seq[j : j + 3] in STOPS:
                hits.append((i, j + 3))
                break
    return hits


def _alphabet(seq: str) -> list[dict]:
    if not seq:
        return [check("alphabet", "Alphabet", "fail", "Empty slot.")]
    bad = [(i, ch) for i, ch in enumerate(seq) if ch not in "AUGC"]
    if bad:
        sample = ", ".join(f"{ch}@{i}" for i, ch in bad[:8])
        return [check("alphabet", "Alphabet", "fail", f"Non-AUGC residues: {sample}.", sample)]
    return [check("alphabet", "Alphabet", "pass", "A/U/G/C.")]


def slot_grammar(name: str, rna: str) -> list[dict]:
    """Required grammar for one slot. No fold, no CAI, no motif lists."""
    if not rna:
        return [check("empty", "Empty slot", "skip", "No sequence in this slot.")]
    out = _alphabet(rna)
    if name == "cds":
        out.extend(_cds_grammar(rna))
    elif name in {"5_utr", "3_utr"}:
        out.extend(_utr_cryptic(name, rna))
    elif name == "polya":
        if rna and set(rna) <= {"A"}:
            out.append(check("polya_homopolymer", "Poly(A)", "pass", f"{len(rna)} trailing A."))
        elif rna:
            out.append(
                check(
                    "polya_homopolymer",
                    "Poly(A)",
                    "fail",
                    "Poly(A) slot contains non-A residues.",
                )
            )
    return out


def _cds_grammar(cds: str) -> list[dict]:
    out = []
    if len(cds) % 3 != 0:
        out.append(check("cds_frame", "CDS frame", "fail", f"CDS length {len(cds)} is not divisible by 3."))
    else:
        out.append(check("cds_frame", "CDS frame", "pass", f"{len(cds)} nt ({len(cds) // 3} codons)."))
    start = cds[:3] if len(cds) >= 3 else ""
    if start == "AUG":
        out.append(check("start_codon", "Start codon", "pass", "CDS begins with AUG.", start))
    else:
        out.append(check("start_codon", "Start codon", "fail", f"CDS begins with {start or '∅'}, not AUG.", start))
    stop = cds[-3:] if len(cds) >= 3 else ""
    if len(cds) % 3 == 0 and stop in STOPS:
        out.append(check("stop_codon", "Stop codon", "pass", f"In-frame stop {stop}.", stop))
    else:
        out.append(check("stop_codon", "Stop codon", "fail", f"CDS does not end on an in-frame stop ({stop or '∅'}).", stop))
    premature = []
    if len(cds) >= 6 and len(cds) % 3 == 0:
        for i in range(0, len(cds) - 3, 3):
            codon = cds[i : i + 3]
            if codon in STOPS:
                premature.append((i, codon))
    if premature:
        sample = ", ".join(f"{c}@{p}" for p, c in premature[:6])
        out.append(
            check(
                "premature_stop",
                "Premature stops",
                "fail",
                f"{len(premature)} in-frame stop(s) before the terminal codon: {sample}.",
                sample,
            )
        )
    else:
        out.append(check("premature_stop", "Premature stops", "pass", "No in-frame stops in the ORF body."))
    return out


def _utr_cryptic(name: str, rna: str) -> list[dict]:
    if name == "5_utr":
        augs = [i for i in range(len(rna) - 2) if rna[i : i + 3] == "AUG"]
        if augs:
            return [
                check(
                    "cryptic_orf",
                    "Cryptic ORF",
                    "fail",
                    f"5′ UTR contains AUG at {augs[:8]}. Excluded; this is not a payload annotation.",
                )
            ]
        return [check("cryptic_orf", "Cryptic ORF", "pass", "No AUG in the 5′ UTR.")]
    hits = cryptic_orfs(rna)
    if hits:
        sample = ", ".join(f"[{a},{b})" for a, b in hits[:6])
        return [
            check(
                "cryptic_orf",
                "Cryptic ORF",
                "fail",
                f"3′ UTR contains AUG-to-stop ORF(s) {sample}. Excluded.",
                sample,
            )
        ]
    return [check("cryptic_orf", "Cryptic ORF", "pass", "No cryptic ORF in the 3′ UTR.")]


def slot_reduction(claimed: bool, identity: bool, checkers: list[dict], rna: str) -> str:
    if not rna:
        return "fail"
    any_fail = any(c["status"] == "fail" for c in checkers)
    if claimed:
        if identity and not any_fail:
            return "match"
        return "fail"
    if any_fail:
        return "fail"
    return "cleared"


def stamp_slot(name: str, raw_seq: str, auth: dict[str, str] | None, claimed: bool) -> dict:
    rna, notes, _ = normalize(raw_seq)
    identity = False
    if claimed and auth is not None and name in auth:
        identity = rna == auth[name]
    checkers = slot_grammar(name, rna)
    if claimed and auth is None:
        checkers.append(
            check(
                "identity",
                "Identity",
                "fail",
                "claimed_id has no packaged auth.seq.",
            )
        )
    elif claimed and not identity:
        checkers.append(
            check(
                "identity",
                "Identity",
                "fail",
                f"{name} does not equal packaged auth.seq for the claimed cassette.",
            )
        )
    elif claimed and identity:
        checkers.append(check("identity", "Identity", "pass", f"{name} equals packaged auth.seq."))
    return {
        "identity": identity,
        "reduction": slot_reduction(claimed, identity, checkers, rna),
        "length": len(rna),
        "checkers": checkers,
        "notes": notes,
    }


def polymer_from_slots(slots: dict) -> str:
    parts = []
    for name in SLOT_NAMES:
        rna, _, _ = normalize(slots.get(name, {}).get("seq") or "")
        parts.append(rna)
    return "".join(parts)


def cds_translation(slots: dict) -> str:
    rna, _, _ = normalize(slots.get("cds", {}).get("seq") or "")
    return translate(rna) if rna else ""
