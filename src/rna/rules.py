"""ORF grammar, composition, structure, and motif lists."""

from __future__ import annotations

import math
import re
from collections import Counter

from .fold import fold, five_prime_paired_fraction, inverted_repeat
from .genetic import (
    STOPS,
    cai,
    kozak,
    rare_codons,
    translate,
)

BASES = "AUGC"
DINUCS = [a + b for a in BASES for b in BASES]
ARE = re.compile(r"AUUUA")
POLYA_SIGNAL = re.compile(r"A[AU]UAAA")
SPLICE_DONOR = re.compile(r"[AG]GGU[AG]AGU")
SHINE_DALGARNO = re.compile(r"AGGAGG")
T7_PROMOTER_RNA = "UAAUACGACUCACUAUA"
T7_PROMOTER_DNA = "TAATACGACTCACTATA"
UG_RICH = re.compile(r"(?:UG){4,}")
PSI_KEYS = ("psi", "pseudour", "m1psi", "n1mpsi", "1-methylpseudo", "n1-methyl")

DEFAULT_BLACKLIST = [
    {
        "id": "T7_promoter_RNA",
        "motif": T7_PROMOTER_RNA,
        "reason": "T7 promoter remnant inside the RNA sequence.",
    },
    {
        "id": "T7_promoter_DNA",
        "motif": T7_PROMOTER_DNA,
        "reason": "T7 promoter DNA leftover (template not fully processed).",
    },
    {
        "id": "polyG8",
        "motif": "GGGGGGGG",
        "reason": "Poly(G) tract: G-quadruplex formation and IVT termination risk.",
    },
    {
        "id": "polyC9",
        "motif": "CCCCCCCCC",
        "reason": "Poly(C) tract: polymerase slippage and unusual structure.",
    },
    {
        "id": "BsaI",
        "motif": "GGUCUC",
        "reason": "BsaI recognition site leftover from a Golden Gate DNA template.",
    },
    {
        "id": "BsmBI",
        "motif": "CGUCUC",
        "reason": "BsmBI recognition site leftover from a Golden Gate DNA template.",
    },
]

CODING_CLASSES = frozenset({"mRNA", "coding"})


def is_coding(rna_class: str) -> bool:
    return rna_class in CODING_CLASSES


def check(id: str, title: str, status: str, detail: str, evidence: str = "") -> dict:
    return {
        "id": id,
        "title": title,
        "status": status,
        "detail": detail,
        "evidence": evidence,
    }


def polya_len(seq: str) -> int:
    n = 0
    for ch in reversed(seq):
        if ch != "A":
            break
        n += 1
    return n


def record_composition(seq: str) -> dict:
    n = len(seq)
    counts = {b: seq.count(b) for b in BASES}
    other = n - sum(counts.values())
    freq = {b: (counts[b] / n if n else 0.0) for b in BASES}
    dinuc = Counter(seq[i : i + 2] for i in range(n - 1))
    dinuc_freq = {d: dinuc.get(d, 0) / (n - 1) if n > 1 else 0.0 for d in DINUCS}
    gc = freq["G"] + freq["C"]
    au = freq["A"] + freq["U"]
    probs = [freq[b] for b in BASES if freq[b] > 0]
    entropy = -sum(p * math.log2(p) for p in probs) if probs else 0.0
    cg = dinuc.get("CG", 0)
    ua = dinuc.get("UA", 0)
    c_count, g_count = counts["C"], counts["G"]
    cpg_oe = (
        (cg / ((c_count * g_count) / n)) if n and c_count and g_count else 0.0
    )
    runs = _homopolymers(seq)
    return {
        "length": n,
        "counts": counts,
        "other": other,
        "frequency": {k: round(v, 4) for k, v in freq.items()},
        "gc_fraction": round(gc, 4),
        "au_fraction": round(au, 4),
        "u_fraction": round(freq["U"], 4),
        "entropy_bits": round(entropy, 4),
        "dinucleotide": {k: round(v, 4) for k, v in dinuc_freq.items()},
        "cpg_count": cg,
        "cpg_oe": round(cpg_oe, 3),
        "upa_count": ua,
        "homopolymers": runs,
    }


def _homopolymers(seq: str) -> dict:
    longest = {b: 0 for b in BASES}
    if not seq:
        return {"longest": longest, "runs": []}
    runs = []
    start = 0
    for i in range(1, len(seq) + 1):
        if i == len(seq) or seq[i] != seq[start]:
            length = i - start
            base = seq[start]
            if base in longest and length > longest[base]:
                longest[base] = length
            if length >= 6:
                runs.append({"base": base, "start": start, "length": length})
            start = i
    return {"longest": longest, "runs": runs}


def record_structure(seq: str) -> dict:
    folded = fold(seq)
    folded["five_prime_paired_fraction"] = round(five_prime_paired_fraction(seq), 4)
    folded["inverted_repeat"] = inverted_repeat(seq)
    return folded


def run_lists(seq: str, raw: str, whitelist: list[str], blacklist: list[dict], use_default: bool) -> dict:
    motifs = list(DEFAULT_BLACKLIST) if use_default else []
    for item in blacklist:
        if isinstance(item, str):
            motifs.append({"id": item, "motif": item, "reason": "Submitted forbidden motif."})
        elif isinstance(item, dict) and item.get("motif"):
            motifs.append(
                {
                    "id": item.get("id") or item["motif"],
                    "motif": item["motif"],
                    "reason": item.get("reason") or "Submitted forbidden motif.",
                }
            )
    black_hits = []
    for item in motifs:
        motif = item["motif"]
        haystack = seq if set(motif.upper()) <= set("AUGC") else raw
        pos = haystack.upper().find(motif.upper())
        if pos >= 0:
            black_hits.append({**item, "start": pos, "hit": True})
    white_hits = []
    for motif in whitelist:
        if not motif:
            continue
        haystack = seq if set(motif.upper()) <= set("AUGC") else raw
        pos = haystack.upper().find(motif.upper())
        white_hits.append(
            {
                "motif": motif,
                "start": pos,
                "hit": pos >= 0,
            }
        )
    return {
        "whitelist": white_hits,
        "blacklist": black_hits,
        "default_blacklist": use_default,
    }


def alphabet_checks(seq: str, notes: list[str], coding: bool = True) -> list[dict]:
    out = []
    if not seq:
        out.append(check("alphabet", "Alphabet", "fail", "No sequence remains after normalization."))
        return out
    bad = [(i, ch) for i, ch in enumerate(seq) if ch not in BASES]
    if bad:
        sample = ", ".join(f"{ch}@{i}" for i, ch in bad[:8])
        out.append(
            check(
                "alphabet",
                "Alphabet",
                "fail",
                f"{len(bad)} non-AUGC residue(s) after T→U normalization: {sample}.",
                sample,
            )
        )
    else:
        out.append(check("alphabet", "Alphabet", "pass", "Sequence is A/U/G/C after normalization."))
    dna_note = next((n for n in notes if n.startswith("Converted")), None)
    if dna_note:
        out.append(check("dna_template", "DNA template residues", "warn", dna_note))
    else:
        out.append(
            check(
                "dna_template",
                "DNA template residues",
                "pass",
                "No T residues were present. Input was already RNA or T-free.",
            )
        )
    n = len(seq)
    if n < 20:
        out.append(check("length", "Length", "fail", f"{n} nt is below a useful 20 nt floor."))
    elif n > 100_000:
        out.append(check("length", "Length", "fail", f"{n} nt exceeds the 100 kb verifier cap."))
    elif n < 60 and coding:
        out.append(check("length", "Length", "warn", f"{n} nt is short for an mRNA cassette; acceptable for ncRNA or mini-ORFs."))
    else:
        out.append(check("length", "Length", "pass", f"{n} nt."))
    if T7_PROMOTER_RNA in seq or T7_PROMOTER_DNA in "".join(notes):
        pass  # leftover T7 is also scored in the motif lists
    if T7_PROMOTER_RNA in seq:
        pos = seq.find(T7_PROMOTER_RNA)
        out.append(
            check(
                "t7_promoter",
                "T7 promoter remnant",
                "fail",
                "T7 promoter (UAAUACGACUCACUAUA) is present in the transcript.",
                f"start={pos}",
            )
        )
    else:
        out.append(check("t7_promoter", "T7 promoter remnant", "pass", "No T7 promoter motif."))
    return out


def coding_checks(seq: str, cds_start: int | None, cds_end: int | None, modifications: list[str]) -> list[dict]:
    out: list[dict] = []
    n = len(seq)
    if cds_start is None or cds_end is None:
        out.append(
            check(
                "cds_bounds",
                "CDS coordinates",
                "fail",
                "No CDS coordinates were supplied.",
            )
        )
        return out
    if not (0 <= cds_start < cds_end <= n):
        out.append(
            check(
                "cds_bounds",
                "CDS coordinates",
                "fail",
                f"CDS [{cds_start}, {cds_end}) is outside the {n} nt sequence.",
            )
        )
        return out
    out.append(
        check(
            "cds_bounds",
            "CDS coordinates",
            "pass",
            f"CDS occupies nucleotides {cds_start}–{cds_end - 1} ({cds_end - cds_start} nt).",
        )
    )
    cds = seq[cds_start:cds_end]
    if len(cds) % 3 != 0:
        out.append(
            check(
                "cds_frame",
                "CDS frame",
                "fail",
                f"CDS length {len(cds)} is not divisible by 3.",
            )
        )
    else:
        out.append(check("cds_frame", "CDS frame", "pass", f"CDS length {len(cds)} nt ({len(cds) // 3} codons)."))

    start = cds[:3] if len(cds) >= 3 else ""
    if start == "AUG":
        out.append(check("start_codon", "Start codon", "pass", "CDS begins with AUG.", start))
    else:
        out.append(
            check(
                "start_codon",
                "Start codon",
                "fail",
                f"CDS begins with {start or '∅'}, not AUG. Near-cognate starts are not accepted for synthetic mRNA.",
                start,
            )
        )

    stop = cds[-3:] if len(cds) >= 3 else ""
    if len(cds) % 3 == 0 and stop in STOPS:
        out.append(check("stop_codon", "Stop codon", "pass", f"In-frame stop {stop} at CDS end.", stop))
    else:
        out.append(
            check(
                "stop_codon",
                "Stop codon",
                "fail",
                f"CDS does not end on an in-frame stop (last codon {stop or '∅'}).",
                stop,
            )
        )

    premature = []
    if len(cds) >= 6 and len(cds) % 3 == 0:
        for i in range(0, len(cds) - 3, 3):
            codon = cds[i : i + 3]
            if codon in STOPS:
                premature.append((cds_start + i, codon))
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
        out.append(check("premature_stop", "Premature stops", "pass", "No in-frame stops inside the ORF body."))

    protein = translate(cds)
    coding_aa = protein[:-1] if protein.endswith("*") else protein
    out.append(
        check(
            "translation",
            "Translation",
            "info",
            f"{len(coding_aa)} amino acids" + (f": {coding_aa}" if len(coding_aa) <= 80 else f": {coding_aa[:80]}…"),
            coding_aa,
        )
    )

    k = kozak(seq, cds_start)
    if not k["complete"]:
        out.append(
            check(
                "kozak",
                "Kozak context",
                "warn",
                "CDS start lacks a 10-nt window (need 6 nt of 5' UTR and 1 nt after AUG).",
            )
        )
    elif k["consensus"]:
        out.append(check("kozak", "Kozak context", "pass", f"Vertebrate Kozak {k['window']}.", k["window"]))
    else:
        missing = []
        if not k["minus3_purine"]:
            missing.append("no purine at −3")
        if not k["plus4_g"]:
            missing.append("no G at +4")
        out.append(
            check(
                "kozak",
                "Kozak context",
                "warn",
                f"Window {k['window']} is a weak Kozak ({', '.join(missing) or 'partial match'}).",
                k["window"],
            )
        )

    utr5 = seq[:cds_start]
    uorfs = [m.start() for m in re.finditer(r"AUG", utr5)]
    if uorfs:
        out.append(
            check(
                "uorf",
                "Upstream AUG",
                "warn",
                f"{len(uorfs)} AUG(s) in the 5' UTR at {uorfs[:8]}. uORFs can suppress the intended start.",
            )
        )
    else:
        out.append(check("uorf", "Upstream AUG", "pass", "No AUG in the 5' UTR."))

    sd = list(SHINE_DALGARNO.finditer(utr5))
    if sd:
        out.append(
            check(
                "shine_dalgarno",
                "Shine–Dalgarno",
                "warn",
                "AGGAGG in the 5' UTR. Unexpected for a eukaryotic mRNA unless a bacterial host is intended.",
            )
        )
    else:
        out.append(check("shine_dalgarno", "Shine–Dalgarno", "pass", "No AGGAGG in the 5' UTR."))

    if seq and seq[0] == "G":
        out.append(
            check(
                "first_nucleotide",
                "5' initiating nucleotide",
                "pass",
                "Transcript starts with G, compatible with T7 initiation and cap addition.",
            )
        )
    else:
        out.append(
            check(
                "first_nucleotide",
                "5' initiating nucleotide",
                "warn",
                f"Transcript starts with {seq[0] if seq else '∅'}. T7 IVT usually initiates at G.",
            )
        )

    gc = (seq.count("G") + seq.count("C")) / n if n else 0
    if gc < 0.35 or gc > 0.70:
        out.append(
            check(
                "gc_content",
                "Global GC",
                "warn",
                f"GC is {gc:.1%}. Therapeutic mRNA is typically 35–70%.",
            )
        )
    else:
        out.append(check("gc_content", "Global GC", "pass", f"GC is {gc:.1%}."))

    cds_gc = (cds.count("G") + cds.count("C")) / len(cds) if cds else 0
    if cds_gc < 0.40 or cds_gc > 0.72:
        out.append(
            check(
                "cds_gc",
                "CDS GC",
                "warn",
                f"CDS GC is {cds_gc:.1%}. Human codon-optimized CDS usually sits near 50–65%.",
            )
        )
    else:
        out.append(check("cds_gc", "CDS GC", "pass", f"CDS GC is {cds_gc:.1%}."))

    cai_value = cai(cds)
    if cai_value >= 0.80:
        out.append(check("cai_human", "Human CAI", "pass", f"CAI {cai_value:.3f} against human CDS usage."))
    elif cai_value >= 0.65:
        out.append(check("cai_human", "Human CAI", "warn", f"CAI {cai_value:.3f} is moderate. Consider further codon optimization for human expression."))
    else:
        out.append(check("cai_human", "Human CAI", "warn", f"CAI {cai_value:.3f} is low for a human-intended CDS."))

    rares = rare_codons(cds)
    if rares:
        sample = ", ".join(f"{c}@{cds_start + i} ({f:.1f}/1000)" for i, c, f in rares[:8])
        out.append(
            check(
                "rare_codons",
                "Rare human codons",
                "warn",
                f"{len(rares)} sense codon(s) below 8 per thousand in human CDS: {sample}.",
                sample,
            )
        )
    else:
        out.append(check("rare_codons", "Rare human codons", "pass", "No rare human sense codons."))

    if len(coding_aa) >= 2 and coding_aa[1] in "RKHFLWYI":
        out.append(
            check(
                "n_degron",
                "N-degron",
                "warn",
                f"Residue 2 is {coding_aa[1]}. Mammalian N-end rule classifies R/K/H/F/L/W/Y/I as destabilizing after Met processing.",
                coding_aa[:4],
            )
        )
    elif len(coding_aa) >= 2:
        out.append(
            check(
                "n_degron",
                "N-degron",
                "pass",
                f"Residue 2 is {coding_aa[1]}, not a primary mammalian N-degron residue.",
                coding_aa[:4],
            )
        )
    else:
        out.append(check("n_degron", "N-degron", "skip", "Peptide too short to score residue 2."))

    tail = polya_len(seq)
    utr3_end = n - tail if tail >= 15 else n
    utr3 = seq[cds_end:utr3_end]
    signals = list(POLYA_SIGNAL.finditer(utr3))
    if signals:
        pos = cds_end + signals[0].start()
        out.append(
            check(
                "polya_signal",
                "Polyadenylation signal",
                "pass",
                f"{signals[0].group()} in the 3' UTR at {pos}.",
                signals[0].group(),
            )
        )
    else:
        if tail >= 40:
            out.append(
                check(
                    "polya_signal",
                    "Polyadenylation signal",
                    "info",
                    "No AAUAAA/AUUAAA in the 3' UTR. Acceptable when a poly(A) tail is encoded on the IVT template.",
                )
            )
        else:
            out.append(
                check(
                    "polya_signal",
                    "Polyadenylation signal",
                    "warn",
                    "No AAUAAA/AUUAAA in the 3' UTR and no long encoded poly(A) tail.",
                )
            )

    ares = list(ARE.finditer(utr3))
    if ares:
        pos = [cds_end + m.start() for m in ares]
        out.append(
            check(
                "are",
                "AU-rich instability",
                "warn",
                f"{len(ares)} AUUUA pentamer(s) in the 3' UTR at {pos[:8]}. Class I AREs destablize mRNA via TTP/HuR.",
            )
        )
    else:
        out.append(check("are", "AU-rich instability", "pass", "No AUUUA in the 3' UTR."))

    if tail >= 80:
        out.append(check("polya_tail", "Poly(A) tail", "pass", f"Encoded 3' poly(A) of {tail} nt."))
    elif tail >= 40:
        out.append(check("polya_tail", "Poly(A) tail", "warn", f"Encoded poly(A) is {tail} nt; IVT mRNA is often 80–120 A."))
    else:
        out.append(
            check(
                "polya_tail",
                "Poly(A) tail",
                "warn",
                f"Only {tail} trailing A residues. If the tail is added enzymatically, record that in modifications.",
            )
        )

    body_end = n - tail if tail >= 15 else n
    body = seq[:body_end]
    runs = _homopolymers(body)["longest"]
    bad_runs = []
    if runs["G"] >= 8:
        bad_runs.append(f"G×{runs['G']}")
    if runs["C"] >= 8:
        bad_runs.append(f"C×{runs['C']}")
    if runs["U"] >= 12:
        bad_runs.append(f"U×{runs['U']}")
    if runs["A"] >= 12:
        bad_runs.append(f"A×{runs['A']}")
    if bad_runs:
        out.append(
            check(
                "homopolymer",
                "Homopolymers",
                "warn",
                "Long homopolymers outside the poly(A) tail: " + ", ".join(bad_runs) + ".",
            )
        )
    else:
        out.append(check("homopolymer", "Homopolymers", "pass", "No long G/C/U/A runs outside the poly(A) tail."))

    g4 = [m for m in re.finditer(r"G{3,}[AUGC]{1,7}G{3,}[AUGC]{1,7}G{3,}[AUGC]{1,7}G{3,}", seq[:cds_end])]
    if g4:
        out.append(
            check(
                "g_quadruplex",
                "G-quadruplex",
                "warn",
                f"{len(g4)} G4-like motif(s) in the 5' UTR or CDS. These can stall scanning ribosomes.",
                g4[0].group(),
            )
        )
    else:
        out.append(check("g_quadruplex", "G-quadruplex", "pass", "No G4 motif in the 5' UTR or CDS."))

    cg = sum(1 for i in range(len(cds) - 1) if cds[i : i + 2] == "CG")
    ua = sum(1 for i in range(len(cds) - 1) if cds[i : i + 2] == "UA")
    cpg_note = f"CDS CpG={cg}, UpA={ua}."
    if cg / max(len(cds) - 1, 1) > 0.08:
        out.append(check("cpg_upa", "CpG / UpA load", "warn", cpg_note + " High CpG can increase innate sensing of unmodified RNA."))
    elif ua / max(len(cds) - 1, 1) > 0.10:
        out.append(check("cpg_upa", "CpG / UpA load", "warn", cpg_note + " High UpA is associated with faster mRNA turnover."))
    else:
        out.append(check("cpg_upa", "CpG / UpA load", "pass", cpg_note))

    u_frac = seq.count("U") / n if n else 0
    mods = " ".join(modifications).lower()
    has_psi = any(k in mods for k in PSI_KEYS)
    if u_frac >= 0.22 and not has_psi:
        out.append(
            check(
                "uridine_load",
                "Uridine load",
                "warn",
                f"U is {u_frac:.1%} and no pseudouridine/m1Ψ is declared. Unmodified U-rich RNA stimulates TLR7/8.",
            )
        )
    elif has_psi:
        out.append(
            check(
                "uridine_load",
                "Uridine load",
                "pass",
                f"U is {u_frac:.1%}. Modified uridine is declared ({', '.join(modifications)}).",
            )
        )
    else:
        out.append(check("uridine_load", "Uridine load", "pass", f"U is {u_frac:.1%}."))

    donors = list(SPLICE_DONOR.finditer(cds))
    if donors:
        pos = cds_start + donors[0].start()
        out.append(
            check(
                "cryptic_splice",
                "Cryptic splice donor",
                "warn",
                f"MAG/GURAGU-like donor at {pos}. Relevant if the construct is ever RNA-spliced (DNA plasmid, nuclear RNA).",
            )
        )
    else:
        out.append(check("cryptic_splice", "Cryptic splice donor", "pass", "No MAG/GURAGU donor in the CDS."))

    alts = []
    if len(cds) % 3 == 0:
        for i in range(3, len(cds) - 3, 3):
            continue
        for i in range(1, len(cds) - 2):
            if i % 3 == 0:
                continue
            if cds[i : i + 3] == "AUG":
                alts.append(cds_start + i)
    if alts:
        out.append(
            check(
                "out_of_frame_aug",
                "Out-of-frame AUG",
                "info",
                f"{len(alts)} out-of-frame AUG(s) in the CDS at {alts[:8]}. Relevant under leaky scanning.",
            )
        )
    else:
        out.append(check("out_of_frame_aug", "Out-of-frame AUG", "pass", "No out-of-frame AUG in the CDS."))

    tandem = seq[cds_end : cds_end + 3]
    if tandem in STOPS:
        out.append(check("tandem_stop", "Tandem stop", "info", f"Additional stop {tandem} immediately after the CDS."))
    else:
        out.append(check("tandem_stop", "Tandem stop", "info", "No immediate tandem stop. Optional for IVT mRNA."))

    ug = list(UG_RICH.finditer(seq))
    if ug:
        out.append(
            check(
                "ug_rich",
                "UG-rich innate motif",
                "warn",
                f"{len(ug)} (UG)≥4 tract(s). GU-rich RNA is a TLR7/8 ligand when unmodified.",
            )
        )
    else:
        out.append(check("ug_rich", "UG-rich innate motif", "pass", "No (UG)4 tract."))

    return out


def list_checks(lists: dict) -> list[dict]:
    out = []
    misses = [w for w in lists["whitelist"] if not w["hit"]]
    if not lists["whitelist"]:
        out.append(check("whitelist", "Required motifs", "skip", "No required motifs were supplied."))
    elif misses:
        motifs = ", ".join(m["motif"] for m in misses)
        out.append(check("whitelist", "Required motifs", "fail", f"Required motif(s) missing: {motifs}."))
    else:
        out.append(check("whitelist", "Required motifs", "pass", f"All {len(lists['whitelist'])} required motif(s) are present."))
    if lists["blacklist"]:
        motifs = ", ".join(f"{h['id']}@{h['start']}" for h in lists["blacklist"])
        out.append(
            check(
                "blacklist",
                "Forbidden motifs",
                "fail",
                f"{len(lists['blacklist'])} forbidden motif(s): {motifs}.",
                motifs,
            )
        )
    else:
        src = "default leftover set" if lists["default_blacklist"] else "submitted forbidden list"
        out.append(check("blacklist", "Forbidden motifs", "pass", f"No hits against the {src}."))
    return out
