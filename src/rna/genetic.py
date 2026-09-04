"""Standard genetic code, human codon usage, translation, and CAI."""

from __future__ import annotations

from collections import defaultdict
from math import exp, log

STANDARD = {
    "UUU": "F",
    "UUC": "F",
    "UUA": "L",
    "UUG": "L",
    "UCU": "S",
    "UCC": "S",
    "UCA": "S",
    "UCG": "S",
    "UAU": "Y",
    "UAC": "Y",
    "UAA": "*",
    "UAG": "*",
    "UGU": "C",
    "UGC": "C",
    "UGA": "*",
    "UGG": "W",
    "CUU": "L",
    "CUC": "L",
    "CUA": "L",
    "CUG": "L",
    "CCU": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "CAU": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGU": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "AUU": "I",
    "AUC": "I",
    "AUA": "I",
    "AUG": "M",
    "ACU": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "AAU": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "AGU": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GUU": "V",
    "GUC": "V",
    "GUA": "V",
    "GUG": "V",
    "GCU": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "GAU": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "GGU": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}

STOPS = frozenset({"UAA", "UAG", "UGA"})

# Homo sapiens CDS codon usage, frequency per thousand codons (Kazusa CUTG).
HUMAN_F = {
    "UUU": 17.6,
    "UUC": 20.3,
    "UUA": 7.7,
    "UUG": 12.9,
    "UCU": 15.2,
    "UCC": 17.7,
    "UCA": 12.2,
    "UCG": 4.4,
    "UAU": 12.2,
    "UAC": 15.3,
    "UAA": 1.0,
    "UAG": 0.8,
    "UGU": 10.6,
    "UGC": 12.6,
    "UGA": 1.6,
    "UGG": 13.2,
    "CUU": 13.2,
    "CUC": 19.6,
    "CUA": 7.2,
    "CUG": 39.6,
    "CCU": 17.5,
    "CCC": 19.8,
    "CCA": 16.9,
    "CCG": 6.9,
    "CAU": 10.9,
    "CAC": 15.1,
    "CAA": 12.3,
    "CAG": 34.2,
    "CGU": 4.5,
    "CGC": 10.4,
    "CGA": 6.2,
    "CGG": 11.4,
    "AUU": 16.0,
    "AUC": 20.8,
    "AUA": 7.5,
    "AUG": 22.0,
    "ACU": 13.1,
    "ACC": 18.9,
    "ACA": 15.1,
    "ACG": 6.1,
    "AAU": 17.0,
    "AAC": 19.1,
    "AAA": 24.4,
    "AAG": 31.9,
    "AGU": 12.1,
    "AGC": 19.5,
    "AGA": 12.2,
    "AGG": 12.0,
    "GUU": 11.0,
    "GUC": 14.5,
    "GUA": 7.1,
    "GUG": 28.1,
    "GCU": 18.4,
    "GCC": 27.7,
    "GCA": 15.8,
    "GCG": 7.4,
    "GAU": 21.8,
    "GAC": 25.1,
    "GAA": 29.0,
    "GAG": 39.6,
    "GGU": 10.8,
    "GGC": 22.2,
    "GGA": 16.5,
    "GGG": 16.5,
}

RARE_PER_THOUSAND = 8.0
N_DEGRON = frozenset("RKHFLWYI")


def _codon_weights() -> dict[str, float]:
    by_aa: dict[str, list[tuple[str, float]]] = defaultdict(list)
    for codon, aa in STANDARD.items():
        if aa == "*":
            continue
        by_aa[aa].append((codon, HUMAN_F[codon]))
    weights: dict[str, float] = {}
    for items in by_aa.values():
        peak = max(freq for _, freq in items)
        for codon, freq in items:
            weights[codon] = freq / peak if peak else 0.0
    return weights


WEIGHTS = _codon_weights()


def translate(cds: str) -> str:
    """Translate an RNA CDS with the standard genetic code. Unknown codons become X."""
    residues = []
    for i in range(0, len(cds) - 2, 3):
        residues.append(STANDARD.get(cds[i : i + 3], "X"))
    return "".join(residues)


def cai(cds: str) -> float:
    """Codon adaptation index against human CDS usage. Stop codons are excluded."""
    weights = []
    for i in range(0, len(cds) - 2, 3):
        codon = cds[i : i + 3]
        aa = STANDARD.get(codon)
        if aa is None or aa == "*":
            continue
        weights.append(max(WEIGHTS[codon], 1e-9))
    if not weights:
        return 0.0
    return exp(sum(log(w) for w in weights) / len(weights))


def rare_codons(cds: str) -> list[tuple[int, str, float]]:
    """Return (nucleotide index, codon, frequency per thousand) for rare sense codons."""
    hits = []
    for i in range(0, len(cds) - 2, 3):
        codon = cds[i : i + 3]
        aa = STANDARD.get(codon)
        if aa is None or aa == "*":
            continue
        freq = HUMAN_F.get(codon, 0.0)
        if freq < RARE_PER_THOUSAND:
            hits.append((i, codon, freq))
    return hits


def longest_orf(seq: str) -> tuple[int, int] | None:
    """Longest AUG-to-stop ORF as 0-based half-open coordinates, or None."""
    best: tuple[int, int] | None = None
    n = len(seq)
    for i in range(0, n - 2):
        if seq[i : i + 3] != "AUG":
            continue
        for j in range(i + 3, n - 2, 3):
            if seq[j : j + 3] in STOPS:
                end = j + 3
                if best is None or (end - i) > (best[1] - best[0]):
                    best = (i, end)
                break
    return best


def kozak(seq: str, start: int) -> dict:
    """Score the vertebrate Kozak window gccRccAUGG around CDS start."""
    if start < 6 or start + 4 > len(seq):
        return {
            "window": "",
            "complete": False,
            "minus3_purine": False,
            "plus4_g": False,
            "consensus": False,
        }
    window = seq[start - 6 : start + 4]
    minus3 = window[3] in "AG"  # R at -3
    plus4 = window[9] == "G"
    gcc = window[0:3] == "GCC" or (window[1] == "C" and window[2] == "C")
    return {
        "window": window,
        "complete": True,
        "minus3_purine": minus3,
        "plus4_g": plus4,
        "consensus": minus3 and plus4 and gcc,
    }
