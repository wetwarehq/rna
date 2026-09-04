"""Approximate RNA structure recorders. Not a Turner nearest-neighbour MFE."""

from __future__ import annotations

import re

PAIRS = {("A", "U"), ("U", "A"), ("G", "C"), ("C", "G"), ("G", "U"), ("U", "G")}
PAIR_ENERGY = {
    ("A", "U"): -2.1,
    ("U", "A"): -2.1,
    ("G", "C"): -3.4,
    ("C", "G"): -3.4,
    ("G", "U"): -0.9,
    ("U", "G"): -0.9,
}
MIN_LOOP = 3
MAX_FOLD = 180
G4 = re.compile(r"G{3,}[AUGC]{1,7}G{3,}[AUGC]{1,7}G{3,}[AUGC]{1,7}G{3,}")
COMP = str.maketrans("AUGC", "UACG")


def reverse_complement(seq: str) -> str:
    return seq.translate(COMP)[::-1]


def _nussinov(seq: str) -> list[list[int]]:
    n = len(seq)
    dp = [[0] * n for _ in range(n)]
    if n == 0:
        return dp
    for span in range(MIN_LOOP + 1, n):
        for i in range(0, n - span):
            j = i + span
            best = dp[i][j - 1]
            if (seq[i], seq[j]) in PAIRS and (j - i - 1) >= MIN_LOOP:
                inner = dp[i + 1][j - 1] if i + 1 <= j - 1 else 0
                if inner + 1 > best:
                    best = inner + 1
            for k in range(i, j):
                right = dp[k + 1][j] if k + 1 <= j else 0
                cand = dp[i][k] + right
                if cand > best:
                    best = cand
            dp[i][j] = best
    return dp


def _traceback(seq: str, dp: list[list[int]], i: int, j: int, pairs: list[tuple[int, int]]) -> None:
    if i >= j or j - i < MIN_LOOP + 1:
        return
    if dp[i][j] == dp[i][j - 1]:
        _traceback(seq, dp, i, j - 1, pairs)
        return
    if (seq[i], seq[j]) in PAIRS and (j - i - 1) >= MIN_LOOP:
        inner = dp[i + 1][j - 1] if i + 1 <= j - 1 else 0
        if dp[i][j] == inner + 1:
            pairs.append((i, j))
            _traceback(seq, dp, i + 1, j - 1, pairs)
            return
    for k in range(i, j):
        right = dp[k + 1][j] if k + 1 <= j else 0
        if dp[i][j] == dp[i][k] + right:
            _traceback(seq, dp, i, k, pairs)
            _traceback(seq, dp, k + 1, j, pairs)
            return


def fold(seq: str) -> dict:
    """Return pairing stats and a dot-bracket for sequences up to MAX_FOLD nt."""
    n = len(seq)
    if n < 8:
        return {
            "method": "nussinov",
            "folded_nt": n,
            "pairs": 0,
            "paired_fraction": 0.0,
            "approx_kcal": 0.0,
            "kcal_per_nt": 0.0,
            "dot_bracket": "." * n,
            "g_quadruplexes": _g4(seq),
            "note": "Sequence too short for a fold.",
        }
    if n > MAX_FOLD:
        windows = []
        for label, start in (("5p", 0), ("mid", max(0, n // 2 - 60)), ("3p", max(0, n - 120))):
            chunk = seq[start : start + 120]
            windows.append({"region": label, "start": start, **_fold_one(chunk)})
        return {
            "method": "nussinov-windows",
            "folded_nt": 120,
            "pairs": max(w["pairs"] for w in windows),
            "paired_fraction": max(w["paired_fraction"] for w in windows),
            "approx_kcal": min(w["approx_kcal"] for w in windows),
            "kcal_per_nt": min(w["kcal_per_nt"] for w in windows),
            "dot_bracket": "",
            "windows": windows,
            "g_quadruplexes": _g4(seq),
            "note": (
                f"Full fold skipped above {MAX_FOLD} nt. "
                "Windows of 120 nt at the 5' end, middle, and 3' end were folded instead. "
                "This is not a ViennaRNA MFE."
            ),
        }
    result = _fold_one(seq)
    result["g_quadruplexes"] = _g4(seq)
    result["note"] = (
        "Nussinov maximum matching with AU/GC/GU pairs. Approximate kcal uses "
        "mean pair enthalpies, not Turner nearest-neighbour parameters. Use ViennaRNA for MFE."
    )
    return result


def _fold_one(seq: str) -> dict:
    n = len(seq)
    dp = _nussinov(seq)
    pairs: list[tuple[int, int]] = []
    if n:
        _traceback(seq, dp, 0, n - 1, pairs)
    db = ["."] * n
    kcal = 0.0
    for i, j in pairs:
        db[i] = "("
        db[j] = ")"
        kcal += PAIR_ENERGY.get((seq[i], seq[j]), 0.0)
    kcal += 1.3 * max(0, len(pairs) // 6)  # coarse loop penalty
    paired = 2 * len(pairs)
    return {
        "method": "nussinov",
        "folded_nt": n,
        "pairs": len(pairs),
        "paired_fraction": (paired / n) if n else 0.0,
        "approx_kcal": round(kcal, 2),
        "kcal_per_nt": round(kcal / n, 4) if n else 0.0,
        "dot_bracket": "".join(db),
    }


def _g4(seq: str) -> list[dict]:
    hits = []
    for match in G4.finditer(seq):
        hits.append({"start": match.start(), "end": match.end(), "motif": match.group()})
    return hits


def five_prime_paired_fraction(seq: str, window: int = 40) -> float:
    chunk = seq[:window]
    if len(chunk) < 8:
        return 0.0
    folded = _fold_one(chunk)
    return folded["paired_fraction"]


def inverted_repeat(seq: str, seed: int = 24, max_mismatch: int = 3) -> dict | None:
    """Search the body for a reverse complement of the 5' seed (local dsRNA risk)."""
    if len(seq) < seed * 2:
        return None
    probe = reverse_complement(seq[:seed])
    best = None
    limit = len(seq) - seed
    for i in range(seed, limit + 1):
        window = seq[i : i + seed]
        mm = sum(a != b for a, b in zip(probe, window))
        if mm <= max_mismatch and (best is None or mm < best["mismatches"]):
            best = {"start": i, "end": i + seed, "mismatches": mm, "window": window}
            if mm == 0:
                break
    return best
