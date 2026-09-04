"""Fixture sequences used by tests and the laboratory UI."""

from __future__ import annotations

UTR5 = "GGGAGACCCAAGCUGGCUAGCGUUUAAACUUAAGCUUGCCACC"
CDS = "AUGGCCAGCAAGGGCGAGGAGCUGUUCACCGGCGUGGUGCCCAUCCUGGUGGAGCUGGACGGCGACUAA"
UTR3 = "GCUGGCCGGCUUCCUCGCCAAUAAAGCUGGCCUUCG"
POLYA = "A" * 80

PASS_MRNA = UTR5 + CDS + UTR3 + POLYA

FAIL_MRNA = (
    "TAATACGACTCACTATAGGGGGGGGGATGCCCTAATGATTTATTTAGCC"
    "CCCCCCCCCTTTTTTTATGAAATAG"
)

NCRNA_TRNA = (
    "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCACCA"
)

EXAMPLES = {
    "pass_mRNA": {
        "name": "mini-ORF cassette",
        "rna_class": "mRNA",
        "host": "Homo sapiens",
        "product": "MASKGEELFTGVVPILVELDGD",
        "purpose": "Human codon-optimized mini-ORF with Kozak, AAUAAA, and encoded poly(A).",
        "modifications": ["m1Psi", "m7G-cap1"],
        "cds_start": len(UTR5),
        "cds_end": len(UTR5) + len(CDS),
        "sequence": PASS_MRNA,
        "whitelist": ["GCCACCAUG"],
        "blacklist": [],
        "use_default_blacklist": True,
    },
    "fail_mRNA": {
        "name": "broken IVT draft",
        "rna_class": "mRNA",
        "host": "Homo sapiens",
        "product": "",
        "purpose": "Negative control: T7 leftover, polyG, premature stops, ARE.",
        "modifications": [],
        "sequence": FAIL_MRNA,
        "whitelist": ["GCCACCAUG"],
        "blacklist": [],
        "use_default_blacklist": True,
    },
    "ncRNA_tRNA": {
        "name": "yeast tRNA-Phe",
        "rna_class": "ncRNA",
        "host": "Saccharomyces cerevisiae",
        "product": "tRNA-Phe",
        "purpose": "Yeast tRNA-Phe. Composition and fold are measured; CDS grammar is not applied.",
        "modifications": [],
        "sequence": NCRNA_TRNA,
        "whitelist": [],
        "blacklist": [],
        "use_default_blacklist": True,
    },
}
