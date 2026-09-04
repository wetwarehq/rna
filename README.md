# RNA Lab

RNA Lab is a verifier for synthetic RNA sequences emitted by agents that operate biological models. An agent fills a Sequence Card with a candidate molecule. The verifier returns the same card with a signed assay: molecular-biology checkers for coding RNA, composition and structure recorders for every class, and optional motif lists. The Python package is named `rna` and is installed with uv.

This document has two audiences. The opening sections are written for a biomedical scientist who needs to know what is being claimed and what is not. The later sections are written for an agentic engineer who needs the card schema, the command line, and the verdict rules.

## Abstract

Synthetic mRNA is a linear ribonucleotide polymer that must be decoded by a ribosome. A well-formed cassette therefore has a constrained 5′ untranslated region, a Kozak-embedded AUG, an open reading frame whose length is a multiple of three, an in-frame stop, a 3′ UTR that does not carry strong instability elements, and usually an encoded or enzymatic poly(A) tail. Innate immune sensors (TLR7/8, RIG-I, MDA5) further constrain uridine density, double-stranded stretches, and GU-rich tracts, which is why therapeutic transcripts are often substituted with N1-methylpseudouridine and why inverted repeats near the 5′ end are a defect rather than a style choice.

Non-coding RNA is not scored against that ORF grammar. For ncRNA the verifier records composition (base counts, GC, dinucleotides, entropy, homopolymers) and an approximate secondary structure (Nussinov pairing, G-quadruplex motifs, 5′ paired fraction). Those recorders run on coding RNA as well. They do not pass or fail a molecule by themselves.

The verifier is a sequence-level assay. It does not replace in-vitro transcription, capping analytics, integrity electrophoresis, expression in cells, or a folding engine such as ViennaRNA. Approximate kcal values are pair-count enthalpies, not Turner nearest-neighbour MFE.

## Sequence Card

A Sequence Card is a JSON object with schema `rna.card.v1`. It has two fills.

The **agent fill** is supplied by the model or the human operator. It names the intended class (`mRNA` / `cRNA` for coding RNA, or an ncRNA class), the host, the product, optional CDS coordinates, claimed chemical modifications, a purpose, and optional whitelist and blacklist motifs.

The **verifier fill** is written by this package. It contains the normalized RNA sequence (whitespace and FASTA headers stripped, `T` transcribed to `U`), resolved CDS coordinates, a translation when coding, a verdict, the checker table, composition, structure, and list hits.

Coordinates are 0-based and half-open on the normalized RNA string. If an mRNA card omits `cds_start` / `cds_end`, the verifier takes the longest AUG-to-stop ORF.

### Agent fill

| Field | Type | Meaning |
| --- | --- | --- |
| `name` | string | Agent or operator identifier. |
| `model` | string | Biological model that emitted the sequence. |
| `sequence` | string | RNA or DNA. FASTA is accepted. |
| `rna_class` | string | `mRNA` or `cRNA` run coding checkers. Any other class is treated as ncRNA. |
| `topology` | string | `linear` (default) or `circular`. |
| `host` | string | Intended expression host. Codon adaptation currently uses human CDS usage. |
| `product` | string | Intended polypeptide or RNA species. |
| `cds_start`, `cds_end` | int | Optional CDS on the normalized sequence. |
| `modifications` | string[] | Declared chemistry, e.g. `m1Psi`, `m7G-cap1`. |
| `purpose` | string | Why this molecule exists. |
| `whitelist` | string[] | Motifs that must occur. A miss is a fail. |
| `blacklist` | string[] or `{id, motif, reason}` | Motifs that must not occur. A hit is a fail. |
| `use_default_blacklist` | bool | Default true. Includes T7 promoter, poly(G)8, poly(C)9, BsaI, BsmBI. |

### Verifier fill

| Field | Meaning |
| --- | --- |
| `verdict` | `pass`, `warn`, or `fail`. Fail if any checker is `fail`. Warn if any checker is `warn` and none fail. |
| `checkers` | Ordered assay rows with `id`, `title`, `status` (`pass` / `warn` / `fail` / `skip` / `info`), `detail`, `evidence`. |
| `composition` | Counts, frequencies, GC, AU, U, Shannon entropy, dinucleotides, CpG O/E, homopolymers. |
| `structure` | Nussinov pairing, approximate kcal, dot-bracket (≤180 nt), G-quadruplexes, 5′ paired fraction, inverted-repeat scan. |
| `lists` | Whitelist and blacklist hits. |
| `translation` | Standard-code protein, `*` for stop. |
| `regions` | 5′ UTR, CDS, 3′ UTR, poly(A) when they can be resolved. |

## Checkers (coding RNA)

These rules run when `rna_class` is `mRNA`, `cRNA`, or `coding`.

**Alphabet and template.** After normalization the polymer must be A/U/G/C. Residual DNA `T` is converted and recorded as a warning. Sequences shorter than 20 nt fail. A T7 promoter (`TAATACGACTCACTATA` / `UAAUACGACUCACUAUA`) in the transcript fails; it is a template leftover, not a 5′ UTR.

**ORF grammar.** The CDS must lie on the sequence, have length divisible by three, begin with AUG, end on UAA/UAG/UGA, and contain no in-frame stop before that terminal codon. Near-cognate starts are rejected for synthetic mRNA. Translation is recorded, not scored.

**Initiation.** A vertebrate Kozak window `gccRccAUGG` is scored. A purine at −3 and G at +4 are expected; a weak window is a warning, not a fail. AUG in the 5′ UTR is reported as a uORF risk. AGGAGG in the 5′ UTR is reported as a bacterial Shine–Dalgarno, unexpected for a eukaryotic host. Transcripts that do not start with G are flagged because T7 initiation and cap addition usually begin at G.

**Composition of the ORF.** Global GC outside 35–70% and CDS GC outside 40–72% warn. Codon adaptation index is computed against human CDS usage (Kazusa CUTG). CAI below 0.80 warns. Sense codons rarer than 8 per thousand warn. Residue 2 in `{R,K,H,F,L,W,Y,I}` warns under the mammalian N-end rule.

**3′ fate.** AAUAAA or AUUAAA in the 3′ UTR passes the polyadenylation-signal check. AUUUA pentamers in the 3′ UTR warn as class I AREs. An encoded poly(A) of ≥80 nt passes; 40–79 nt warns; shorter tails warn unless the operator has declared enzymatic tailing in `modifications`.

**Structure and innate sensing.** Homopolymers of G/C ≥8 or A/U ≥12 outside the tail warn. G-quadruplex motifs in the 5′ UTR or CDS warn because they stall scanning. High CDS CpG or UpA warns. U-fraction ≥22% without a declared pseudouridine/m1Ψ substitution warns (TLR7/8). MAG/GURAGU-like splice donors in the CDS warn if the construct will ever see a nucleus. A reverse complement of the 5′ 24 nt with ≤2 mismatches warns as local dsRNA (RIG-I/MDA5). A Nussinov paired fraction ≥75% in the first 40 nt warns as 5′ structure. UG-rich tracts (`UG`×4) warn.

**Lists.** A missing whitelist motif fails. A blacklist hit fails. The default blacklist is T7 promoter (RNA and DNA), G8, C9, BsaI `GGUCUC`, and BsmBI `CGUCUC`.

## Recorders (all RNA)

Recorders always run. They never set the verdict by themselves.

Composition is nucleotide counts and frequencies, GC and AU fractions, uridine fraction, Shannon entropy on the four-letter alphabet, the sixteen dinucleotides, CpG observed/expected, UpA count, and homopolymer runs.

Structure is a Nussinov maximum matching on A–U, G–C, and G–U pairs with a minimum hairpin loop of three. Sequences longer than 180 nt are folded as 120-nt windows at the 5′ end, middle, and 3′ end. Reported kcal is a sum of mean pair enthalpies plus a coarse loop penalty. G-quadruplexes are the standard `G3+ N1–7` ×4 pattern. The 5′ paired fraction and a 5′ inverted-repeat scan are attached so coding checkers can read them.

## Verdict

`fail` if any checker has status `fail`. `warn` if any checker has status `warn` and none fail. `pass` otherwise. `info` and `skip` rows are documented but do not change the verdict. The process exit code is 0 for pass/warn and 2 for fail, so an agent can gate on the card without parsing every row.

## Install and call

The package has no runtime dependencies.

```text
uv pip install -e .
```

From Python:

```python
from rna import verify

card = verify({
    "agent": {
        "name": "design-agent",
        "model": "example-bio-model",
        "sequence": "GCCACCAUGGCCUAA",
        "rna_class": "mRNA",
        "host": "Homo sapiens",
        "modifications": ["m1Psi"],
        "whitelist": ["GCCACCAUG"],
    }
})
print(card["verifier"]["verdict"])
```

A raw string is accepted and treated as mRNA. FASTA is accepted. DNA letters are transcribed to RNA.

From the shell:

```text
rna verify card.json
rna verify --sequence AUGGCCUAA --class mRNA
printf '{"agent":{"sequence":"AUGGCCUAA"}}' | rna verify -
```

The filled card is written to stdout as JSON.

## Spec for agents

1. Emit a Sequence Card with at least `agent.sequence` and `agent.rna_class`.
2. Call `verify`.
3. Read `verifier.verdict`. Do not ship a `fail` card into a wet-lab queue.
4. Treat `warn` as a design defect to be addressed or explicitly accepted in `agent.notes`.
5. If you designed a CDS, set `cds_start` and `cds_end` yourself. Auto-ORF is a fallback, not a reading-frame oracle.
6. Declare chemistry in `modifications` (`m1Psi`, `m7G-cap1`, enzymatic polyA) so uridine and tail checkers can credit it.
7. Put barcodes, restriction sites, or Kozak windows you require on the whitelist. Put promoter leftovers and toxin motifs on the blacklist.
8. Do not treat approximate kcal, CAI, or G4 calls as experimental measurements.

## Limitations

Human codon usage is the only adaptation table shipped in v0.1. Bacterial, yeast, and murine hosts are recorded on the card but not rescored. Circular topology is stored and not given special pairing. Selenocysteine, suppressor tRNAs, IRES-driven initiation, and programmed frameshifts are out of scope. The fold is not ViennaRNA. The assay does not inspect capping, dsRNA contamination of an IVT reaction, endotoxin, or protein expression.

## License

MIT. See `LICENSE`.
