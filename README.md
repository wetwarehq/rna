# RNA Lab

Synthetic RNA designed by a model is not yet a molecule. It is a string that still has to satisfy the grammar of translation, the chemistry of innate sensing, and the housekeeping of an in-vitro transcript. RNA Lab scores that string. The package `rna` accepts a candidate sequence and returns a signed evaluation: whether a coding cassette is well-formed, what the polymer is made of, how it is likely to fold, and whether any required or forbidden motifs are present.

## Abstract

A eukaryotic mRNA that will be decoded by a ribosome is not an arbitrary A/U/G/C polymer. It needs a 5′ untranslated region that does not divert initiation, a Kozak-embedded AUG, an open reading frame whose length is a multiple of three, an in-frame stop, a 3′ UTR free of strong AU-rich instability elements, and usually a poly(A) tail. Innate sensors (TLR7/8, RIG-I, MDA5) further constrain uridine density, double-stranded stretches, and GU-rich tracts. That is why therapeutic transcripts are often substituted with N1-methylpseudouridine, and why a 5′ inverted repeat is a defect rather than a design flourish.

Non-coding RNA is not judged against that open-reading-frame grammar. For every input, including mRNA, the assay still reports composition — base counts, GC, dinucleotides, entropy, homopolymers — and an approximate secondary structure: Nussinov pairing, G-quadruplex motifs, and the paired fraction of the 5′ end. Those measurements describe the polymer. They do not, by themselves, pass or fail it.

This is a sequence-level evaluation. It does not replace in-vitro transcription, cap analysis, integrity electrophoresis, cellular expression, or a folding engine such as ViennaRNA. Reported free energies are pair-count enthalpies, not Turner nearest-neighbour MFE.

## Model details

RNA Lab is the Python package `rna`, version 0.1.0. It is a rule-based evaluator for synthetic RNA, with no runtime dependencies, released under MIT. Each call returns a JSON document of schema `rna.card.v1` that contains the submitted design and the evaluation. The intended host for codon statistics is *Homo sapiens*; other organisms may be named on the document and are not rescored.

## Intended use

The primary use is to gate synthetic mRNA and other RNA designs produced by models that operate on biological sequence, before anyone spends an in-vitro transcription reaction. A human operator can submit the same document. The evaluation is meant for research constructs and therapeutic-style cassettes aimed at human expression.

It is not a clinical diagnostic, not a manufacturing release assay, and not a substitute for wet-lab quality control. Selenocysteine, suppressor tRNAs, IRES-driven initiation, programmed frameshifts, circular pairing, host-specific codon tables other than human, capping, double-stranded RNA contamination of an IVT reaction, endotoxin, and protein expression are outside the scope of this version.

## Input

The caller submits a JSON object, or a raw sequence that is treated as mRNA. FASTA is accepted. DNA T is transcribed to U. Coordinates are 0-based and half-open on the normalised RNA. If a coding design omits CDS bounds, the longest AUG-to-stop ORF is used as a fallback, not as a reading-frame oracle.

The design may name the caller and the biological model that emitted the sequence, the RNA class, topology, intended host, intended product, optional CDS coordinates, declared chemistry such as N1-methylpseudouridine or an m7G cap, and a short purpose. Class `mRNA` or `cRNA` invokes coding rules; any other class is treated as non-coding. Motifs that must occur, motifs that must not occur, and a switch for the default forbidden set (T7 promoter, poly(G)₈, poly(C)₉, BsaI, BsmBI) may be supplied.

## Evaluation

Coding designs are scored against the molecular biology of a eukaryotic cassette.

The polymer must be A/U/G/C after normalisation. Residual T is noted. Sequences shorter than 20 nt fail. A T7 promoter in the transcript fails; it is template leftover, not a 5′ UTR.

The CDS must lie on the sequence, have length divisible by three, begin with AUG, end on UAA, UAG or UGA, and contain no in-frame stop before that codon. Near-cognate starts are rejected for synthetic mRNA. The translation is reported, not scored.

Initiation is scored at the vertebrate Kozak window gccRccAUGG. A weak window warns. AUG in the 5′ UTR is reported as a uORF risk. AGGAGG in the 5′ UTR is a bacterial Shine–Dalgarno, unexpected for a eukaryotic host. Transcripts that do not start with G are flagged because T7 initiation and cap addition usually begin at G.

Codon use is compared with human CDS frequencies (Kazusa CUTG). A codon adaptation index below 0.80 warns. Sense codons rarer than eight per thousand warn. Residue 2 in R, K, H, F, L, W, Y or I warns under the mammalian N-end rule. Global GC outside 35–70% and CDS GC outside 40–72% warn.

In the 3′ UTR, AAUAAA or AUUAAA satisfies the polyadenylation-signal check. AUUUA pentamers warn as class I AREs. An encoded poly(A) of 80 nt or more passes; 40–79 nt warns; shorter tails warn unless enzymatic tailing is declared.

Homopolymers of G or C of length eight or more, or of A or U of length twelve or more outside the tail, warn. G-quadruplex motifs in the 5′ UTR or CDS warn because they can stall scanning ribosomes. High CDS CpG or UpA warns. A uridine fraction of 22% or more without declared pseudouridine or m1Ψ warns as a TLR7/8 ligand. MAG/GURAGU-like splice donors in the CDS warn if the construct will ever see a nucleus. A reverse complement of the first 24 nt with two or fewer mismatches warns as local double-stranded RNA. A Nussinov paired fraction of 75% or more in the first 40 nt warns as 5′ structure. UG-rich tracts warn.

A missing required motif fails. A forbidden motif fails.

## Reported measurements

Every input, coding or not, receives composition and structure. Composition is nucleotide counts and frequencies, GC and AU, uridine fraction, Shannon entropy on the four-letter alphabet, the sixteen dinucleotides, CpG observed/expected, UpA count, and homopolymer runs. Structure is a Nussinov maximum matching on A–U, G–C and G–U pairs with a minimum hairpin loop of three. Sequences longer than 180 nt are folded as 120-nt windows at the 5′ end, the middle, and the 3′ end. Reported kcal is a sum of mean pair enthalpies plus a coarse loop penalty. G-quadruplexes follow the standard G₃₊ N₁–₇ repeated four times.

These measurements never set the decision by themselves.

## Decision rule

The document returns fail if any coding rule or motif constraint fails, warn if any rule warns and none fail, and pass otherwise. Informational and skipped rows are recorded and do not change the decision. The process exits 0 on pass or warn and 2 on fail, so a caller can gate a wet-lab queue without parsing every row.

The returned document includes the normalised sequence, resolved CDS, translation when coding, the decision, the ordered evaluation rows, composition, structure, motif hits, and annotated regions (5′ UTR, CDS, 3′ UTR, poly(A)) when they can be resolved.

## How to run

There are no runtime dependencies.

```text
uv pip install -e .
```

```python
from rna import verify

result = verify({
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
print(result["verifier"]["verdict"])
```

A raw string is treated as mRNA. From the shell:

```text
rna verify card.json
rna verify --sequence AUGGCCUAA --class mRNA
```

Do not advance a failing design into a wet-lab queue. Treat a warning as a defect to fix, or accept it explicitly in the notes. If a CDS was designed, supply its coordinates; automatic ORF finding is a fallback. Declare chemistry so uridine and tail rules can credit it. Do not treat approximate kcal, codon adaptation, or G-quadruplex calls as experimental measurements.

## Limitations

Human codon usage is the only adaptation table in v0.1. Bacterial, yeast, and murine hosts are recorded and not rescored. Circular topology is stored and not given special pairing. The fold is not ViennaRNA.

## License

MIT. See `LICENSE`.
