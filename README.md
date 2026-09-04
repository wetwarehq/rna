# RNA Lab

Synthetic RNA designed by a model is not yet a molecule. It is four slots that still have to be the cassette they claim to be, and a CDS that still has to satisfy the grammar of translation. RNA Lab stamps that card. The package `rna` (CLI `rn01`) accepts a card, leaves `slots.*.seq` as the agent wrote them, and writes `stamps` beside those slots: identity against packaged `auth.seq` when a cassette is claimed, fail-closed ORF grammar on the CDS, cryptic-ORF excludes on the UTRs, and composition plus fold as recorders that do not decide.

## Abstract

A eukaryotic mRNA that will be decoded by a ribosome is not an arbitrary A/U/G/C polymer. It needs a 5′ untranslated region that does not divert initiation, a Kozak-embedded AUG, an open reading frame whose length is a multiple of three, an in-frame stop, a 3′ UTR free of strong AU-rich instability elements, and usually a poly(A) tail. Innate sensors further constrain uridine density and double-stranded stretches. Those facts belong in the record. They are not identity.

Identity is exact equality of a filled slot with the packaged `auth.seq` for a `claimed_id`. A cassette with the right Kozak substring and the wrong 5′ UTR is not cassette A. Substring search is not identity. The longest AUG-to-stop ORF in a blob is not the payload; that annotator is not used.

This is a sequence-level evaluation. It does not replace in-vitro transcription, cap analysis, integrity electrophoresis, cellular expression, or a folding engine such as ViennaRNA. Reported free energies are pair-count enthalpies, not Turner nearest-neighbour MFE. Structure is recorded. It does not set the stamp.

## Model details

RNA Lab is the Python package `rna`, version 0.2.0. It is a rule-based evaluator for synthetic RNA, with no runtime dependencies, released under MIT. Each call returns a JSON document of schema `rna.card.v1` whose canvas is `slots` plus `stamps`. Packaged identities live in `rna.auth` as `auth/<claimed_id>.fa`. The intended host for codon statistics is *Homo sapiens*; other organisms may be named and are not rescored. There is no `authorised` field.

## Intended use

The primary use is to gate synthetic mRNA designs produced by models that operate on biological sequence, before anyone spends an in-vitro transcription reaction. When `claimed_id` is set, the stamp answers whether the polymer is the named cassette. When it is not set, the stamp answers whether fail-closed ORF grammar cleared. A human operator can submit the same document.

It is not a clinical diagnostic, not a manufacturing release assay, and not a substitute for wet-lab quality control. Selenocysteine, suppressor tRNAs, IRES-driven initiation, programmed frameshifts, circular pairing, host-specific codon tables other than human, capping, double-stranded RNA contamination of an IVT reaction, endotoxin, and protein expression are outside the scope of this version.

## Input

The agent fills `slots.5_utr.seq`, `slots.cds.seq`, `slots.3_utr.seq`, and `slots.polya.seq`. Those strings are stored as submitted. Stamps never overwrite them. Optional `claimed_id` names a packaged cassette. Agent metadata (name, model, class, host, product, declared chemistry, purpose) may sit beside the slots.

A raw sequence or FASTA is accepted as a legacy blob. DNA T is transcribed to U for stamping. A blob does not invent a CDS. If the caller also supplies explicit CDS bounds, the blob may be sliced into slots; there is no ORF search.

Class `mRNA` or `cRNA` invokes CDS grammar on the CDS slot. Any other class skips that grammar. Motifs that must occur, motifs that must not occur, and the default leftover set (T7 promoter, poly(G)₈, poly(C)₉, BsaI, BsmBI) remain fail-closed manufacturing constraints. They are not identity.

Two fixture cards ship in `cards/`: `gold.json` (claim A, polymer A) and `idiot.json` (claim A, 5′ UTR is B).

## Evaluation

Each slot is stamped with `identity: true|false` and a reduction of `match`, `cleared`, or `fail`.

If `claimed_id` is set, identity is true only when the normalised slot equals the packaged `auth.seq` for that cassette and slot. A missing or unknown `claimed_id` file fails identity. Claim A with polymer B fails. Kozak `GCCACCAUG` present in a swapped 5′ UTR does not save it.

Fail-closed grammar, independent of identity: the CDS length is divisible by three, begins with AUG, ends on UAA/UAG/UGA, and contains no in-frame stop before that codon. The 5′ UTR must not contain AUG. The 3′ UTR must not contain an AUG-to-stop ORF. A T7 promoter in the polymer fails. Near-cognate starts are rejected.

Heuristics (codon adaptation, N-degron, uridine load, splice-donor regex, GC, AREs) are recorded and do not set the stamp. 5′ paired fraction and inverted-repeat calls live only in the structure recorder.

## Reported measurements

Every input receives composition and structure under `stamps`. Composition is nucleotide counts and frequencies, GC and AU, uridine fraction, Shannon entropy, dinucleotides, CpG observed/expected, UpA, and homopolymer runs. Structure is a Nussinov maximum matching on A–U, G–C and G–U pairs. Sequences longer than 180 nt are folded as 120-nt windows. These measurements never set the reduction.

## Decision rule

Per slot: `match` if a cassette was claimed and identity holds and fail-closed rules pass; `fail` if identity was claimed and missed, or if fail-closed rules fail; `cleared` if nothing was claimed and fail-closed rules pass.

The card reduction is `fail` if any slot fails or a leftover motif hits; `match` if a cassette was claimed and every slot identity is true; otherwise `cleared`. Informational rows do not change it. `rn01` exits 0 on match or cleared and 2 on fail.

## How to run

There are no runtime dependencies.

```text
uv pip install -e .
rn01 cards/gold.json
rn01 cards/idiot.json
```

```python
from rna import verify
import json

gold = json.load(open("cards/gold.json"))
print(verify(gold)["stamps"]["reduction"])  # match

idiot = json.load(open("cards/idiot.json"))
print(verify(idiot)["stamps"]["reduction"])  # fail
print(verify(idiot)["stamps"]["5_utr"]["identity"])  # False
```

Do not advance a failing card into a wet-lab queue. Do not treat a `cleared` well-formedness stamp as identity. Do not treat approximate kcal or codon adaptation as experimental measurements.

## Limitations

Human codon usage is the only adaptation table in v0.2. Bacterial, yeast, and murine hosts are recorded and not rescored. Circular topology is stored and not given special pairing. The fold is not ViennaRNA. Only cassettes with a file in `rna.auth` can match.

## License

MIT. See `LICENSE`.
