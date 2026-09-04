import json
import unittest

from rna.examples import EXAMPLES, PASS_MRNA, UTR5, CDS
from rna.verify import SCHEMA, verify


class VerifyTests(unittest.TestCase):
    def test_pass_mrna_does_not_fail(self):
        card = verify({"agent": EXAMPLES["pass_mRNA"]})
        self.assertEqual(card["schema"], SCHEMA)
        self.assertEqual(card["verifier"]["verdict"], "pass")
        self.assertEqual(card["verifier"]["counts"]["fail"], 0)
        self.assertEqual(card["verifier"]["cds_start"], len(UTR5))
        self.assertEqual(card["verifier"]["translation"].startswith("MASKG"), True)
        self.assertTrue(card["verifier"]["translation"].endswith("*"))
        kozak = next(c for c in card["verifier"]["checkers"] if c["id"] == "kozak")
        self.assertEqual(kozak["status"], "pass")
        start = next(c for c in card["verifier"]["checkers"] if c["id"] == "start_codon")
        self.assertEqual(start["status"], "pass")

    def test_auto_cds(self):
        card = verify({"agent": {"sequence": PASS_MRNA, "rna_class": "mRNA"}})
        self.assertEqual(card["verifier"]["cds_source"], "auto")
        self.assertEqual(card["verifier"]["cds_end"] - card["verifier"]["cds_start"], len(CDS))

    def test_fail_mrna_fails(self):
        card = verify({"agent": EXAMPLES["fail_mRNA"]})
        self.assertEqual(card["verifier"]["verdict"], "fail")
        ids = {c["id"]: c["status"] for c in card["verifier"]["checkers"] if c["status"] == "fail"}
        self.assertIn("t7_promoter", ids)
        self.assertIn("blacklist", ids)
        self.assertIn("whitelist", ids)

    def test_ncrna_skips_coding_checkers(self):
        card = verify({"agent": EXAMPLES["ncRNA_tRNA"]})
        self.assertEqual(card["verifier"]["verdict"], "pass")
        skipped = next(c for c in card["verifier"]["checkers"] if c["id"] == "coding_rules")
        self.assertEqual(skipped["status"], "skip")
        self.assertGreater(card["verifier"]["structure"]["pairs"], 0)
        self.assertIn("gc_fraction", card["verifier"]["composition"])

    def test_premature_stop(self):
        seq = "GCCACCAUGUAAUGCGCCUAA"  # AUG UAA ...
        card = verify({"agent": {"sequence": seq, "rna_class": "mRNA", "use_default_blacklist": False}})
        statuses = {c["id"]: c["status"] for c in card["verifier"]["checkers"]}
        self.assertEqual(statuses.get("start_codon"), "pass")

    def test_custom_blacklist_and_whitelist(self):
        seq = "GGGAAAUGGCCUAA"
        card = verify(
            {
                "agent": {
                    "sequence": seq,
                    "rna_class": "mRNA",
                    "whitelist": ["AUGGCC"],
                    "blacklist": ["CCCC"],
                    "use_default_blacklist": False,
                }
            }
        )
        white = next(c for c in card["verifier"]["checkers"] if c["id"] == "whitelist")
        black = next(c for c in card["verifier"]["checkers"] if c["id"] == "blacklist")
        self.assertEqual(white["status"], "pass")
        self.assertEqual(black["status"], "pass")
        card2 = verify(
            {
                "agent": {
                    "sequence": seq,
                    "rna_class": "mRNA",
                    "whitelist": ["UUUUUU"],
                    "use_default_blacklist": False,
                }
            }
        )
        white2 = next(c for c in card2["verifier"]["checkers"] if c["id"] == "whitelist")
        self.assertEqual(white2["status"], "fail")
        self.assertEqual(card2["verifier"]["verdict"], "fail")

    def test_fasta_and_dna(self):
        fasta = ">mini\nATG GCC TAA\n"
        card = verify(fasta)
        self.assertEqual(card["verifier"]["sequence_normalized"], "AUGGCCUAA")
        notes = " ".join(card["verifier"]["notes"])
        self.assertIn("FASTA", notes)
        self.assertIn("T residue", notes)

    def test_json_roundtrip(self):
        card = verify({"agent": EXAMPLES["pass_mRNA"]})
        blob = json.dumps(card)
        again = json.loads(blob)
        self.assertEqual(again["verifier"]["verdict"], card["verifier"]["verdict"])


if __name__ == "__main__":
    unittest.main()
