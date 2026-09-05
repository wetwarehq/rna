import json
import unittest
from pathlib import Path

from rna.auth.load import load_auth
from rna.examples import FAIL_MRNA, NCRNA_TRNA, UTR5, CDS
from rna.stamps import SLOT_NAMES
from rna.verify import SCHEMA, verify

ROOT = Path(__file__).resolve().parents[1]
GOLD = json.loads((ROOT / "cards" / "gold.json").read_text())
IDIOT = json.loads((ROOT / "cards" / "idiot.json").read_text())


class SlotStampTests(unittest.TestCase):
    def test_gold_matches(self):
        card = verify(GOLD)
        self.assertEqual(card["schema"], SCHEMA)
        self.assertEqual(card["claimed_id"], "mini-orf-a")
        self.assertEqual(card["stamps"]["reduction"], "match")
        for name in SLOT_NAMES:
            self.assertEqual(card["slots"][name]["seq"], GOLD["slots"][name]["seq"])
            self.assertTrue(card["stamps"][name]["identity"])
            self.assertEqual(card["stamps"][name]["reduction"], "match")
        self.assertTrue(card["stamps"]["translation"].startswith("MASKG"))
        ids = [c["id"] for name in SLOT_NAMES for c in card["stamps"][name]["checkers"]]
        self.assertNotIn("inverted_repeat", ids)
        self.assertNotIn("five_prime_structure", ids)

    def test_idiot_cassette_swap_fails(self):
        card = verify(IDIOT)
        self.assertEqual(card["stamps"]["reduction"], "fail")
        self.assertFalse(card["stamps"]["5_utr"]["identity"])
        self.assertTrue(card["stamps"]["cds"]["identity"])
        self.assertEqual(card["stamps"]["5_utr"]["reduction"], "fail")
        self.assertEqual(card["slots"]["5_utr"]["seq"], IDIOT["slots"]["5_utr"]["seq"])
        ident = next(c for c in card["stamps"]["5_utr"]["checkers"] if c["id"] == "identity")
        self.assertEqual(ident["status"], "fail")

    def test_stamps_do_not_overwrite_seq(self):
        before = GOLD["slots"]["cds"]["seq"]
        card = verify(GOLD)
        self.assertEqual(card["slots"]["cds"]["seq"], before)
        self.assertNotIn("seq", card["stamps"]["cds"])

    def test_no_orf_search_on_polymer(self):
        card = verify({"sequence": UTR5 + CDS, "rna_class": "mRNA"})
        self.assertEqual(card["slots"]["cds"]["seq"], "")
        self.assertIn("Slots were empty", " ".join(card["stamps"]["notes"]))
        self.assertNotEqual(card["stamps"]["reduction"], "match")

    def test_premature_stop_fails(self):
        card = verify(
            {
                "slots": {
                    "cds": {"seq": "AUGUGCUAAGCCUAA"},
                },
                "agent": {"rna_class": "mRNA", "use_default_blacklist": False},
            }
        )
        statuses = {c["id"]: c["status"] for c in card["stamps"]["cds"]["checkers"]}
        self.assertEqual(statuses["start_codon"], "pass")
        self.assertEqual(statuses["premature_stop"], "fail")
        self.assertEqual(card["stamps"]["cds"]["reduction"], "fail")
        self.assertEqual(card["stamps"]["reduction"], "fail")

    def test_cryptic_orf_excluded_on_utr(self):
        card = verify(
            {
                "slots": {"5_utr": {"seq": "GGGAUGUAAACC"}},
                "agent": {"rna_class": "mRNA", "use_default_blacklist": False},
            }
        )
        cryptic = next(c for c in card["stamps"]["5_utr"]["checkers"] if c["id"] == "cryptic_orf")
        self.assertEqual(cryptic["status"], "fail")
        self.assertEqual(card["stamps"]["reduction"], "fail")

    def test_fail_blob_t7(self):
        card = verify(
            {
                "agent": {
                    "sequence": FAIL_MRNA,
                    "rna_class": "mRNA",
                    "whitelist": ["GCCACCAUG"],
                    "use_default_blacklist": True,
                }
            }
        )
        self.assertEqual(card["stamps"]["reduction"], "fail")

    def test_ncrna_blob_cleared(self):
        card = verify({"agent": {"sequence": NCRNA_TRNA, "rna_class": "ncRNA"}})
        self.assertEqual(card["stamps"]["reduction"], "cleared")
        self.assertGreater(card["stamps"]["structure"]["pairs"], 0)
        self.assertIn("gc_fraction", card["stamps"]["composition"])

    def test_whitelist_is_not_identity(self):
        """Kozak substring present, 5′ UTR swapped, claimed A → fail."""
        card = verify(IDIOT)
        self.assertFalse(card["stamps"]["5_utr"]["identity"])
        polymer = "".join(card["slots"][n]["seq"] for n in SLOT_NAMES)
        self.assertIn("GCCACCAUG", polymer.replace("T", "U"))

    def test_auth_packaged(self):
        auth = load_auth("mini-orf-a")
        self.assertIsNotNone(auth)
        self.assertEqual(auth["5_utr"], GOLD["slots"]["5_utr"]["seq"])
        self.assertIsNone(load_auth("../etc/passwd"))
        self.assertIsNone(load_auth("no-such-cassette"))

    def test_empty_unnamed_fails(self):
        card = verify({"agent": {"rna_class": "mRNA", "use_default_blacklist": False}})
        self.assertEqual(card["stamps"]["reduction"], "fail")
        for name in SLOT_NAMES:
            self.assertEqual(card["stamps"][name]["reduction"], "fail")
            statuses = {c["id"]: c["status"] for c in card["stamps"][name]["checkers"]}
            self.assertEqual(statuses["empty"], "skip")
        self.assertNotEqual(card["stamps"]["reduction"], "cleared")
        self.assertNotIn("skip", [card["stamps"][n]["reduction"] for n in SLOT_NAMES])

    def test_empty_named_fails(self):
        card = verify({"claimed_id": "mini-orf-a", "slots": {}})
        self.assertEqual(card["stamps"]["reduction"], "fail")
        self.assertFalse(card["stamps"]["5_utr"]["identity"])
        self.assertEqual(card["stamps"]["cds"]["reduction"], "fail")

    def test_empty_string_fails(self):
        card = verify("")
        self.assertEqual(card["stamps"]["reduction"], "fail")
        self.assertEqual(card["stamps"]["composition"]["length"], 0)

    def test_fasta_and_dna_blob(self):
        card = verify(">mini\nATG GCC TAA\n")
        self.assertIn("Slots were empty", " ".join(card["stamps"]["notes"]))
        self.assertIn("T residue", " ".join(card["stamps"]["notes"]))


if __name__ == "__main__":
    unittest.main()
