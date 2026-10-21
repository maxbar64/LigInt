#!/usr/bin/env python3
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

spec=importlib.util.spec_from_file_location("shared",Path(__file__).resolve().parents[1]/"tools"/"gpcrdb_shared_fingerprints.py")
shared=importlib.util.module_from_spec(spec)
spec.loader.exec_module(shared)

class SharedGpcrTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.d=Path(self.temp.name)
        self.a=self.d/"a.tsv";self.b=self.d/"b.tsv"
        self.a.write_text("chain\tresi\taa\tgpcrdb\tsegment\nR\t113\tD\t3x32\tTM3\nR\t114\tV\t3x33\tTM3\n")
        self.b.write_text("chain\tresi\taa\tgpcrdb\tsegment\nA\t110\tD\t3x32\tTM3\nA\t286\tW\t6x48\tTM6\n")
    def test_intersection_comparable(self):
        a=shared.create([self.a,self.b],self.d/"out")
        self.assertEqual(a["positions"],["3x32"])
        self.assertIn("ASP",(self.d/"out"/"a.pocket.tsv").read_text())
        self.assertIn("3x32",(self.d/"out"/"b.pocket.tsv").read_text())
    def test_union_missing_rejected_and_audited(self):
        with self.assertRaisesRegex(ValueError,"Missing generic"):
            shared.create([self.a,self.b],self.d/"out",mode="union")
        report=json.loads((self.d/"out"/"coverage_audit.json").read_text())
        self.assertEqual(report["positions"],["3x32","3x33","6x48"])
        self.assertEqual(report["structures"][0]["missing"],["6x48"])
    def test_partial_is_not_silently_zeroed(self):
        report=shared.create([self.a,self.b],self.d/"out",mode="union",missing_policy="report")
        self.assertEqual(report["structures"][1]["missing"],["3x33"])
        self.assertNotIn("3x33",(self.d/"out"/"b.pocket.tsv").read_text())

if __name__=="__main__":
    unittest.main()
