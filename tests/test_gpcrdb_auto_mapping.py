#!/usr/bin/env python3
"""Offline regression tests for GPCRdb mapping (no network required)."""
import importlib.util
import json
import pathlib
import tempfile
import unittest
from unittest.mock import patch
from io import BytesIO

SCRIPT = pathlib.Path(__file__).resolve().parents[1] / "tools" / "gpcrdb_auto_mapping.py"
spec = importlib.util.spec_from_file_location("gpcrdb_auto_mapping", SCRIPT)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)

class AutoGpcrdbTests(unittest.TestCase):
    def test_mapping_from_cached_reference(self):
        with tempfile.TemporaryDirectory() as temp:
            directory = pathlib.Path(temp)
            seq = "ARNDCEQGHILKMFPSTWYV"*4
            reverse = {aa: code for code, aa in module.AA.items() if len(code)==3}
            pdb = directory/"example.pdb"
            pdb.write_text("".join(
                f"ATOM  {i+1:5d}   CA {reverse[aa]:>3s} R{i+101:4d}    {float(i):8.3f}{0.:8.3f}{0.:8.3f}{1.:6.2f}{20.:6.2f}           C\n"
                for i,aa in enumerate(seq)
            )+"END\n")
            cache = directory/"cache"
            cache.mkdir()
            ref=[{"sequence_number":i+1,"amino_acid":aa,
                  "display_generic_number":f"3.{i%80+10:02d}x{i%80+10:02d}",
                  "protein_segment":"TM3"} for i,aa in enumerate(seq)]
            (cache/"adrb2_human.residues.json").write_text(json.dumps(ref))
            out, audit = directory/"mapping.tsv", directory/"audit.json"
            report = module.create(pdb,"R","adrb2_human",cache,out,audit)
            self.assertEqual(report["exact_matches"],80)
            self.assertEqual(report["numbered_positions"],80)
            self.assertEqual(len(out.read_text().splitlines()),81)
            self.assertEqual(json.loads(audit.read_text())["coverage"],1.0)

    def test_detect_receptor_chain(self):
        with tempfile.TemporaryDirectory() as temp:
            root = pathlib.Path(temp)
            seq = "ARNDCEQGHILKMFPSTWYV"*4
            reverse = {aa: code for code, aa in module.AA.items() if len(code)==3}
            lines = []
            for i, aa in enumerate(seq):
                lines.append(f"ATOM  {i+1:5d}   CA {reverse[aa]:>3s} R{i+101:4d}    {float(i):8.3f}{0.:8.3f}{0.:8.3f}{1.:6.2f}{20.:6.2f}           C\\n")
            for i, aa in enumerate("G"*50):
                lines.append(f"ATOM  {100+i:5d}   CA {reverse[aa]:>3s} A{i+1:4d}    {float(i):8.3f}{0.:8.3f}{0.:8.3f}{1.:6.2f}{20.:6.2f}           C\\n")
            pdb = root/"complex.pdb"
            pdb.write_text("".join(lines))
            reference = [(i+1, aa, "3x32", "TM3") for i,aa in enumerate(seq)]
            self.assertEqual(module.detect_chain(pdb, reference), "R")

    def test_pdb_receptor_lookup(self):
        with patch.object(module, "urlopen", return_value=BytesIO(
            json.dumps({"protein": {"entry_name": "adrb2_human"}}).encode())):
            self.assertEqual(module.identify_entry_from_pdb("3SN6"), "adrb2_human")

    def test_pdb_lookup_rejects_unrecognized_schema(self):
        with patch.object(module, "urlopen", return_value=BytesIO(b'{"protein":"unknown-value"}')):
            self.assertEqual(module.identify_entry_from_pdb("3SN6"), "unknown-value")
        with patch.object(module, "urlopen", return_value=BytesIO(b'{"result":"missing"}')):
            with self.assertRaises(ValueError):
                module.identify_entry_from_pdb("3SN6")

    def test_reject_mismatched_reference(self):
        with self.assertRaises(ValueError):
            module.align("A"*3000, "A"*3000)

if __name__ == "__main__":
    unittest.main()
