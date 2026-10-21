#!/usr/bin/env python3
"""Offline CCD bond parser checks; no remote traffic."""
import importlib.util
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

path=Path(__file__).resolve().parents[1]/"tools"/"ccd_preview.py"
spec=importlib.util.spec_from_file_location("ccd_preview",path)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)

CCD="""data_TEST
#
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
_chem_comp_bond.pdbx_aromatic_flag
TEST C1 N1 SING N
TEST N1 C2 DOUB N
TEST C2 C3 SING Y
#
"""
class CcdTests(unittest.TestCase):
    def test_order_and_atom_names(self):
        self.assertEqual(m.parse_component(CCD,"TEST"),
                         [["C1","N1",1],["N1","C2",2],["C2","C3",4]])
    def test_refuses_missing_bond_table(self):
        with self.assertRaises(ValueError):
            m.parse_component("data_TEST\n#\n","TEST")
    def test_cached_component(self):
        with tempfile.TemporaryDirectory() as folder:
            (Path(folder)/"TEST.cif").write_text(CCD)
            with patch.object(m,"urlopen",side_effect=AssertionError("network attempted")):
                self.assertEqual(len(m.get_component("TEST",Path(folder))),3)
if __name__=="__main__":
    unittest.main()
