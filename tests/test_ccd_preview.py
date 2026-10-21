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
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
TEST C1 C
TEST N1 N
TEST C2 C
TEST C3 C
TEST H1 H
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
TEST N1 H1 SING N
#
"""
class CcdTests(unittest.TestCase):
    def test_order_and_atom_names(self):
        self.assertEqual(m.parse_component(CCD,"TEST"),
                         [["C1","N1",1],["N1","C2",2],["C2","C3",4],["N1","H1",1]])
    def test_wrapped_cif_loop_rows(self):
        wrapped=CCD.replace("TEST N1 C2 DOUB N", "TEST N1\nC2 DOUB N")
        self.assertEqual(len(m.parse_component(wrapped,"TEST")),4)
        self.assertEqual(m.parse_atom_elements(wrapped)["N1"],"N")
    def test_refuses_missing_bond_table(self):
        with self.assertRaises(ValueError):
            m.parse_component("data_TEST\n#\n","TEST")
    def test_cached_component(self):
        with tempfile.TemporaryDirectory() as folder:
            (Path(folder)/"TEST.cif").write_text(CCD)
            with patch.object(m,"urlopen",side_effect=AssertionError("network attempted")):
                bonds, sites, elements=m.get_component("TEST",Path(folder))
                self.assertEqual(len(bonds),4)
                self.assertEqual(sites,{"N1":1})
                self.assertEqual(elements["N1"],"N")
if __name__=="__main__":
    unittest.main()
