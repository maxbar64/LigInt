import json
import tempfile
import unittest
from pathlib import Path
from pocket_benchmark import coverage, evaluate


def record(serial, name, res, chain, seq, x, y, z):
    return (f"HETATM{serial:5d} {name:4s} {res:3s} {chain}{seq:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{75.0:6.2f}          {'C':>2s}\n")


class BenchmarkTest(unittest.TestCase):
    def test_coverage_uses_xyz_not_pocket_id(self):
        lig = [(("B", "902"), (1.0, 2.0, 3.0))]
        self.assertEqual(coverage(lig, [(1.0, 2.0, 3.0)], 2.0), 1.0)
        self.assertEqual(coverage(lig, [(9.0, 8.0, 7.0)], 2.0), 0.0)

    def test_pdb_json_integration(self):
        with tempfile.TemporaryDirectory() as td:
            base = Path(td)
            crystal = base / "crystal.pdb"
            pocket = base / "pocket.pdb"
            metadata = base / "result.json"
            crystal.write_text(record(1, "C1", "MK1", "B", 902, 1, 2, 3), encoding="utf-8")
            pocket.write_text(record(1, "C", "PKT", "P", 1, 1.5, 2, 3), encoding="utf-8")
            metadata.write_text(json.dumps({"pockets": [{"id": 1, "points": 1,
                "volume_a3": 1.0, "heuristic_score": 0.5}]}), encoding="utf-8")
            result = evaluate(str(crystal), "MK1", str(pocket), str(metadata),
                              "B", "902", 3.0)
            self.assertEqual(result["ligand_heavy_atoms"], 1)
            self.assertEqual(result["best_rank"], 1)
            self.assertEqual(result["top1_coverage"], 1.0)


if __name__ == "__main__":
    unittest.main()
