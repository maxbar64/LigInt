import unittest
from pocket_group_candidates import candidates


class GroupingTests(unittest.TestCase):
    def test_adjacent_groups_and_volume_limit(self):
        # 2x2 face-connected blocks: deterministic grouping, no ligand input.
        a=[(x,y,0.0) for x in (0.0,1.0) for y in (0.0,1.0)]
        b=[(x,y,0.0) for x in (2.0,3.0) for y in (0.0,1.0)]
        result=candidates({1:a,2:b},1.0,min_faces=2,min_fraction=.01,max_volume=100)
        self.assertEqual(result["groups"][0]["members"],[1,2])
        self.assertEqual(result["pair_candidates"][0]["members"],[1,2])
        restricted=candidates({1:a,2:b},1.0,min_faces=2,min_fraction=.01,max_volume=7)
        self.assertEqual(restricted["groups"],[])

    def test_disconnected_not_grouped(self):
        a=[(0.0,0.0,0.0)]
        b=[(9.0,0.0,0.0)]
        self.assertEqual(candidates({1:a,2:b},1.0,min_faces=1,min_fraction=0)["groups"],[])


if __name__=="__main__":
    unittest.main()
