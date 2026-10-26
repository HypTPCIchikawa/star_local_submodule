import unittest
from recenter_survey import selected_tables


class SurveySelectionTests(unittest.TestCase):
    def test_unassessed_slice_and_its_nodes_are_removed_together(self):
        slices = [['a', 'TAG_D', 'd', '88', '1', '2', '1', '2'],
                  ['b', 'TAG_T', 't', '88', '3', '4', '3', '4']]
        nodes = [['a', 'd', '1', '1.2', '.1'], ['b', 't', '1', '3.2', '.2']]
        selected, shapes = selected_tables(slices, nodes, {'a': dict(fitLo=.5, fitHi=2.5)})
        self.assertEqual(len(selected), 1)
        self.assertEqual(shapes, nodes[:1])
        self.assertEqual(selected[0][4:6], [.5, 2.5])
        self.assertEqual(slices[0][4:6], ['1', '2'])

    def test_unknown_or_empty_selection_is_rejected(self):
        for changes in ({}, {'b': dict(fitLo=1, fitHi=2)}):
            with self.assertRaises(ValueError):
                selected_tables([['a']], [['a']], changes)

    def test_missing_nodes_cannot_silently_fit(self):
        with self.assertRaises(ValueError):
            selected_tables([['a', 0, 0, 0, 0, 0]], [], {'a': dict(fitLo=1, fitHi=2)})
