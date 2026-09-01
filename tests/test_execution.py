"""
Unit tests for skipalignments.execution.ExecutionTree's interleaving logic.

These pin down the linearization mechanism that underlies
ExecutionManager.coninciding_agns: coinciding alignments differ only by
reordering concurrent (And) moves, never by padding in extra moves. See the
interleaving operator diamond and its worked example in the paper's
preliminaries.

Run with:
    python -m unittest tests.test_execution -v
"""
import unittest

from skipalignments.execution import ExecutionTree


def _et():
    # all_order_preserving_shuffles/assign don't touch self.execution or
    # self.children, so a bare instance is enough to call them.
    return ExecutionTree(None, [], None)


class TestOrderPreservingShuffles(unittest.TestCase):

    def test_paper_worked_example(self):
        # <a,b> diamond <c> = { <a,b,c>, <a,c,b>, <c,a,b> }
        result = {tuple(r) for r in _et().all_order_preserving_shuffles(['a', 'b'], ['c'])}
        expected = {('a', 'b', 'c'), ('a', 'c', 'b'), ('c', 'a', 'b')}
        self.assertEqual(result, expected)

    def test_two_pairs_interleaving_count_and_order_preservation(self):
        result = _et().all_order_preserving_shuffles(['a', 'b'], ['c', 'd'])
        # binom(4,2) = 6 distinct interleavings
        self.assertEqual(len(result), 6)
        for r in result:
            self.assertEqual([x for x in r if x in ('a', 'b')], ['a', 'b'])
            self.assertEqual([x for x in r if x in ('c', 'd')], ['c', 'd'])
        # all interleavings are distinct
        self.assertEqual(len({tuple(r) for r in result}), 6)

    def test_three_way_interleaving_preserves_each_branch_order(self):
        result = _et().all_order_preserving_shuffles(['a', 'b'], ['c'], ['d', 'e'])
        for r in result:
            self.assertEqual([x for x in r if x in ('a', 'b')], ['a', 'b'])
            self.assertEqual([x for x in r if x in ('d', 'e')], ['d', 'e'])
        # multiset of elements is preserved, only order varies
        for r in result:
            self.assertEqual(sorted(r), ['a', 'b', 'c', 'd', 'e'])

    def test_single_element_paths_produce_only_two_orderings(self):
        result = {tuple(r) for r in _et().all_order_preserving_shuffles(['x'], ['y'])}
        self.assertEqual(result, {('x', 'y'), ('y', 'x')})


if __name__ == '__main__':
    unittest.main()
