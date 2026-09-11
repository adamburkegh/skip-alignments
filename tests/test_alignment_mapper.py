"""
Regression coverage for Mapper.node_to_index, written before switching its
implementation from a linear list.index() scan to the O(1) reverse_lookup
dict that was already being built (but left unused) alongside it -- see
CHANGELOG.md and the process-voids performance investigation this
followed.

Run with:
    python -m unittest tests.test_alignment_mapper -v
"""
import unittest

from skipalignments.alignment import Mapper
from skipalignments.processtree import Activity, Sequence, Xor, And, Loop


def _build_mixed_tree():
    a = Activity(None, 'a', 100000); a.id = 'a'
    b = Activity(None, 'b', 100000); b.id = 'b'
    c = Activity(None, 'c', 100000); c.id = 'c'
    d = Activity(None, 'd', 100000); d.id = 'd'
    e = Activity(None, 'e', 100000); e.id = 'e'
    f = Activity(None, 'f', 100000); f.id = 'f'

    xor_cd = Xor(None, [c, d]); c.set_parent(xor_cd); d.set_parent(xor_cd)
    loop_ef = Loop(None, [e, f]); e.set_parent(loop_ef); f.set_parent(loop_ef)
    and_node = And(None, [xor_cd, loop_ef]); xor_cd.set_parent(and_node); loop_ef.set_parent(and_node)
    tree = Sequence(None, [a, b, and_node]); a.set_parent(tree); b.set_parent(tree); and_node.set_parent(tree)
    return tree, [tree, a, b, and_node, xor_cd, c, d, loop_ef, e, f]


class TestNodeToIndexMatchesTraversalOrder(unittest.TestCase):

    def setUp(self):
        self.tree, self.all_nodes = _build_mixed_tree()
        self.mapper = Mapper(self.tree)

    def test_every_node_has_a_distinct_index(self):
        indices = [self.mapper.node_to_index(n) for n in self.all_nodes]
        self.assertEqual(len(indices), len(set(indices)), "expected every node to map to a distinct index")

    def test_indices_are_a_contiguous_range(self):
        indices = sorted(self.mapper.node_to_index(n) for n in self.all_nodes)
        self.assertEqual(indices, list(range(len(self.all_nodes))))

    def test_index_to_node_round_trips(self):
        for node in self.all_nodes:
            idx = self.mapper.node_to_index(node)
            self.assertIs(self.mapper.index_to_node(idx), node)

    def test_root_is_always_index_zero(self):
        # State.__init__ relies on this: state.state[0] is the root's slot.
        self.assertEqual(self.mapper.node_to_index(self.tree), 0)

    def test_repeated_lookups_are_consistent(self):
        for node in self.all_nodes:
            self.assertEqual(self.mapper.node_to_index(node), self.mapper.node_to_index(node))


if __name__ == '__main__':
    unittest.main()
