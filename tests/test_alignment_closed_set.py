"""
Regression coverage for the closed-set dedup key used by
Aligner.align_normal_form, written before replacing the O(n) linear scan
(closedset as a list, scanned with State.matches_without_cost per lookup)
with an O(1) average dict lookup keyed on _closed_key -- see CHANGELOG.md
and the process-voids performance investigation this followed.

_closed_key must be an exact hashable equivalent of
State.matches_without_cost: two states are "the same" for closed-set
purposes iff their node-state vectors and remaining traces are equal,
regardless of cost, path, or object identity.

Run with:
    python -m unittest tests.test_alignment_closed_set -v
"""
import unittest

from skipalignments.alignment import Aligner, _closed_key
from skipalignments.processtree import Activity, Sequence, Xor


def _build_tree():
    a = Activity(None, 'a', 100000); a.id = 'a'
    b = Activity(None, 'b', 100000); b.id = 'b'
    c = Activity(None, 'c', 100000); c.id = 'c'
    choice = Xor(None, [b, c]); b.set_parent(choice); c.set_parent(choice)
    tree = Sequence(None, [a, choice]); a.set_parent(tree); choice.set_parent(tree)
    return tree


class TestClosedKeyMatchesWithoutCostSemantics(unittest.TestCase):

    def setUp(self):
        self.tree = _build_tree()
        self.aligner = Aligner(self.tree)

    def _initial_state(self, trace):
        from skipalignments.alignment import State
        return State.initial_state(self.tree, trace, self.aligner.mapper)

    def test_key_is_hashable(self):
        state = self._initial_state(['a', 'b'])
        hash(_closed_key(state))  # must not raise

    def test_two_states_with_identical_state_and_trace_share_a_key(self):
        state1 = self._initial_state(['a', 'b'])
        state2 = self._initial_state(['a', 'b'])
        self.assertTrue(state1.matches_without_cost(state2))
        self.assertEqual(_closed_key(state1), _closed_key(state2))

    def test_states_differing_only_in_cost_or_path_still_share_a_key(self):
        # matches_without_cost is, by name and by its own implementation,
        # indifferent to acc_costs/path -- the dedup key must be too.
        state1 = self._initial_state(['a', 'b'])
        state2 = self._initial_state(['a', 'b'])
        state2.acc_costs = 12345
        state2.path = [('x', 'y')]
        self.assertTrue(state1.matches_without_cost(state2))
        self.assertEqual(_closed_key(state1), _closed_key(state2))

    def test_states_with_different_remaining_trace_have_different_keys(self):
        state1 = self._initial_state(['a', 'b'])
        state2 = self._initial_state(['a', 'c'])
        self.assertFalse(state1.matches_without_cost(state2))
        self.assertNotEqual(_closed_key(state1), _closed_key(state2))

    def test_states_with_different_node_state_vector_have_different_keys(self):
        state1 = self._initial_state(['a', 'b'])
        state2 = self._initial_state(['a', 'b'])
        state2.state = list(state2.state)
        state2.state[0] = state2.state[1]  # perturb the vector
        self.assertFalse(state1.matches_without_cost(state2))
        self.assertNotEqual(_closed_key(state1), _closed_key(state2))


if __name__ == '__main__':
    unittest.main()
