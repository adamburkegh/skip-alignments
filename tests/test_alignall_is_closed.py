"""
Unit tests for is_closed's dedup correctness after the perf fix that
changed `closed` from a list to a set (see CHANGELOG.md,
tests/test_alignall_search_performance.py): the search's membership
check must consider two states "the same" under exactly the same rule as
before -- equal (alignment, marking) -- and must NOT start treating
structurally-equal-but-object-distinct states as different (which would
silently disable dedup and re-explore states, the opposite failure mode
of the perf fix) or treat genuinely different alignments as the same
(which would incorrectly skip real work).

Uses plain fake state objects (duck-typing __reconstruct_alignment's
`.p`/`.t` walk) rather than real pm4py SearchTuple/net machinery, so this
is decoupled from any actual search -- a fast, precise, isolated pin on
is_closed's own equality semantics.

Run with:
    python -m unittest tests.test_alignall_is_closed -v
"""
import unittest

from skipalignments.alignall import is_closed


class _FakeState:
    """Minimal stand-in for pm4py's SearchTuple: is_closed (via
    __reconstruct_alignment) only ever reads .p, .t, .m, and .g (cost,
    unused by is_closed's own dedup key but read into the returned dict)."""
    def __init__(self, p, t, m):
        self.p = p
        self.t = t
        self.m = m
        self.g = 0


def _chain(*transitions, marking):
    """Builds a linked _FakeState chain (root has no parent/transition)
    ending in a state with the given trailing marking -- matching how
    __reconstruct_alignment walks from a state up to the search root."""
    root = _FakeState(None, None, None)
    state = root
    for t in transitions:
        state = _FakeState(state, t, marking)
    return state


class TestIsClosedDedupSemantics(unittest.TestCase):

    def test_a_state_not_yet_in_closed_is_not_closed(self):
        closed = set()
        state = _chain('t1', 't2', marking='m')
        self.assertFalse(is_closed(state, closed, False))

    def test_an_equal_but_distinct_object_graph_is_recognized_as_closed(self):
        # this is the case the list->set change must not break: two
        # *different* _FakeState chains that reconstruct to the same
        # (alignment, marking) must still dedup, exactly as the old
        # list-based `in` check did
        closed = {(('t1', 't2'), 'm')}
        state = _chain('t1', 't2', marking='m')
        self.assertTrue(is_closed(state, closed, False))

    def test_a_different_alignment_with_the_same_marking_is_not_closed(self):
        closed = {(('t1', 't2'), 'm')}
        state = _chain('t1', 't3', marking='m')
        self.assertFalse(is_closed(state, closed, False))

    def test_the_same_alignment_with_a_different_marking_is_not_closed(self):
        closed = {(('t1', 't2'), 'm')}
        state = _chain('t1', 't2', marking='m2')
        self.assertFalse(is_closed(state, closed, False))

    def test_alignment_order_matters(self):
        # (t1, t2) and (t2, t1) must be treated as distinct -- this is
        # exactly what all_optimal is enumerating ties over
        closed = {(('t1', 't2'), 'm')}
        state = _chain('t2', 't1', marking='m')
        self.assertFalse(is_closed(state, closed, False))


if __name__ == '__main__':
    unittest.main()
