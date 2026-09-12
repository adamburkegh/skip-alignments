"""
Regression coverage for EbiOccurance.trace_probs's dedup mechanism, written
before switching its membership check from an O(n) linear scan over a
plain list (checked_ids) to an O(1) average set lookup -- the same
anti-pattern already fixed in alignment.py's Mapper.node_to_index and
Aligner's closed-set (see CHANGELOG.md), recurring here in the
per-unique-model-path dedup that gates ebi_trace_prob subprocess calls.

trace_probs already deduplicates correctly (two agns sharing the same
model path only trigger one ebi_trace_prob call) -- what's being fixed is
only how that membership check is done, not whether it's done. These
tests assert the call-count behavior directly via mocking ebi_trace_prob,
so they hold before and after the fix.

Run with:
    python -m unittest tests.test_probabilities_trace_probs_dedup -v
"""
import unittest
from concurrent.futures import ThreadPoolExecutor
from unittest.mock import patch

from skipalignments.probabilities import EbiOccurance

# trace_probs submits self.ebi_trace_prob to a ProcessPoolExecutor, which
# must pickle the callable to send to a worker -- a mocked bound method
# (unittest.mock.MagicMock) isn't picklable. Swapping in ThreadPoolExecutor
# for these tests keeps everything in-process (no pickling at all) without
# changing trace_probs' own logic; it's patched at the concurrent.futures
# module level since trace_probs does `from concurrent.futures import
# ProcessPoolExecutor` fresh, inside the method, on every call.
_no_process_pool = patch('concurrent.futures.ProcessPoolExecutor', ThreadPoolExecutor)


class _FakeNode:
    """Minimal stand-in for a ProcessTree leaf: trace_probs only reads .id."""
    def __init__(self, id):
        self.id = id


# Real ProcessTree leaves are shared, identity-stable objects (one per tree
# node) -- trace_counts keys on the node objects themselves (no custom
# __eq__/__hash__ on ProcessTree), so two agns referencing "the same path"
# must reuse the same node instances to behave like real usage. This
# registry mirrors that: same id -> same object, every time.
_NODES = {}


def _node(id):
    if id not in _NODES:
        _NODES[id] = _FakeNode(id)
    return _NODES[id]


def _agn(*ids):
    # A fake agn: a list of (log_move, model_move) pairs, all synchronous
    # (no '>>' model moves).
    return [(i, _node(i)) for i in ids]


def _path(*ids):
    return tuple(_node(i) for i in ids)


class TestTraceProbsDedupesByModelPath(unittest.TestCase):

    def setUp(self):
        self.ebi = EbiOccurance()

    def test_shared_model_path_across_variants_calls_ebi_once(self):
        # 3 agns total, all representing the same model path ('a','b') --
        # only 1 ebi_trace_prob call should ever be dispatched.
        agns = {
            'var1': [_agn('a', 'b'), _agn('a', 'b')],
            'var2': [_agn('a', 'b')],
        }
        with _no_process_pool, patch.object(EbiOccurance, 'ebi_trace_prob', return_value=(0.5, 123)) as mock_call:
            trace_probs_d, trace_counts, total_time = self.ebi.trace_probs(agns, model='fake.slpn')

        self.assertEqual(mock_call.call_count, 1)
        self.assertEqual(trace_probs_d[('a', 'b')], 0.5)

    def test_distinct_model_paths_each_call_ebi_once(self):
        agns = {
            'var1': [_agn('a', 'b'), _agn('a', 'c')],
            'var2': [_agn('a', 'c')],
        }
        with _no_process_pool, patch.object(EbiOccurance, 'ebi_trace_prob', return_value=(0.25, 1)) as mock_call:
            trace_probs_d, trace_counts, total_time = self.ebi.trace_probs(agns, model='fake.slpn')

        self.assertEqual(mock_call.call_count, 2)
        self.assertEqual(set(trace_probs_d.keys()), {('a', 'b'), ('a', 'c')})

    def test_trace_counts_unaffected_by_dedup_mechanism(self):
        # trace_counts is computed in a separate loop, before any dedup --
        # this must stay exactly 2 for ('a','b') within var1 regardless of
        # how the ebi-call dedup below it is implemented.
        agns = {'var1': [_agn('a', 'b'), _agn('a', 'b')]}
        with _no_process_pool, patch.object(EbiOccurance, 'ebi_trace_prob', return_value=(0.5, 1)):
            _, trace_counts, _ = self.ebi.trace_probs(agns, model='fake.slpn')

        self.assertEqual(trace_counts['var1'][_path('a', 'b')], 2)

    def test_many_duplicate_agns_still_call_ebi_once_per_unique_path(self):
        # Larger fan-in, to make an O(n) vs O(1) membership-check distinction
        # meaningful in principle (not timed here, just correctness) -- 50
        # agns, all the same path.
        agns = {'varN': [_agn('a', 'b') for _ in range(50)]}
        with _no_process_pool, patch.object(EbiOccurance, 'ebi_trace_prob', return_value=(0.9, 1)) as mock_call:
            trace_probs_d, trace_counts, _ = self.ebi.trace_probs(agns, model='fake.slpn')

        self.assertEqual(mock_call.call_count, 1)
        self.assertEqual(trace_counts['varN'][_path('a', 'b')], 50)


if __name__ == '__main__':
    unittest.main()
