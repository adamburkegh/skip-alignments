"""
Regression test for __search's tie-count scaling: __reconstruct_alignment
rebuilt its alignment-so-far list via repeated `[parent.t] + alignment`
while walking up to the root -- O(depth) per step, O(depth) steps, so
O(depth^2) per single reconstruction -- and `closed` was a plain list, so
`is_closed`'s `search_tuple in closed` was an O(n) linear scan. Together,
on a real log with thousands of tied optimal alignments (see
test_alignall_tie_explosion.py's rtfm.xes report), this was the actual
dominant cost, not anything inherent to the search itself: n=7 (5040
ties) on the fixture below took 43s before the fix, 1.36s after, with an
*identical* tie count both times (this is a pure performance fix, zero
change in output -- see tests/perf_search_scaling.py for the full
before/after numbers up to n=8/40320 ties, too slow to run as part of the
normal suite).

This test uses n=7 specifically because it's the smallest case in that
scaling probe where the pre-fix runtime (43s) would have failed CI/been
obviously wrong, while the post-fix runtime (~1-2s) comfortably fits a
generous ceiling -- large enough to catch a reintroduced O(depth^2) or
O(n) regression, small enough to run in the normal suite.

Run with:
    python -m unittest tests.test_alignall_search_performance -v
"""
import math
import time
import unittest

from skipalignments.alignall import align_pn_all
from skipalignments.probabilities import EbiOccurance

from tests.test_alignall_tie_explosion import _tree_with_n_concurrent_optional_branches

N = 7
CEILING_SECONDS = 10  # ~7x the ~1.4s observed post-fix; 43s pre-fix would fail this easily


class TestSearchScalesWithTieCountNotQuadratically(unittest.TestCase):

    def test_seven_concurrent_branches_completes_well_under_ceiling_with_correct_tie_count(self):
        tree = _tree_with_n_concurrent_optional_branches(N)
        net, im, fm, activity_to_id, tau_ids, id_loop_list = EbiOccurance().build_petri_net(tree)

        start = time.time()
        result = align_pn_all([activity_to_id['act0']], net, im, fm, id_loop_list, timeout=60, tau_ids=tau_ids)
        elapsed = time.time() - start

        _, (opt_agns, timed_out, _) = result[0]
        self.assertEqual(timed_out, 0)
        self.assertEqual(len(opt_agns), math.factorial(N))
        self.assertLess(elapsed, CEILING_SECONDS)


if __name__ == '__main__':
    unittest.main()
