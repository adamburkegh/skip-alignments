"""
Unit tests for alignall's tied-optimal-alignment ceiling (LARGE_TIE_COUNT_THRESHOLD/
MAX_TIE_COUNT/TieExplosionError), the classical-alignment-path counterpart
to execution.py's LARGE_SHUFFLE_COUNT_THRESHOLD/MAX_SHUFFLE_COUNT/
ShuffleExplosionError (see test_execution.py's TestConfigurableThresholds
for the same pattern on that path).

Root cause this guards against: a classical Petri net has no notion of
"these two moves are concurrent, their relative order doesn't matter" the
way the process-tree-native path's And node does -- state_equation_a_star's
search counts every interleaving of concurrent branches as a genuinely
distinct optimal alignment. Found via a process-voids report on rtfm.xes:
a 2-activity variant against a heavily And/Xor(Tau,_)-nested discovered
tree produced 4115+ tied optimal alignments well before the per-trace
timeout was reached (confirmed directly: id_loop_list was empty for that
tree -- no loops at all -- and is_cycling never fired once in 33,116
calls, ruling out the loop cycle-guard as either cause or fix here).

Run with:
    python -m unittest tests.test_alignall_tie_explosion -v
"""
import pickle
import unittest

from skipalignments.alignall import (
    TieExplosionError, align_pn_all, get_large_tie_count_threshold, get_max_tie_count,
    set_large_tie_count_threshold, set_max_tie_count,
)
from skipalignments.probabilities import EbiOccurance
from skipalignments.processtree import Activity, And, Tau, Xor


def _tree_with_n_concurrent_optional_branches(n):
    """
    And(Xor(Tau, act1), Xor(Tau, act2), ..., Xor(Tau, actN)): n
    independent, concurrent, individually-optional branches. Against an
    empty trace, the cheapest (cost 0) resolution is every branch taking
    its Tau alternative -- but those n concurrent zero-cost firings have
    n! distinct interleavings, each counted as its own tied optimal
    alignment by the classical aligner (no concurrency-equivalence
    collapsing, unlike the process-tree-native And node's shuffle()).
    """
    branches = []
    for i in range(n):
        act = Activity(None, f'act{i}', 100000)
        act.id = f"a{i}"
        tau = Tau(None, f'tau{i}', 0)
        tau.id = f"t{i}"
        choice = Xor(None, [act, tau])
        act.set_parent(choice)
        tau.set_parent(choice)
        choice.id = f"c{i}"
        branches.append(choice)

    tree = And(None, branches)
    for b in branches:
        b.set_parent(tree)
    tree.id = "root"
    return tree


class TestConfigurableTieThresholds(unittest.TestCase):

    def setUp(self):
        self._orig_max = get_max_tie_count()
        self._orig_warn = get_large_tie_count_threshold()

    def tearDown(self):
        # module-global state -- must not leak between tests
        set_max_tie_count(self._orig_max)
        set_large_tie_count_threshold(self._orig_warn)

    def test_get_set_round_trip(self):
        set_max_tie_count(42)
        self.assertEqual(get_max_tie_count(), 42)
        set_large_tie_count_threshold(7)
        self.assertEqual(get_large_tie_count_threshold(), 7)

    def test_set_max_tie_count_rejects_non_positive(self):
        with self.assertRaises(ValueError):
            set_max_tie_count(0)
        with self.assertRaises(ValueError):
            set_max_tie_count(-5)

    def test_set_large_tie_count_threshold_rejects_non_positive(self):
        with self.assertRaises(ValueError):
            set_large_tie_count_threshold(0)

    def test_lowering_max_tie_count_raises_on_a_real_concurrent_tie_explosion(self):
        # 4 concurrent optional branches, one resolved via a real sync
        # move (act0), the rest via their tau alternative -- still 4
        # concurrent, unordered branch firings -> 4! = 24 tied optimal
        # alignments, well above a ceiling of 5. (An empty trace hits an
        # unrelated pm4py variant-extraction edge case, so a one-event
        # trace is used instead -- see test_default_ceiling_... below.)
        set_max_tie_count(5)
        tree = _tree_with_n_concurrent_optional_branches(4)
        net, im, fm, activity_to_id, tau_ids, id_loop_list = EbiOccurance().build_petri_net(tree)

        with self.assertRaises(TieExplosionError) as ctx:
            align_pn_all([activity_to_id['act0']], net, im, fm, id_loop_list, timeout=30, tau_ids=tau_ids)
        self.assertIn("MAX_TIE_COUNT=5", str(ctx.exception))
        self.assertGreater(ctx.exception.tie_count, 5)

    def test_tie_explosion_error_survives_a_pickle_round_trip(self):
        # the __reduce__ override this pins directly: Exception's default
        # pickling reconstructs via type(self)(*self.args), and self.args
        # is just (the formatted message,) after __init__'s
        # super().__init__(message) call -- not the (partial_agns,
        # tie_count, ceiling) __init__ actually needs. Without __reduce__,
        # this raises a confusing unrelated TypeError instead of
        # reconstructing the exception -- exactly what happens when
        # align_pn_all_multi's workers pickle it back across a
        # ProcessPoolExecutor boundary (see test_alignall_multiprocessing.py,
        # a slower end-to-end confirmation of the same fix).
        original = TieExplosionError(partial_agns=['a', 'b'], tie_count=7, ceiling=5)
        restored = pickle.loads(pickle.dumps(original))

        self.assertEqual(restored.partial_agns, ['a', 'b'])
        self.assertEqual(restored.tie_count, 7)
        self.assertEqual(restored.ceiling, 5)
        self.assertEqual(str(restored), str(original))

    def test_explicit_max_tie_count_override_beats_the_global(self):
        # per-call max_tie_count doesn't require touching the global at
        # all -- confirms the parameter path independent of
        # set_max_tie_count, which the tests above already cover
        tree = _tree_with_n_concurrent_optional_branches(4)  # 4! = 24 ties
        net, im, fm, activity_to_id, tau_ids, id_loop_list = EbiOccurance().build_petri_net(tree)
        unchanged_global = get_max_tie_count()

        with self.assertRaises(TieExplosionError) as ctx:
            align_pn_all([activity_to_id['act0']], net, im, fm, id_loop_list, timeout=30, tau_ids=tau_ids,
                         max_tie_count=5)
        self.assertIn("MAX_TIE_COUNT=5", str(ctx.exception))
        self.assertEqual(get_max_tie_count(), unchanged_global, "the global must be untouched by a per-call override")

    def test_default_ceiling_does_not_interfere_with_a_small_case(self):
        # 2 concurrent optional branches -> only 2! = 2 ties, nowhere near
        # the default ceiling -- must complete normally, not raise
        tree = _tree_with_n_concurrent_optional_branches(2)
        net, im, fm, activity_to_id, tau_ids, id_loop_list = EbiOccurance().build_petri_net(tree)

        result = align_pn_all([activity_to_id['act0']], net, im, fm, id_loop_list, timeout=30, tau_ids=tau_ids)
        _, (opt_agns, timed_out, _) = result[0]
        self.assertEqual(timed_out, 0)
        self.assertEqual(len(opt_agns), 2)

    def test_large_tie_count_threshold_logs_a_warning(self):
        # threshold well below max, so the search completes normally
        # (24 ties from 4 branches) but crosses the warning threshold on
        # the way -- the warning fires exactly once, when the count
        # equals the threshold precisely (see __search's `elif ==` check)
        set_max_tie_count(1000)
        set_large_tie_count_threshold(3)
        tree = _tree_with_n_concurrent_optional_branches(4)
        net, im, fm, activity_to_id, tau_ids, id_loop_list = EbiOccurance().build_petri_net(tree)

        with self.assertLogs('skipalignments.alignall', level='WARNING') as ctx:
            result = align_pn_all([activity_to_id['act0']], net, im, fm, id_loop_list, timeout=30, tau_ids=tau_ids)

        _, (opt_agns, timed_out, _) = result[0]
        self.assertEqual(timed_out, 0)
        self.assertEqual(len(opt_agns), 24)
        self.assertTrue(any("accumulated 3 tied optimal alignments" in msg for msg in ctx.output))

    def test_below_large_tie_count_threshold_logs_nothing(self):
        # 2 ties, well under a threshold of 3 -- no warning should fire
        set_max_tie_count(1000)
        set_large_tie_count_threshold(3)
        tree = _tree_with_n_concurrent_optional_branches(2)
        net, im, fm, activity_to_id, tau_ids, id_loop_list = EbiOccurance().build_petri_net(tree)

        with self.assertNoLogs('skipalignments.alignall', level='WARNING'):
            result = align_pn_all([activity_to_id['act0']], net, im, fm, id_loop_list, timeout=30, tau_ids=tau_ids)

        _, (opt_agns, timed_out, _) = result[0]
        self.assertEqual(timed_out, 0)
        self.assertEqual(len(opt_agns), 2)


if __name__ == '__main__':
    unittest.main()
