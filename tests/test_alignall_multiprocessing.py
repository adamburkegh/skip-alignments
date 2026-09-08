"""
Unit tests for align_pn_all_multi (the multiprocessing path), which had
zero coverage for any of this session's fixes: tau_ids/id_loop_list
threading, and TieExplosionError specifically.

The suspected risk going in was pickling: TieExplosionError carries pm4py
Transition objects in .partial_agns, raised inside a worker process and
sent back across the ProcessPoolExecutor boundary to the caller. What was
actually broken was more fundamental, and confirmed directly (see
CHANGELOG.md): MAX_TIE_COUNT/LARGE_TIE_COUNT_THRESHOLD are plain module
globals, mutated in-process by set_max_tie_count/
set_large_tie_count_threshold. On Windows (the default 'spawn' start
method), each ProcessPoolExecutor worker is a *fresh* interpreter that
re-imports alignall from scratch -- it never saw the parent's mutated
value, only the module's default, so set_max_tie_count(1) in the caller
had *zero effect* on align_pn_all_multi's actual worker-side enforcement.
Fixed by resolving get_max_tie_count()/get_large_tie_count_threshold()
once in the parent (apply_multiprocessing) and passing the concrete
values through explicitly to every worker's apply_trace call, rather than
relying on each worker's own copy of the global.

Run with:
    python -m unittest tests.test_alignall_multiprocessing -v
"""
import unittest

import pandas as pd

from skipalignments.alignall import TieExplosionError, align_pn_all_multi, set_max_tie_count
from skipalignments.probabilities import EbiOccurance

from tests.test_alignall_tie_explosion import _tree_with_n_concurrent_optional_branches


def _one_case_log(activity_id):
    return pd.DataFrame({
        'case:concept:name': ['1'],
        'concept:name': [activity_id],
        'time:timestamp': [pd.Timestamp(year=2000, month=1, day=1)],
    })


class TestAlignPnAllMultiThreading(unittest.TestCase):

    def setUp(self):
        self._orig_max = None  # set per-test only where needed

    def test_threads_tau_ids_and_id_loop_list_to_the_correct_cost_and_tie_count(self):
        tree = _tree_with_n_concurrent_optional_branches(2)
        net, im, fm, activity_to_id, tau_ids, id_loop_list = EbiOccurance().build_petri_net(tree)
        log = _one_case_log(activity_to_id['act0'])

        result = align_pn_all_multi(log, net, im, fm, id_loop_list, tree=tree, timeout=30, tau_ids=tau_ids)

        self.assertEqual(len(result), 1)
        (time_end, (opt_agns, timed_out, _)), = result.values()
        self.assertEqual(timed_out, 0)
        self.assertEqual(len(opt_agns), 2)
        # cost 0 confirms tau_ids reached the worker process correctly --
        # without it, act1's Tau alternative would be mispriced at 100000
        # (see test_alignall_tau_cost.py)
        self.assertEqual(opt_agns[0]['cost'], 0)


class TestAlignPnAllMultiTieExplosionPropagation(unittest.TestCase):

    def setUp(self):
        self._orig_max = 100_000
        from skipalignments.alignall import get_max_tie_count
        self._orig_max = get_max_tie_count()

    def tearDown(self):
        set_max_tie_count(self._orig_max)

    def test_tie_explosion_error_survives_the_process_boundary(self):
        set_max_tie_count(1)
        tree = _tree_with_n_concurrent_optional_branches(3)  # 3! = 6 > 1
        net, im, fm, activity_to_id, tau_ids, id_loop_list = EbiOccurance().build_petri_net(tree)
        log = _one_case_log(activity_to_id['act0'])

        with self.assertRaises(TieExplosionError) as ctx:
            align_pn_all_multi(log, net, im, fm, id_loop_list, tree=tree, timeout=30, tau_ids=tau_ids)
        self.assertIn("MAX_TIE_COUNT=1", str(ctx.exception))


if __name__ == '__main__':
    unittest.main()
