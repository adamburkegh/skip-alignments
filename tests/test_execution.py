"""
Unit tests for skip_alignments.execution.ExecutionTree's interleaving logic.

These pin down the linearization mechanism that underlies
ExecutionManager.coninciding_agns: coinciding alignments differ only by
reordering concurrent (And) moves, never by padding in extra moves. See the
interleaving operator diamond and its worked example in the paper's
preliminaries.

Also covers the two mitigations added for the And-node combinatorial blowup
found via process-voids' BPI2013 Incidents dose-response run (a single
lengths=[15,11] interleaving call recurred for 24+ minutes with no
completion, and a lengths=[13,9] call generated 41.78 million candidates to
keep 659):
  - TestSyncOnlyMerge / TestConstrainedAndShuffle: when every element being
    interleaved is a synchronous move, its final order is already uniquely
    determined by the log's own order (see remove_log_moves' position
    suffixing), so there is exactly one valid interleaving -- no
    combinatorial search is needed at all. _sync_only_merge computes that
    directly; ExecutionTree.shuffle's And-branch uses it whenever every
    branch is purely synchronous, falling back to the unconstrained
    all_order_preserving_shuffles (unchanged) whenever a model-only move is
    present.
  - TestShuffleExplosionCeiling: MAX_SHUFFLE_COUNT is a hard, catchable stop
    for whatever combinatorial search remains, so a pathological trace can't
    silently consume an entire run's time budget with no feedback.

Run with:
    python -m unittest tests.test_execution -v
"""
import contextlib
import io
import math
import unittest

from skip_alignments.execution import (
    ExecutionManager, ExecutionTree, ShuffleExplosionError, _predicted_shuffle_count, _sync_only_merge,
    get_large_shuffle_count_threshold, get_max_shuffle_count, set_large_shuffle_count_threshold,
    set_max_shuffle_count,
)
from skip_alignments.processtree import Activity, And, Execution, Sequence


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


class TestPredictedShuffleCount(unittest.TestCase):
    """_predicted_shuffle_count must equal the real output size of
    all_order_preserving_shuffles, computed analytically (no enumeration) --
    both the ceiling (Part A) and the warning depend on this being exact."""

    def test_matches_real_output_size_two_paths(self):
        for l1, l2 in [(1, 1), (2, 3), (3, 3), (0, 4)]:
            paths = [list(range(l1)), list(range(100, 100 + l2))]
            with self.subTest(l1=l1, l2=l2):
                real = len(_et().all_order_preserving_shuffles(*paths))
                self.assertEqual(_predicted_shuffle_count(paths), real)

    def test_matches_real_output_size_three_paths(self):
        paths = [['a', 'b'], ['c'], ['d', 'e']]
        real = len(_et().all_order_preserving_shuffles(*paths))
        self.assertEqual(_predicted_shuffle_count(paths), real)

    def test_matches_multinomial_formula(self):
        # multinomial(9; 3,4,2) = 9! / (3! 4! 2!) = 1260
        paths = [list(range(3)), list(range(4)), list(range(2))]
        self.assertEqual(_predicted_shuffle_count(paths), 1260)


class TestSyncOnlyMerge(unittest.TestCase):
    """
    _sync_only_merge(paths, sync_rank): when every (log, model) element
    across all paths has a log-side entry present in sync_rank (i.e. is a
    synchronous move, not a model-only '>>' move), the merged order is fully
    determined by sync_rank -- returns that single merge. Returns None the
    moment any element's log side is absent from sync_rank (a model-only
    move), signalling the caller to fall back to the general algorithm.
    """

    def test_all_sync_two_paths_returns_the_single_correct_merge(self):
        path1 = [('a0', 'a0'), ('b1', 'b1')]
        path2 = [('c2', 'c2'), ('d3', 'd3')]
        sync_rank = {'a0': 0, 'b1': 1, 'c2': 2, 'd3': 3}
        result = _sync_only_merge([path1, path2], sync_rank)
        self.assertEqual(result, [[('a0', 'a0'), ('b1', 'b1'), ('c2', 'c2'), ('d3', 'd3')]])

    def test_all_sync_result_is_order_preserving_per_branch(self):
        # ranks interleaved across branches, not just concatenated
        path1 = [('a0', 'a0'), ('c2', 'c2')]
        path2 = [('b1', 'b1'), ('d3', 'd3')]
        sync_rank = {'a0': 0, 'b1': 1, 'c2': 2, 'd3': 3}
        result = _sync_only_merge([path1, path2], sync_rank)
        self.assertEqual(result, [[('a0', 'a0'), ('b1', 'b1'), ('c2', 'c2'), ('d3', 'd3')]])

    def test_result_is_a_member_of_the_unconstrained_output(self):
        # cross-check against the trusted, already-tested primitive: the
        # fast-path answer must be one of the outputs the general algorithm
        # would also produce (a genuine subset, not something impossible)
        path1 = [('a0', 'a0'), ('b1', 'b1'), ('e4', 'e4')]
        path2 = [('c2', 'c2'), ('d3', 'd3')]
        sync_rank = {'a0': 0, 'b1': 1, 'c2': 2, 'd3': 3, 'e4': 4}
        fast = _sync_only_merge([path1, path2], sync_rank)
        general = _et().all_order_preserving_shuffles(path1, path2)
        self.assertEqual(len(fast), 1)
        self.assertIn(fast[0], general)

    def test_model_only_move_falls_back_to_none(self):
        path1 = [('a0', 'a0'), ('>>', 'b_model')]
        path2 = [('c2', 'c2'), ('d3', 'd3')]
        sync_rank = {'a0': 0, 'c2': 1, 'd3': 2}
        self.assertIsNone(_sync_only_merge([path1, path2], sync_rank))

    def test_empty_paths_return_empty_merge(self):
        result = _sync_only_merge([[], []], {})
        self.assertEqual(result, [[]])

    def test_three_paths_all_sync(self):
        path1 = [('a0', 'a0')]
        path2 = [('b1', 'b1')]
        path3 = [('c2', 'c2')]
        sync_rank = {'a0': 0, 'b1': 1, 'c2': 2}
        result = _sync_only_merge([path1, path2, path3], sync_rank)
        self.assertEqual(result, [[('a0', 'a0'), ('b1', 'b1'), ('c2', 'c2')]])


def _leaf(node, start, stop):
    return ExecutionTree(None, [], Execution(node, start, stop))


def _seq_tree(seq_node, children, start, stop):
    tree = ExecutionTree(None, list(children), Execution(seq_node, start, stop))
    for c in tree.children:
        c.set_parent(tree)
    return tree


class _FakeState:
    """ExecutionTree.shuffle only ever reads state.path -- no need for a
    real alignment.State (which requires a Mapper) to test it directly."""
    def __init__(self, path):
        self.path = path


class TestConstrainedAndShuffle(unittest.TestCase):
    """
    Integration-level check of ExecutionTree.shuffle's And-branch: built by
    hand (bypassing the full alignment search) so the And case can be
    exercised directly, matching how ExecutionManager.build_execution_tree
    would have wired it up for a real state.
    """

    def _and_of_two_pairs(self):
        # And(Seq(a,b), Seq(c,d)), matching path positions 0..3
        a = Activity(None, 'a', 100)
        b = Activity(None, 'b', 100)
        c = Activity(None, 'c', 100)
        d = Activity(None, 'd', 100)
        seq1_node = Sequence(None, [a, b])
        seq2_node = Sequence(None, [c, d])
        and_node = And(None, [seq1_node, seq2_node])

        leaf_a, leaf_b = _leaf(a, 0, 1), _leaf(b, 1, 2)
        leaf_c, leaf_d = _leaf(c, 2, 3), _leaf(d, 3, 4)
        seq1 = _seq_tree(seq1_node, [leaf_a, leaf_b], 0, 2)
        seq2 = _seq_tree(seq2_node, [leaf_c, leaf_d], 2, 4)
        and_tree = ExecutionTree(None, [seq1, seq2], Execution(and_node, 0, 4))
        seq1.set_parent(and_tree)
        seq2.set_parent(and_tree)
        return and_tree

    def test_all_sync_uses_fast_path_and_matches_unconstrained_result(self):
        and_tree = self._and_of_two_pairs()
        path = [('a0', 'a0'), ('b1', 'b1'), ('c2', 'c2'), ('d3', 'd3')]
        state = _FakeState(path)
        sync_rank = {'a0': 0, 'b1': 1, 'c2': 2, 'd3': 3}

        fast_result = and_tree.shuffle(state, sync_rank=sync_rank)
        unconstrained_result = and_tree.shuffle(state, sync_rank=None)

        # exactly one candidate from the fast path, and it's a member of
        # what the unconstrained (6-candidate) algorithm would produce
        self.assertEqual(len(fast_result), 1)
        self.assertIn(fast_result[0], unconstrained_result)
        self.assertEqual(len(unconstrained_result), 6)
        # and it must be the one whose sync order matches sync_rank exactly
        self.assertEqual(fast_result[0], path)

    def test_fast_path_result_equals_the_unconstrained_result_after_the_real_filter(self):
        # stronger than membership: replicate ExecutionManager.shuffle's own
        # filter criterion (log side, '>>' entries excluded, must equal the
        # true sync order) against the UNCONSTRAINED candidate set, and
        # assert that filtered subset is exactly the fast path's answer --
        # not just that the fast answer happens to be one of six candidates
        and_tree = self._and_of_two_pairs()
        path = [('a0', 'a0'), ('b1', 'b1'), ('c2', 'c2'), ('d3', 'd3')]
        state = _FakeState(path)
        sync_rank = {'a0': 0, 'b1': 1, 'c2': 2, 'd3': 3}
        sync_moves_log = ['a0', 'b1', 'c2', 'd3']

        fast_result = and_tree.shuffle(state, sync_rank=sync_rank)
        unconstrained_result = and_tree.shuffle(state, sync_rank=None)

        filtered_unconstrained = {
            tuple(agn) for agn in unconstrained_result
            if [l for l, _ in agn if l != '>>'] == sync_moves_log
        }
        self.assertEqual(filtered_unconstrained, {tuple(fast_result[0])})

    def test_model_only_move_falls_back_and_still_contains_correct_answer(self):
        and_tree = self._and_of_two_pairs()
        # b is a model-only move now
        path = [('a0', 'a0'), ('>>', 'b1'), ('c2', 'c2'), ('d3', 'd3')]
        state = _FakeState(path)
        sync_rank = {'a0': 0, 'c2': 1, 'd3': 2}

        result = and_tree.shuffle(state, sync_rank=sync_rank)
        unconstrained_result = and_tree.shuffle(state, sync_rank=None)

        # falls back to the general algorithm -- same result set either way
        self.assertEqual(result, unconstrained_result)
        self.assertEqual(len(result), 6)
        # the original path itself must still be one of the candidates
        self.assertIn(path, result)

    def test_default_sync_rank_none_preserves_old_unconstrained_behaviour(self):
        and_tree = self._and_of_two_pairs()
        path = [('a0', 'a0'), ('b1', 'b1'), ('c2', 'c2'), ('d3', 'd3')]
        state = _FakeState(path)
        # calling without sync_rank at all (positional/default) must behave
        # exactly as before this change -- no silent behaviour change for
        # any caller that doesn't opt in
        self.assertEqual(and_tree.shuffle(state), and_tree.shuffle(state, sync_rank=None))
        self.assertEqual(len(and_tree.shuffle(state)), 6)


class TestExecutionManagerShuffleWithFastPath(unittest.TestCase):
    """
    Closes the gap the above tests leave: ExecutionManager.shuffle is the
    real call site (used by coninciding_agns) that builds sync_rank from
    sync_moves_log and does log-move merging -- everything above calls
    ExecutionTree.shuffle directly with a hand-built sync_rank. This
    exercises the actual production path end to end, and cross-checks it
    against the same computation with the fast path forced off (via
    ExecutionTree.shuffle(state, sync_rank=None) + the identical filter
    ExecutionManager.shuffle itself applies), so the "on" and "off" paths
    are compared through the real filtering/merging logic, not just at the
    ExecutionTree.shuffle level.
    """

    def _and_tree_and_state(self):
        a = Activity(None, 'a', 100)
        b = Activity(None, 'b', 100)
        c = Activity(None, 'c', 100)
        d = Activity(None, 'd', 100)
        seq1_node = Sequence(None, [a, b])
        seq2_node = Sequence(None, [c, d])
        and_node = And(None, [seq1_node, seq2_node])

        leaf_a, leaf_b = _leaf(a, 0, 1), _leaf(b, 1, 2)
        leaf_c, leaf_d = _leaf(c, 2, 3), _leaf(d, 3, 4)
        seq1 = _seq_tree(seq1_node, [leaf_a, leaf_b], 0, 2)
        seq2 = _seq_tree(seq2_node, [leaf_c, leaf_d], 2, 4)
        and_tree = ExecutionTree(None, [seq1, seq2], Execution(and_node, 0, 4))
        seq1.set_parent(and_tree)
        seq2.set_parent(and_tree)

        path = [('a0', 'a0'), ('b1', 'b1'), ('c2', 'c2'), ('d3', 'd3')]
        return and_tree, _FakeState(path)

    def test_real_call_site_yields_the_single_correct_result(self):
        and_tree, state = self._and_tree_and_state()
        log_moves = []
        log_path = ['a0', 'b1', 'c2', 'd3']  # no log-only moves: log_path == sync_moves_log

        yielded = list(ExecutionManager().shuffle(state, log_moves, log_path, and_tree))

        self.assertEqual(yielded, [[('a0', 'a0'), ('b1', 'b1'), ('c2', 'c2'), ('d3', 'd3')]])

    def test_real_call_site_matches_fast_path_forced_off(self):
        and_tree, state = self._and_tree_and_state()
        log_moves = []
        log_path = ['a0', 'b1', 'c2', 'd3']
        sync_moves_log = [l for l in log_path if l not in log_moves]

        # the real path: ExecutionManager.shuffle builds sync_rank itself
        # and uses it (fast path engaged)
        on_result = list(ExecutionManager().shuffle(state, log_moves, log_path, and_tree))

        # replicate the same filtering/merging ExecutionManager.shuffle
        # does, but against ExecutionTree.shuffle(state, sync_rank=None) --
        # i.e. the fast path forced off, going through the unconstrained
        # all_order_preserving_shuffles for every candidate
        off_candidates = and_tree.shuffle(state, sync_rank=None)
        off_result = [agn for agn in off_candidates
                      if [l for l, _ in agn if l != '>>'] == sync_moves_log
                      and [l for l, _ in agn if l != '>>'] == log_path]

        self.assertEqual(on_result, off_result)
        self.assertEqual(len(on_result), 1)


def _and_of_four_leaves():
    # And(a, b, c, d): four leaves directly under one And node, each a
    # single-move branch -- 4! = 24 raw interleavings, of which exactly one
    # is order-preserving-consistent with any given total order over the
    # four moves. Wider than the 2x2 case used elsewhere in this file, per
    # a review of the fast-path/ceiling design asking specifically for
    # coinciding-set equality on a 3-4 leaf And, not just a 2-branch one.
    nodes = [Activity(None, name, 100) for name in ('a', 'b', 'c', 'd')]
    and_node = And(None, nodes)
    leaves = [_leaf(n, i, i + 1) for i, n in enumerate(nodes)]
    and_tree = ExecutionTree(None, leaves, Execution(and_node, 0, 4))
    for lf in leaves:
        lf.set_parent(and_tree)
    path = [('a0', 'a0'), ('b1', 'b1'), ('c2', 'c2'), ('d3', 'd3')]
    return and_tree, path


class TestFourLeafAndNodeCoincidingSetEquality(unittest.TestCase):
    """
    Directly answers the concern raised reviewing the fast path: not just
    that _sync_only_merge's answer is *a* valid candidate, but that the
    coinciding *set* -- the thing |C[state]| in coninciding_agns actually
    counts, which _skip_agn_probs sums a term over once per member -- is
    identical with the fast path on vs. off, on an And wider than the 2x2
    case used elsewhere in this file (four leaves, not two branches of two).
    """

    def test_fast_path_set_equals_filtered_unconstrained_set(self):
        and_tree, path = _and_of_four_leaves()
        state = _FakeState(path)
        sync_rank = {log: i for i, (log, _) in enumerate(path)}
        sync_moves_log = [log for log, _ in path]

        fast_result = and_tree.shuffle(state, sync_rank=sync_rank)
        unconstrained_result = and_tree.shuffle(state, sync_rank=None)

        self.assertEqual(len(unconstrained_result), 24)  # 4! raw interleavings
        filtered_unconstrained = {
            tuple(agn) for agn in unconstrained_result
            if [l for l, _ in agn if l != '>>'] == sync_moves_log
        }
        self.assertEqual(len(fast_result), 1)
        self.assertEqual(filtered_unconstrained, {tuple(fast_result[0])})

    def test_real_call_site_coinciding_set_matches_fast_path_forced_off(self):
        and_tree, path = _and_of_four_leaves()
        state = _FakeState(path)
        log_moves = []
        log_path = [log for log, _ in path]
        sync_moves_log = [l for l in log_path if l not in log_moves]

        on_result = list(ExecutionManager().shuffle(state, log_moves, log_path, and_tree))

        off_candidates = and_tree.shuffle(state, sync_rank=None)
        off_result = [agn for agn in off_candidates
                      if [l for l, _ in agn if l != '>>'] == sync_moves_log
                      and [l for l, _ in agn if l != '>>'] == log_path]

        self.assertEqual(on_result, off_result)
        self.assertEqual(len(on_result), 1)


def _two_leaf_and_tree(big_n):
    # And(a-run-of-length-big_n, b-run-of-length-big_n), each leaf a run of
    # big_n distinctly-labelled moves -- used to force a specific,
    # controlled multinomial coefficient (comb(2*big_n, big_n)) without
    # needing branches as long as the real BPI2013 case.
    a_node = Activity(None, 'a', 100)
    b_node = Activity(None, 'b', 100)
    and_node = And(None, [a_node, b_node])
    leaf_a = ExecutionTree(None, [], Execution(a_node, 0, big_n))
    leaf_b = ExecutionTree(None, [], Execution(b_node, big_n, 2 * big_n))
    and_tree = ExecutionTree(None, [leaf_a, leaf_b], Execution(and_node, 0, 2 * big_n))
    leaf_a.set_parent(and_tree)
    leaf_b.set_parent(and_tree)
    path_a = [(f'a{i}', f'a{i}') for i in range(big_n)]
    path_b = [(f'b{i}', f'b{i}') for i in range(big_n)]
    return and_tree, path_a, path_b


def _smallest_n_exceeding(count):
    n = 1
    while math.comb(2 * n, n) <= count:
        n += 1
    return n


class TestShuffleExplosionCeiling(unittest.TestCase):
    """
    MAX_SHUFFLE_COUNT: a hard, catchable stop for an And-branch interleaving
    whose predicted output size is too large to ever reasonably finish --
    added after a BPI2013 Incidents run stalled 24+ minutes with no
    completion and no way to bound or catch it (see execution.py).
    """

    def test_predicted_count_within_ceiling_does_not_raise(self):
        paths = [list(range(3)), list(range(2))]
        self.assertLessEqual(_predicted_shuffle_count(paths), get_max_shuffle_count())
        # should complete normally
        result = _et().all_order_preserving_shuffles(*paths)
        self.assertEqual(len(result), math.comb(5, 2))

    def test_and_branch_raises_shuffle_explosion_error_above_ceiling(self):
        big_n = _smallest_n_exceeding(get_max_shuffle_count())
        and_tree, path_a, path_b = _two_leaf_and_tree(big_n)
        state = _FakeState(path_a + path_b)
        # no sync_rank given (or a sync_rank that still leaves ambiguity is
        # irrelevant here -- with none given the fast path never applies)
        with self.assertRaises(ShuffleExplosionError) as ctx:
            and_tree.shuffle(state)
        self.assertEqual(ctx.exception.predicted_count, math.comb(2 * big_n, big_n))

    def test_fast_path_bypasses_the_ceiling_entirely(self):
        # the whole point of the sync-only fast path: a purely-synchronous
        # And-branch never even computes a combinatorial candidate count,
        # so it can't trip the ceiling no matter how long the branches are
        big_n = _smallest_n_exceeding(get_max_shuffle_count())
        and_tree, path_a, path_b = _two_leaf_and_tree(big_n)
        combined = path_a + path_b
        sync_rank = {log: i for i, (log, _) in enumerate(combined)}
        state = _FakeState(combined)

        result = and_tree.shuffle(state, sync_rank=sync_rank)
        self.assertEqual(result, [combined])


class TestConfigurableThresholds(unittest.TestCase):
    """
    MAX_SHUFFLE_COUNT/LARGE_SHUFFLE_COUNT_THRESHOLD are process-wide config,
    not per-call context (unlike sync_rank) -- get_/set_ functions mutate
    the same module globals shuffle()'s And-branch already reads directly,
    so a set_*() call takes effect for every subsequent call at any
    recursion depth, with nothing threaded through any call signature.
    """

    def setUp(self):
        self._orig_max = get_max_shuffle_count()
        self._orig_warn = get_large_shuffle_count_threshold()

    def tearDown(self):
        # module-global state -- must not leak between tests
        set_max_shuffle_count(self._orig_max)
        set_large_shuffle_count_threshold(self._orig_warn)

    def test_get_set_round_trip(self):
        set_max_shuffle_count(42)
        self.assertEqual(get_max_shuffle_count(), 42)
        set_large_shuffle_count_threshold(7)
        self.assertEqual(get_large_shuffle_count_threshold(), 7)

    def test_set_max_shuffle_count_rejects_non_positive(self):
        with self.assertRaises(ValueError):
            set_max_shuffle_count(0)
        with self.assertRaises(ValueError):
            set_max_shuffle_count(-5)

    def test_set_large_shuffle_count_threshold_rejects_non_positive(self):
        with self.assertRaises(ValueError):
            set_large_shuffle_count_threshold(0)

    def test_lowering_max_shuffle_count_changes_and_branch_behaviour(self):
        # comb(6,3) = 20 -- ordinarily well within the default ceiling, but
        # not once it's lowered below 20, and no threading was needed to
        # make shuffle() see the new value
        set_max_shuffle_count(10)
        and_tree, path_a, path_b = _two_leaf_and_tree(3)
        state = _FakeState(path_a + path_b)
        with self.assertRaises(ShuffleExplosionError) as ctx:
            and_tree.shuffle(state)
        self.assertEqual(ctx.exception.predicted_count, 20)

    def test_raised_error_message_reflects_current_ceiling(self):
        set_max_shuffle_count(10)
        and_tree, path_a, path_b = _two_leaf_and_tree(3)
        state = _FakeState(path_a + path_b)
        with self.assertRaises(ShuffleExplosionError) as ctx:
            and_tree.shuffle(state)
        self.assertIn("MAX_SHUFFLE_COUNT=10", str(ctx.exception))


class TestConincidingAgnsEmptyInputs(unittest.TestCase):
    """
    ExecutionManager.coninciding_agns divides by len(ratio_per_var) (all
    variants) and, per variant, by len(states) -- both zero-able and
    previously unguarded. Reported from a real process-voids run: an empty
    skip_dict (no variants at all, e.g. a degraded-to-nothing log) crashed
    with ZeroDivisionError at the final compression-ratio print.
    """

    def test_empty_skip_dict_does_not_raise(self):
        with contextlib.redirect_stdout(io.StringIO()):
            C, global_C = ExecutionManager().coninciding_agns({})
        self.assertEqual(C, {})
        self.assertEqual(global_C, {})

    def test_variant_with_zero_states_does_not_raise(self):
        with contextlib.redirect_stdout(io.StringIO()):
            C, global_C = ExecutionManager().coninciding_agns({'some_var': []})
        self.assertEqual(C, {})
        self.assertEqual(global_C, {'some_var': []})

if __name__ == '__main__':
    unittest.main()
