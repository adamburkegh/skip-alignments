"""
Unit tests for skipalignments.skips.Skipper's per-node skip/non-skip
counting.

Run with:
    python -m unittest tests.test_skips -v
"""
import os
import shutil
import tempfile
import unittest

import pandas as pd

from skipalignments import DerivationPipeline
from skipalignments.alignment import Aligner
from skipalignments.processtree import Activity, Loop, Sequence
from skipalignments.skips import Skipper

from tests.test_run_example import build_example_tree

COST = 100000


class TestSkipperLoopDoChildCounting(unittest.TestCase):
    """
    Skipper.count_skip_executions/count_non_skip_executions check
    `m.node == tree` by exact identity (skips.py), with no awareness that a
    leaf can be masked inside a coarser Skip/TauPath wrapper placed on one
    of its ancestors. This is fine when the aligner localizes a Skip
    directly onto the leaf (as it does for an isolated gap flanked by
    present neighbours), but when a whole compound subtree goes
    unwitnessed as a block, the Skip wraps the ancestor once instead -- and
    the leaf's own skip/non-skip counts silently read as 0/0.

    Confirmed structurally (see probe_missing_subtree.py in this session):
    the same masking happens for Loop's do-child and for a nested Sequence
    alike, so this isn't Loop-specific -- it's any leaf under a wholly
    unwitnessed multi-node compound subtree.
    """

    @classmethod
    def setUpClass(cls):
        cls.tree = build_example_tree()
        seq1, loop = cls.tree.children
        cls.a, cls.choice = seq1.children
        cls.c, cls.d = cls.choice.children
        cls.do_child, cls.redo_child = loop.children  # 'e', 'f'

        Aligner.set_level_incentive(0)
        aligner = Aligner(cls.tree)
        # trace has 'a' and 'c' but no occurrence of the loop's activities at all
        states, _ = aligner.align2(['a', 'c'], [COST, COST], all_optimal=True, timeout=30)
        cls.states = states

    def test_loop_do_child_counted_as_skipped_when_whole_loop_unwitnessed(self):
        skipper = Skipper()
        for state in self.states:
            prob = skipper._conditional_skip_prob(self.do_child, state)
            # the loop's do-child never occurs in the trace at all, and the
            # whole loop is masked inside one lump Skip on the Loop node
            # itself -- _conditional_skip_prob should recognise the
            # do-child as contained within that Skip and treat it as fully
            # skipped, rather than as "never reached".
            self.assertEqual(prob, 1)

    def test_loop_redo_child_also_counted_as_skipped(self):
        # same reasoning as the do-child, for the loop's other child
        skipper = Skipper()
        for state in self.states:
            prob = skipper._conditional_skip_prob(self.redo_child, state)
            self.assertEqual(prob, 1)


class TestSkipperNestedSequenceCounting(unittest.TestCase):
    """
    The masking isn't Loop-specific: a nested Sequence that goes entirely
    unwitnessed gets the same single lump Skip on the Sequence node, and
    its children need the same containment-aware fallback in
    _conditional_skip_prob to be recognised as skipped.
    """

    @classmethod
    def setUpClass(cls):
        o = Activity(None, 'o', COST); o.id = 'o'
        sched = Activity(None, 'sched', COST); sched.id = 'sched'
        p = Activity(None, 'p', COST); p.id = 'p'
        cls.x = Activity(None, 'x', COST); cls.x.id = 'x'
        cls.y = Activity(None, 'y', COST); cls.y.id = 'y'
        inner = Sequence(None, [cls.x, cls.y])
        cls.x.set_parent(inner); cls.y.set_parent(inner)
        inner.id = 'inner_seq'

        tree = Sequence(None, [o, inner, sched, p])
        o.set_parent(tree); inner.set_parent(tree); sched.set_parent(tree); p.set_parent(tree)
        tree.id = 'root'
        cls.tree = tree

        Aligner.set_level_incentive(0)
        aligner = Aligner(tree)
        # o, sched, p are present; the whole inner Sequence(x, y) is absent
        states, _ = aligner.align2(['o', 'sched', 'p'], [COST, COST, COST], all_optimal=True, timeout=30)
        cls.states = states

    def test_nested_sequence_children_counted_as_skipped(self):
        skipper = Skipper()
        for state in self.states:
            for leaf in (self.x, self.y):
                prob = skipper._conditional_skip_prob(leaf, state)
                self.assertEqual(prob, 1)


class TestSkipProbsEndToEndLoopMasking(unittest.TestCase):
    """
    End-to-end regression for the process-voids-reported case: dv.skip_probs
    for a Loop's do-child, when the whole loop goes unwitnessed while
    surrounding siblings are present. DerivationPipeline.prob_per_variant_and_node
    and prob_per_node (derivation.py) had their own, separate instance of the
    same 0/0 misinterpretation that Skipper._conditional_skip_prob's fix
    alone didn't reach, since both called count_skip_executions/
    count_non_skip_executions directly instead of going through it.
    """

    @classmethod
    def setUpClass(cls):
        o = Activity(None, 'o', COST); o.id = 'o'
        cls.a = Activity(None, 'a', COST); cls.a.id = 'a'
        e = Activity(None, 'e', COST); e.id = 'e'
        sched = Activity(None, 'sched', COST); sched.id = 'sched'
        p = Activity(None, 'p', COST); p.id = 'p'
        loop = Loop(None, [cls.a, e])
        cls.a.set_parent(loop); e.set_parent(loop)
        loop.id = 'loop'
        tree = Sequence(None, [o, loop, sched, p])
        o.set_parent(tree); loop.set_parent(tree); sched.set_parent(tree); p.set_parent(tree)
        tree.id = 'root'

        log = pd.DataFrame({
            'case:concept:name': [1, 1, 1],
            'concept:name': ['o', 'sched', 'p'],
            'time:timestamp': [pd.Timestamp(year=2020, month=1, day=i + 1) for i in range(3)],
        })
        # single trace variant, so pl/pn_measure only need to cover it and
        # its one realized model path: o, then the loop's cheapest
        # execution (one occurrence of do-child 'a'), then sched, p
        pl = {('o', 'sched', 'p'): 1.0}
        pn_measure = {('o', 'a', 'sched', 'p'): 1.0}

        derivation = DerivationPipeline(tree, log, pl=pl, pn_measure=pn_measure)
        cls._original_cwd = os.getcwd()
        cls._tmp_dir = tempfile.mkdtemp(prefix="skipalignments_test_")
        os.chdir(cls._tmp_dir)
        derivation.compute("out")
        os.chdir(cls._original_cwd)
        cls.derivation = derivation

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls._tmp_dir, ignore_errors=True)

    def test_loop_do_child_skip_prob_end_to_end(self):
        self.assertAlmostEqual(self.derivation.skip_probs[self.a], 1.0, places=6)


if __name__ == '__main__':
    unittest.main()
