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
    of its ancestors. When a whole compound subtree goes unwitnessed as a
    block, the Skip wraps the ancestor once instead -- and the leaf's own
    skip/non-skip counts read as 0/0.

    That 0/0 is correct, not a bug: per "Skip Probabilities for
    Subprocesses" (ICPM 2025)'s own Definition of Execution and worked
    example, a node masked inside an ancestor's lump has NO EXECUTION at
    all in that alignment -- not "reached and skipped" -- even when other
    activities genuinely synchronize elsewhere in the same alignment
    ("No execution of the children N10 and N11 of N5 exists in sagn_21:
    Since N5 is skipped, the nodes N10 and N11 are not traversed in
    sagn_21", despite sagn_21 = <a/a, >>/s(N5), >>/s(N6), f/f, g/g> having
    genuine synchronous moves a, f, g elsewhere). A prior fix here
    (commit b840ce9) treated masked-but-alignment-has-other-syncs as
    "skipped" -- confirmed against the paper text to be wrong; see
    CHANGELOG.md for the full writeup.
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

    def test_loop_do_child_has_no_execution_when_whole_loop_unwitnessed(self):
        skipper = Skipper()
        for state in self.states:
            prob = skipper._conditional_skip_prob(self.do_child, state)
            # the loop's do-child never occurs in the trace at all, and the
            # whole loop is masked inside one lump Skip on the Loop node
            # itself -- it has no execution in this alignment (per the
            # paper), which _conditional_skip_prob represents as 0, not 1.
            self.assertEqual(prob, 0)

    def test_loop_redo_child_also_has_no_execution(self):
        # same reasoning as the do-child, for the loop's other child
        skipper = Skipper()
        for state in self.states:
            prob = skipper._conditional_skip_prob(self.redo_child, state)
            self.assertEqual(prob, 0)


class TestSkipperNestedSequenceCounting(unittest.TestCase):
    """
    The masking isn't Loop-specific: a nested Sequence that goes entirely
    unwitnessed gets the same single lump Skip on the Sequence node. Its
    children have no execution in that alignment either -- see
    TestSkipperLoopDoChildCounting's docstring for the paper citation.
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

    def test_nested_sequence_children_have_no_execution(self):
        skipper = Skipper()
        for state in self.states:
            for leaf in (self.x, self.y):
                prob = skipper._conditional_skip_prob(leaf, state)
                self.assertEqual(prob, 0)


class TestSkipProbsEndToEndLoopMasking(unittest.TestCase):
    """
    End-to-end: dv.skip_probs for a Loop's do-child, in a log where the
    whole loop is masked in its one and only trace variant (no other
    variant ever genuinely executes the loop). Per the paper (see
    TestSkipperLoopDoChildCounting's docstring), the do-child has no
    execution anywhere in this log -- prob_per_node's existing "no variant
    reaches this node" fallback correctly defaults that to 0, the same
    fallback it already uses when a node is absent from the tree/log
    entirely, not to 1.
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
        self.assertAlmostEqual(self.derivation.skip_probs[self.a], 0.0, places=6)


if __name__ == '__main__':
    unittest.main()
