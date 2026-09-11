"""
Additional regression tests for the masked-execution skip_probs fix (see
CHANGELOG.md, Skipper.node_reached), contributed by process-voids' house
critic review (2026-09-11) of that fix. Each case is designed to catch a
specific failure mode the original fix's own tests couldn't distinguish.

Run with:
    python -m unittest tests.test_masked_execution_review_cases -v
"""
import os
import tempfile
import unittest

import pandas as pd

from skipalignments import DerivationPipeline
from skipalignments.processtree import Activity, Sequence, Xor

from tests.test_run_example import build_example_tree


def _run_in_tmp_dir(derivation, out="out"):
    tmp = tempfile.mkdtemp(prefix="skipalignments_test_")
    orig = os.getcwd()
    os.chdir(tmp)
    try:
        derivation.compute(out)
    finally:
        os.chdir(orig)
    return derivation


class TestMaskingTwoLevelsDeepWithRenormalisation(unittest.TestCase):
    """
    Reviewer's case 1 -- "the one the existing tests don't cover": a
    fixture where every masked node's *expected* value happens to be 0
    can't distinguish a correctly-restricted denominator (masked traces
    dropped from both numerator and denominator) from an unrestricted one
    (masked traces dropped from the numerator but wrongly kept in the
    denominator, diluting it) -- both give 0 either way. This fixture has
    a masked node (Y) whose correct value is a nonzero 1/3, which an
    unrestricted denominator would instead compute as 1/4.

    Model: seq(a, X, e) with X = seq(b, Y), Y = seq(c, d).
    Log: <a,b,c,d,e> x2, <a,b,e> x1, <a,e> x1.

    <a,b,c,d,e>: all synchronous -- X and Y both executed, non-skip.
    <a,b,e>: (a,a)(b,b)(>>,s(Y))(e,e) -- Y skipped (direct, own move);
             c and d masked (no execution at all).
    <a,e>: (a,a)(>>,s(X))(e,e) -- X skipped (direct, own move); b, Y, c,
           d all masked (no execution at all).
    """

    @classmethod
    def setUpClass(cls):
        a = Activity(None, 'a', 100000); a.id = 'a'
        b = Activity(None, 'b', 100000); b.id = 'b'
        c = Activity(None, 'c', 100000); c.id = 'c'
        d = Activity(None, 'd', 100000); d.id = 'd'
        e = Activity(None, 'e', 100000); e.id = 'e'
        Y = Sequence(None, [c, d]); c.set_parent(Y); d.set_parent(Y); Y.id = 'Y'
        X = Sequence(None, [b, Y]); b.set_parent(X); Y.set_parent(X); X.id = 'X'
        tree = Sequence(None, [a, X, e]); a.set_parent(tree); X.set_parent(tree); e.set_parent(tree); tree.id = 'root'

        cls.tree, cls.a, cls.b, cls.c, cls.d, cls.e, cls.X, cls.Y = tree, a, b, c, d, e, X, Y

        log = pd.DataFrame({
            'case:concept:name': [1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 4, 4],
            'concept:name': ['a', 'b', 'c', 'd', 'e'] * 2 + ['a', 'b', 'e'] + ['a', 'e'],
            'time:timestamp': pd.to_datetime(
                [f'2020-01-0{i}' for i in range(1, 6)] * 2
                + [f'2020-01-0{i}' for i in range(1, 4)]
                + [f'2020-01-0{i}' for i in range(1, 3)]
            ),
        })
        pl = {('a', 'b', 'c', 'd', 'e'): 0.5, ('a', 'b', 'e'): 0.25, ('a', 'e'): 0.25}
        pn_measure = {('a', 'b', 'c', 'd', 'e'): 1.0}

        derivation = DerivationPipeline(tree, log, pl=pl, pn_measure=pn_measure)
        _run_in_tmp_dir(derivation)
        cls.derivation = derivation

    def test_X_skip_prob_is_one_quarter(self):
        self.assertAlmostEqual(self.derivation.skip_probs[self.X], 0.25, places=6)

    def test_Y_skip_prob_is_one_third_not_one_quarter(self):
        # 1/3: Y has an execution in 3 of the 4 cases (2x fully-synchronous,
        # 1x directly skipped), 1 skip among them -- the <a,e> case (where
        # Y is masked, not skipped) must be excluded from the denominator
        # too. An unrestricted denominator (masked traces kept in the
        # denominator) would instead give 1/4.
        self.assertAlmostEqual(self.derivation.skip_probs[self.Y], 1 / 3, places=6)

    def test_masked_descendants_and_root_are_zero(self):
        for node, name in [(self.tree, 'root'), (self.a, 'a'), (self.b, 'b'), (self.c, 'c'), (self.d, 'd'), (self.e, 'e')]:
            with self.subTest(node=name):
                self.assertAlmostEqual(self.derivation.skip_probs[node], 0.0, places=6)


class TestPaperExampleWithNaturalLogFrequencies(unittest.TestCase):
    """
    Reviewer's case 2: the paper's own running example (fig:details /
    fig:en in "Skip Probabilities for Subprocesses", ICPM 2025), but
    driven by natural log-frequency `pl` (derived from case counts) rather
    than an explicit override -- the existing paper-fixture test
    (test_run_example.test_skip_prob_matches_paper_worked_example) always
    passes an explicit `pl` override, which never exercises
    DerivationPipeline's own frequency-counting path.

    10 cases: 1x <b>, 2x <a,f,e>, 7x <a,c,e> -- giving exact natural
    frequencies P_L(sigma1)=0.1, P_L(sigma2)=0.2, P_L(sigma3)=0.7,
    matching the paper's own numbers exactly (no rounding needed).

    Paper: "P_{N_5}(sigma_2) = P_L(sigma_2) * 1/(P_L(sigma_2)+P_L(sigma_3))
    = 0.22" and P_{N_5}(sk) = P_{N_5}(e_21)+P_{N_5}(e_22) = 0.22 (rounded);
    exactly, 0.2/(0.2+0.7) = 2/9. N_1 (the tree root) is skipped in
    exactly the <b> cases (log move on 'b', model move s(N_1)) -- its
    skip_prob is the plain P_L(sigma1) = 0.1, since the root has an
    execution (skipped or not) in every case.
    """

    @classmethod
    def setUpClass(cls):
        cls.tree = build_example_tree()
        sequence, loop = cls.tree.children
        a, choice = sequence.children
        cls.choice = choice  # N_5

        cases = [1] + [cid for cid in range(2, 4) for _ in range(3)] + [cid for cid in range(4, 11) for _ in range(3)]
        events = ['b'] + ['a', 'f', 'e'] * 2 + ['a', 'c', 'e'] * 7
        log = pd.DataFrame({
            'case:concept:name': cases,
            'concept:name': events,
            'time:timestamp': [pd.Timestamp(year=1000 + i, month=1, day=1) for i in range(len(cases))],
        })
        model_dist = {
            ('4', '8', '6', '7', '6'): 0.1,
            ('4', '9', '6', '7', '6'): 0.1,
            ('4', '8', '6'): 0.3,
            ('4', '9', '6'): 0.3,
        }

        derivation = DerivationPipeline(cls.tree, log, pn_measure=model_dist)
        _run_in_tmp_dir(derivation)
        cls.derivation = derivation

    def test_natural_pl_matches_expected_case_frequencies(self):
        # sanity check on the fixture itself before trusting downstream
        # numbers: pl must genuinely be derived from case counts, not
        # something we accidentally still overrode.
        total = sum(self.derivation.pl.values())
        self.assertAlmostEqual(total, 1.0, places=6)
        self.assertAlmostEqual(self.derivation.pl.get(('b',), 0.0), 0.1, places=6)
        self.assertAlmostEqual(self.derivation.pl.get(('a', 'f', 'e'), 0.0), 0.2, places=6)
        self.assertAlmostEqual(self.derivation.pl.get(('a', 'c', 'e'), 0.0), 0.7, places=6)

    def test_choice_node_skip_prob_is_two_ninths(self):
        self.assertAlmostEqual(self.derivation.skip_probs[self.choice], 2 / 9, places=6)

    def test_root_skip_prob_is_one_tenth(self):
        self.assertAlmostEqual(self.derivation.skip_probs[self.tree], 0.1, places=6)


class TestDistinctLeavesWithSameLabelEq9Restriction(unittest.TestCase):
    """
    Reviewer's case 3 (marked optional): exercises Eq. (9)'s restriction of
    a trace's skip alignments to those that give a node an execution, using
    two distinct model leaves that share a log activity label ('a').

    Model: xor(seq(a1,b), seq(a2,c)) -- two separate 'a' leaves, one per
    branch. Log: <a>, <a,b>, <a,c>, one case each.

    <a>: whichever xor branch is chosen, its own trailing activity (b or c)
    is a genuine model-only move (skipped), and the OTHER branch is a
    genuine log-only/model-skip that never executes at all -- so this case
    contributes a skip to whichever one of b/c belongs to the branch that
    "won" the optimal alignment, and nothing to the other (masked, no
    execution).
    <a,b>: branch 1 taken, b executes (non-skip); branch 2 (a2,c) has no
    execution at all (masked -- the whole xor branch not chosen).
    <a,c>: branch 2 taken, c executes (non-skip); branch 1 has no execution.

    Each of b, c has an execution in exactly 2 of the 3 cases (once
    skipped, once genuinely executed) -- skip_prob(b) = skip_prob(c) = 1/2.
    """

    @classmethod
    def setUpClass(cls):
        a1 = Activity(None, 'a', 100000); a1.id = 'a1'
        b = Activity(None, 'b', 100000); b.id = 'b'
        a2 = Activity(None, 'a', 100000); a2.id = 'a2'
        c = Activity(None, 'c', 100000); c.id = 'c'
        branch1 = Sequence(None, [a1, b]); a1.set_parent(branch1); b.set_parent(branch1); branch1.id = 'branch1'
        branch2 = Sequence(None, [a2, c]); a2.set_parent(branch2); c.set_parent(branch2); branch2.id = 'branch2'
        tree = Xor(None, [branch1, branch2]); branch1.set_parent(tree); branch2.set_parent(tree); tree.id = 'root'

        cls.tree, cls.b, cls.c = tree, b, c

        log = pd.DataFrame({
            'case:concept:name': [1, 2, 2, 3, 3],
            'concept:name': ['a', 'a', 'b', 'a', 'c'],
            'time:timestamp': pd.to_datetime(['2020-01-01', '2020-01-01', '2020-01-02', '2020-01-01', '2020-01-02']),
        })

        pn_measure = {('a1', 'b'): 0.5, ('a2', 'c'): 0.5}
        derivation = DerivationPipeline(tree, log, pn_measure=pn_measure)
        _run_in_tmp_dir(derivation)
        cls.derivation = derivation

    def test_b_skip_prob_is_one_half(self):
        self.assertAlmostEqual(self.derivation.skip_probs[self.b], 0.5, places=6)

    def test_c_skip_prob_is_one_half(self):
        self.assertAlmostEqual(self.derivation.skip_probs[self.c], 0.5, places=6)


if __name__ == '__main__':
    unittest.main()
