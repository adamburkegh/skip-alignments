"""
End-to-end skip_probs scenarios that were previously only exercised by
process-voids' own integration test suite -- brought into skip-alignments'
own suite per the house critic review (2026-09-11) "also worth owning"
table, so this package's own test suite is the source of truth for its own
core calculation rather than relying on a downstream consumer's coverage.

Run with:
    python -m unittest tests.test_end_to_end_scenarios -v
"""
import os
import tempfile
import unittest

import pandas as pd

from skipalignments import DerivationPipeline
from skipalignments.processtree import Activity, Sequence, Loop


def _run_in_tmp_dir(derivation, out="out"):
    tmp = tempfile.mkdtemp(prefix="skipalignments_test_")
    orig = os.getcwd()
    os.chdir(tmp)
    try:
        derivation.compute(out)
    finally:
        os.chdir(orig)
    return derivation


class TestSequenceOfSequenceWithFullySkippedCase(unittest.TestCase):
    """
    seq(a, seq(b, c)), log <a,b,c>, <z> (a log-only trace unrelated to the
    model -- forces a fully log-move case where the whole tree is
    skipped).
    """

    @classmethod
    def setUpClass(cls):
        a = Activity(None, 'a', 100000); a.id = 'a'
        b = Activity(None, 'b', 100000); b.id = 'b'
        c = Activity(None, 'c', 100000); c.id = 'c'
        inner = Sequence(None, [b, c]); b.set_parent(inner); c.set_parent(inner); inner.id = 'inner'
        tree = Sequence(None, [a, inner]); a.set_parent(tree); inner.set_parent(tree); tree.id = 'root'
        cls.tree, cls.inner = tree, inner

        log = pd.DataFrame({
            'case:concept:name': [1, 1, 1, 2],
            'concept:name': ['a', 'b', 'c', 'z'],
            'time:timestamp': pd.to_datetime(['2020-01-01', '2020-01-02', '2020-01-03', '2020-01-01']),
        })
        pn_measure = {('a', 'b', 'c'): 1.0}
        derivation = DerivationPipeline(tree, log, pn_measure=pn_measure)
        _run_in_tmp_dir(derivation)
        cls.derivation = derivation

    def test_inner_sequence_has_no_execution_in_the_fully_skipped_case(self):
        # <z> skips the whole root as one block -- inner is masked there,
        # not skipped, so it's excluded from both numerator and
        # denominator (see Skipper.node_reached). That leaves only the
        # <a,b,c> case, where inner genuinely executes and is never
        # skipped -- skip_prob is 0.0, not 0.5.
        self.assertAlmostEqual(self.derivation.skip_probs[self.inner], 0.0, places=6)


class TestSequenceWithSkippedMiddleSubprocess(unittest.TestCase):
    """
    seq(o, seq(b, c), p), log <o,b,c,p> (full execution), <o,p> (middle
    subprocess entirely skipped, o and p still genuinely execute).
    """

    @classmethod
    def setUpClass(cls):
        o = Activity(None, 'o', 100000); o.id = 'o'
        b = Activity(None, 'b', 100000); b.id = 'b'
        c = Activity(None, 'c', 100000); c.id = 'c'
        p = Activity(None, 'p', 100000); p.id = 'p'
        inner = Sequence(None, [b, c]); b.set_parent(inner); c.set_parent(inner); inner.id = 'inner'
        tree = Sequence(None, [o, inner, p])
        o.set_parent(tree); inner.set_parent(tree); p.set_parent(tree)
        tree.id = 'root'
        cls.tree, cls.inner, cls.o, cls.p = tree, inner, o, p

        log = pd.DataFrame({
            'case:concept:name': [1, 1, 1, 1, 2, 2],
            'concept:name': ['o', 'b', 'c', 'p', 'o', 'p'],
            'time:timestamp': pd.to_datetime(
                ['2020-01-01', '2020-01-02', '2020-01-03', '2020-01-04', '2020-01-01', '2020-01-02']
            ),
        })
        pn_measure = {('o', 'b', 'c', 'p'): 1.0}
        derivation = DerivationPipeline(tree, log, pn_measure=pn_measure)
        _run_in_tmp_dir(derivation)
        cls.derivation = derivation

    def test_inner_sequence_skip_prob_is_one_half(self):
        self.assertAlmostEqual(self.derivation.skip_probs[self.inner], 0.5, places=6)

    def test_outer_activities_never_skipped(self):
        self.assertAlmostEqual(self.derivation.skip_probs[self.o], 0.0, places=6)
        self.assertAlmostEqual(self.derivation.skip_probs[self.p], 0.0, places=6)


class TestLoopWithFullySkippedCase(unittest.TestCase):
    """
    loop(a, e), log <a> (one isolated do-iteration, non-skip), <z> (a
    log-only trace -- the whole loop skipped).
    """

    @classmethod
    def setUpClass(cls):
        a = Activity(None, 'a', 100000); a.id = 'a'
        e = Activity(None, 'e', 100000); e.id = 'e'
        loop = Loop(None, [a, e]); a.set_parent(loop); e.set_parent(loop); loop.id = 'loop'
        cls.loop, cls.a = loop, a

        log = pd.DataFrame({
            'case:concept:name': [1, 2],
            'concept:name': ['a', 'z'],
            'time:timestamp': pd.to_datetime(['2020-01-01', '2020-01-01']),
        })
        pn_measure = {('a',): 1.0}
        derivation = DerivationPipeline(loop, log, pn_measure=pn_measure)
        _run_in_tmp_dir(derivation)
        cls.derivation = derivation

    def test_loop_skip_prob_is_one_half(self):
        self.assertAlmostEqual(self.derivation.skip_probs[self.loop], 0.5, places=6)

    def test_loop_do_child_has_no_execution_in_the_fully_skipped_case(self):
        # <z> skips the whole loop as one block -- its do-part 'a' is
        # masked there (no execution at all), not itself skipped, so
        # that case is excluded from both numerator and denominator.
        # That leaves only the <a> case, where 'a' genuinely executes and
        # is never skipped -- skip_prob is 0.0, not 0.5.
        self.assertAlmostEqual(self.derivation.skip_probs[self.a], 0.0, places=6)


class TestSixTracePaymentExample(unittest.TestCase):
    """
    A small payment-process model with six distinct traces, covering a mix
    of fully-synchronous, partially-skipped, and fully-skipped cases in one
    log -- the closest analogue in this suite to the kind of realistic,
    multi-trace fixture process-voids exercises via its own end-to-end
    integration tests.

    Model: seq(order, seq(pay, confirm), ship).
    Traces (1 case each):
      <order,pay,confirm,ship>  -- full synchronous execution
      <order,pay,confirm,ship>  -- (repeat, different case)
      <order,ship>              -- payment subprocess entirely skipped
      <order,pay,ship>          -- confirm alone masked inside a skipped
                                    pay/confirm move on 'confirm' (log-only
                                    'pay' still synchronises)
      <x>                       -- unrelated log-only trace: whole tree
                                    skipped
      <order,pay,confirm,ship>  -- (repeat, third full execution)
    """

    @classmethod
    def setUpClass(cls):
        order = Activity(None, 'order', 100000); order.id = 'order'
        pay = Activity(None, 'pay', 100000); pay.id = 'pay'
        confirm = Activity(None, 'confirm', 100000); confirm.id = 'confirm'
        ship = Activity(None, 'ship', 100000); ship.id = 'ship'
        payment = Sequence(None, [pay, confirm]); pay.set_parent(payment); confirm.set_parent(payment)
        payment.id = 'payment'
        tree = Sequence(None, [order, payment, ship])
        order.set_parent(tree); payment.set_parent(tree); ship.set_parent(tree)
        tree.id = 'root'
        cls.tree, cls.payment, cls.order, cls.ship = tree, payment, order, ship

        cases = (
            [1] * 4 + [2] * 4 + [3] * 2 + [4] * 3 + [5] * 1 + [6] * 4
        )
        events = (
            ['order', 'pay', 'confirm', 'ship'] * 2
            + ['order', 'ship']
            + ['order', 'pay', 'ship']
            + ['x']
            + ['order', 'pay', 'confirm', 'ship']
        )
        log = pd.DataFrame({
            'case:concept:name': cases,
            'concept:name': events,
            'time:timestamp': pd.to_datetime(
                [f'2020-01-0{i}' for i in range(1, 5)] * 2
                + [f'2020-01-0{i}' for i in range(1, 3)]
                + [f'2020-01-0{i}' for i in range(1, 4)]
                + ['2020-01-01']
                + [f'2020-01-0{i}' for i in range(1, 5)]
            ),
        })
        pn_measure = {('order', 'pay', 'confirm', 'ship'): 1.0}
        derivation = DerivationPipeline(tree, log, pn_measure=pn_measure)
        _run_in_tmp_dir(derivation)
        cls.derivation = derivation

    def test_order_and_ship_never_skipped(self):
        # order/ship execute (genuinely or as a direct skip) in every case
        # including <x>, where the whole root is skipped as one block --
        # order/ship are masked there, not counted, so their skip_prob is
        # driven only by the 5 cases where the tree itself has an
        # execution, none of which skip order/ship individually.
        self.assertAlmostEqual(self.derivation.skip_probs[self.order], 0.0, places=6)
        self.assertAlmostEqual(self.derivation.skip_probs[self.ship], 0.0, places=6)

    def test_payment_subprocess_skip_prob(self):
        # payment has an execution in 5 of 6 cases (masked, not skipped,
        # in the <x> case where the whole root is skipped as one block);
        # of those 5, it is skipped outright in exactly 1 (<order,ship>).
        self.assertAlmostEqual(self.derivation.skip_probs[self.payment], 1 / 5, places=6)


if __name__ == '__main__':
    unittest.main()
