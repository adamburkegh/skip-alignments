"""
Unit test for a redundant computation in DerivationPipeline.__init__,
found via a process-voids performance report: constructing a
DerivationPipeline over a 103,987-case log took ~24.75s, entirely before
DerivationPipeline.compute() was even called, and independent of the
mass-term alignment stage (compute()'s own "1/4" stage) that's actually
supposed to do the expensive work.

Root cause: __init__ calls self.get_variant_dict(self.aligned_log) once
directly (line ~70, to set self.variants), and -- whenever the caller
doesn't pass pl explicitly, the common case -- ALSO calls
self.variant_prob_dist(self.aligned_log), which internally calls
self.get_variant_dict(log) *again* on the same log to compute the same
thing a second time. get_variant_dict wraps
pm4py.statistics.variants.log.get.get_variants_from_log_trace_idx, an
O(log size) full pass grouping every case into its variant -- doing that
twice roughly doubles that part of construction for no reason, since the
result is identical both times (same log, same method, no side effects
in between).

This doesn't fully account for the reported ~24.75s on its own (a
103,987-case log's variant extraction may simply be expensive once), but
it's a real, confirmed, easily-avoided 2x on whatever that cost is.

Run with:
    python -m unittest tests.test_derivation_construction_redundancy -v
"""
import unittest
from unittest.mock import patch

from skipalignments.derivation import DerivationPipeline

from tests.test_run_example import build_example_log, build_example_tree


class TestConstructorDoesNotComputeVariantsTwice(unittest.TestCase):

    def test_get_variant_dict_called_once_when_pl_not_supplied(self):
        tree = build_example_tree()
        log = build_example_log()

        with patch.object(DerivationPipeline, 'get_variant_dict', wraps=DerivationPipeline.get_variant_dict, autospec=True) as spy:
            DerivationPipeline(tree, log, pn_measure={})

        self.assertEqual(spy.call_count, 1)

    def test_pl_still_correct_when_not_supplied(self):
        # the fix must not change the actual computed distribution --
        # only how many times the underlying pass over the log happens
        tree = build_example_tree()
        log = build_example_log()

        derivation = DerivationPipeline(tree, log, pn_measure={})

        total = sum(derivation.variants.values())
        expected_pl = {k: v / total for k, v in derivation.variants.items()}
        self.assertEqual(derivation.pl, expected_pl)

    def test_explicit_pl_is_still_used_as_is(self):
        tree = build_example_tree()
        log = build_example_log()
        explicit_pl = {('a',): 1.0}

        derivation = DerivationPipeline(tree, log, pl=explicit_pl, pn_measure={})

        self.assertEqual(derivation.pl, explicit_pl)


if __name__ == '__main__':
    unittest.main()
