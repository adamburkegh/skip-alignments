"""
Unit tests for Aligner.align2()'s rename to align_normal_form(): align2()
was a confusing name (it doesn't distinguish itself from classical
alignment, or hint at the normal-form subtree lumping it performs -- see
alignall.align_pn_all/align_pn_all_multi for the classical alternative with
no lumping). align2() is kept as a backward-compatible alias (pending
deprecation, not deprecated outright) so collaborators depending on an
ancestor project's API don't break.

Run with:
    python -m unittest tests.test_alignment_deprecation -v
"""
import unittest
import warnings

from skipalignments.alignment import Aligner

from tests.test_run_example import build_example_tree

COST = 100000


class TestAlign2DeprecatedAlias(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.tree = build_example_tree()
        Aligner.set_level_incentive(0)

    def test_align2_warns_pending_deprecation(self):
        aligner = Aligner(self.tree)
        with self.assertWarns(PendingDeprecationWarning):
            aligner.align2(['a', 'c'], [COST, COST], all_optimal=True, timeout=30)

    def test_align_normal_form_does_not_warn(self):
        aligner = Aligner(self.tree)
        with warnings.catch_warnings():
            warnings.simplefilter("error", PendingDeprecationWarning)
            aligner.align_normal_form(['a', 'c'], [COST, COST], all_optimal=True, timeout=30)

    def test_align2_delegates_to_align_normal_form(self):
        # same inputs, same aligner state (each gets its own instance since
        # align2/align_normal_form share no cross-call state) -- results
        # must match exactly, not just "both succeed"
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", PendingDeprecationWarning)
            old_states, _ = Aligner(self.tree).align2(['a', 'c'], [COST, COST], all_optimal=True, timeout=30)
        new_states, _ = Aligner(self.tree).align_normal_form(['a', 'c'], [COST, COST], all_optimal=True, timeout=30)

        self.assertEqual(len(old_states), len(new_states))
        self.assertEqual(
            sorted(str(s.path) for s in old_states),
            sorted(str(s.path) for s in new_states),
        )


if __name__ == '__main__':
    unittest.main()
