"""
Unit tests for EbiOccurance.build_petri_net's transition-label
substitution loop, specifically the latent bug fixed alongside the
insert_cycle_checks wiring (see CHANGELOG.md, tests/test_alignall_cycle_guard.py):
every transition label with no match in the original tree (pm4py's own
invisible routing transitions, label None) was aliased through the same
`None` dict key -- harmless while None was the only such label (aliasing
None to None is a no-op), but insert_cycle_checks' new "TAU_entry_"/
"TAU_exit_" sentinel labels are a second, distinct kind of unmatched
label: the old code would alias the second one encountered to whichever
unmatched label was seen *first*, silently corrupting it.

Run with:
    python -m unittest tests.test_probabilities_build_petri_net -v
"""
import unittest

from skipalignments.probabilities import EbiOccurance
from skipalignments.processtree import Activity, Loop, Tau

from tests.test_alignall_cycle_guard import _tree_with_free_loop_cycle


class TestUnmatchedTransitionLabelsStayDistinct(unittest.TestCase):

    def test_entry_and_exit_sentinels_keep_their_own_distinct_labels(self):
        tree = _tree_with_free_loop_cycle()
        net, im, fm, activity_to_id, tau_ids, id_loop_list = EbiOccurance().build_petri_net(tree)
        self.assertEqual(len(id_loop_list), 1, "fixture must actually trigger insert_cycle_checks")
        loop_id, _ = id_loop_list[0]
        entry_label = "TAU_entry_" + loop_id
        exit_label = "TAU_exit_" + loop_id

        labels = [t.label for t in net.transitions]

        # under the old bug, whichever of these two was encountered first
        # in net.transitions would "win", and the second would be
        # silently overwritten to match it -- so this pins both halves of
        # the corruption directly: the second sentinel actually exists...
        self.assertIn(entry_label, labels)
        self.assertIn(exit_label, labels)
        # ...and neither one exists more than once (i.e. the other wasn't
        # aliased onto it)
        self.assertEqual(labels.count(entry_label), 1)
        self.assertEqual(labels.count(exit_label), 1)

    def test_genuine_pm4py_invisible_transitions_stay_none_labelled(self):
        # Loop(Activity, Tau): the do-child isn't free (a real activity),
        # so this doesn't trigger insert_cycle_checks (no TAU_entry/
        # TAU_exit sentinels) -- but pm4py's own Petri net conversion of
        # a Loop still inserts its own invisible (label=None) helper
        # transitions regardless, confirmed empirically for this shape.
        # Those must stay label=None (not get aliased to some other
        # unmatched label) for net_model_move's `t.label is None` check
        # to keep working.
        act = Activity(None, 'act', 100000)
        act.id = "1"
        redo_tau = Tau(None, 'redo_tau', 0)
        redo_tau.id = "2"
        loop = Loop(None, [act, redo_tau])
        act.set_parent(loop)
        redo_tau.set_parent(loop)
        loop.id = "3"

        net, im, fm, activity_to_id, tau_ids, id_loop_list = EbiOccurance().build_petri_net(loop)
        self.assertEqual(id_loop_list, [], "fixture must NOT trigger insert_cycle_checks -- testing plain pm4py invisibles")

        none_labelled = [t for t in net.transitions if t.label is None]
        self.assertGreater(len(none_labelled), 0, "fixture must actually produce a pm4py invisible transition")


if __name__ == '__main__':
    unittest.main()
