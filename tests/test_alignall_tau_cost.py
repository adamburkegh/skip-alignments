"""
Unit tests for alignall's classical (Petri-net-level) all-optimal-alignment
functions pricing a genuine ProcessTree Tau leaf correctly.

net_model_move (repeated in align_pn_all/align_pn_all_multi/
align_pn_all_for_one/align_pn_one_for_one) decides a model move is free via
`t.label is None or t.label.startswith("TAU")`. That convention matches
ProcessTree.to_pm4py(use_ids=False), where a Tau leaf is labelled
"TAU_"+name, and matches pm4py's own auto-inserted invisible routing
transitions (label None). It does NOT survive
EbiOccurance.build_petri_net, which always calls to_pm4py(use_ids=True):
there, a Tau leaf keeps its own opaque tree id (e.g. "3"), same as any
Activity, so net_model_move silently prices choosing the model's own free
Tau branch at the full model-move cost (100000) instead of 0. Found via a
process-voids report (running-example fixture, payment_approval's
schedule_choice = Xor(s, Tau)): a variant that legitimately took the tau
branch showed 1 deficit rather than 0.

Run with:
    python -m unittest tests.test_alignall_tau_cost -v
"""
import unittest

from skipalignments.alignall import align_pn_all
from skipalignments.probabilities import EbiOccurance
from skipalignments.processtree import Activity, Sequence, Tau, Xor


def _tree_with_tau_branch():
    """Sequence(s, Xor(a, Tau)): trace ['s'] can only be completed by the
    model firing its own Tau branch of the Xor (the other branch, 'a', has
    no matching log event and costs strictly more to fire)."""
    s = Activity(None, 's', 100000)
    s.id = "1"
    a = Activity(None, 'a', 100000)
    a.id = "2"
    tau = Tau(None, 'tau', 0)
    tau.id = "3"

    choice = Xor(None, [a, tau])
    a.set_parent(choice)
    tau.set_parent(choice)
    choice.id = "4"

    tree = Sequence(None, [s, choice])
    s.set_parent(tree)
    choice.set_parent(tree)
    tree.id = "5"
    return tree


class TestAlignPnAllTauLeafCost(unittest.TestCase):

    def test_choosing_the_models_own_tau_branch_costs_zero(self):
        tree = _tree_with_tau_branch()
        net, im, fm, activity_to_id, tau_ids, id_loop_list = EbiOccurance().build_petri_net(tree)

        result = align_pn_all([activity_to_id['s']], net, im, fm, [], timeout=30, tau_ids=tau_ids)
        _, (opt_agns, timed_out, _) = result[0]

        self.assertEqual(timed_out, 0)
        self.assertGreater(len(opt_agns), 0)
        # the model legitimately took its own free tau branch -- deficit
        # (cost) must be 0, not 100000 (the mispriced "fire a real
        # activity" cost that net_model_move currently assigns to it)
        self.assertEqual(opt_agns[0]['cost'], 0)

    def test_without_tau_ids_the_old_mispricing_still_applies(self):
        # documents the deliberate backward-compatible default: a caller
        # that doesn't pass tau_ids (e.g. building its own net by hand,
        # with pm4py's "TAU_"-prefixed label convention rather than
        # build_petri_net's opaque ids) keeps getting the original
        # string-sniffing behaviour, unchanged.
        tree = _tree_with_tau_branch()
        net, im, fm, activity_to_id, _tau_ids, _id_loop_list = EbiOccurance().build_petri_net(tree)

        result = align_pn_all([activity_to_id['s']], net, im, fm, [], timeout=30)
        _, (opt_agns, timed_out, _) = result[0]

        self.assertEqual(timed_out, 0)
        self.assertEqual(opt_agns[0]['cost'], 100000)


if __name__ == '__main__':
    unittest.main()
