"""
Unit tests for alignall's loop cycle-guard (insert_cycle_checks /
id_loop_list), applied to nets built the only way this codebase actually
builds them: EbiOccurance.build_petri_net's to_pm4py(use_ids=True).

Root cause (a process-voids report on a real 104k-case log: 42/48 trace
variants clustered at the ~100s per-variant timeout regardless of variant
size/rarity -- the exact signature of an unbounded free-tau-cycle search):
insert_cycle_checks/does_allow_tau_path/is_cycling_exec were never wired
into build_petri_net, AND, independently, all three string-sniff a "TAU_"
label prefix (ProcessTree.to_pm4py(use_ids=False)'s convention) to
recognize a tau. build_petri_net always uses use_ids=True, under which a
genuine Tau leaf keeps its own opaque tree id instead -- the same
mismatch already found and fixed for net_model_move's cost function (see
test_alignall_tau_cost.py), but here it's the *cycle guard itself* that
silently never engages, not just a cost. Fixed by threading tau_ids
(build_petri_net's own tree walk) through
ProcessTree.from_pm4py/does_allow_tau_path/insert_cycle_checks (to decide
which loops need guarding) and through is_cycling/is_cycling_exec (to
recognize a tau'd do/redo execution during search), the same tau_ids
already returned by build_petri_net for the cost fix.

Run with:
    python -m unittest tests.test_alignall_cycle_guard -v
"""
import unittest

from skipalignments.alignall import get_leafs, is_cycling_exec
from skipalignments.probabilities import EbiOccurance
from skipalignments.processtree import Activity, Loop, Tau, Xor


def _tree_with_free_loop_cycle():
    """
    Loop(do=Xor(act, Tau), redo=Tau): both the do-child and the redo-child
    can execute for free, so the model can alternate do-tau/redo-tau
    forever at zero marginal cost -- the exact does_allow_tau_path
    precondition the cycle guard exists to bound. An empty trace forces
    the search to consider this free cycle, since firing 'act' would only
    ever add cost for no reason.
    """
    act = Activity(None, 'act', 100000)
    act.id = "1"
    do_tau = Tau(None, 'do_tau', 0)
    do_tau.id = "2"
    do_choice = Xor(None, [act, do_tau])
    act.set_parent(do_choice)
    do_tau.set_parent(do_choice)
    do_choice.id = "3"

    redo_tau = Tau(None, 'redo_tau', 0)
    redo_tau.id = "4"

    loop = Loop(None, [do_choice, redo_tau])
    do_choice.set_parent(loop)
    redo_tau.set_parent(loop)
    loop.id = "5"
    return loop


class TestBuildPetriNetDetectsRiskyLoop(unittest.TestCase):
    """Pins the does_allow_tau_path/from_pm4py fix directly and cheaply,
    independent of any actual search timing."""

    def test_id_loop_list_is_non_empty_for_a_free_cycle_loop(self):
        tree = _tree_with_free_loop_cycle()
        net, im, fm, activity_to_id, tau_ids, id_loop_list = EbiOccurance().build_petri_net(tree)
        self.assertEqual(len(id_loop_list), 1)

    def test_id_loop_list_is_empty_when_no_loop_can_go_fully_free(self):
        # redo-child is a real (non-free) activity -- alternating do/redo
        # forever always costs more, so no guard is needed and none
        # should be inserted (matches does_allow_tau_path's own
        # precondition: both do and at least one redo must allow tau)
        act = Activity(None, 'act', 100000)
        act.id = "1"
        do_tau = Tau(None, 'do_tau', 0)
        do_tau.id = "2"
        do_choice = Xor(None, [act, do_tau])
        act.set_parent(do_choice)
        do_tau.set_parent(do_choice)
        do_choice.id = "3"

        redo_act = Activity(None, 'redo_act', 100000)
        redo_act.id = "4"

        loop = Loop(None, [do_choice, redo_act])
        do_choice.set_parent(loop)
        redo_act.set_parent(loop)
        loop.id = "5"

        net, im, fm, activity_to_id, tau_ids, id_loop_list = EbiOccurance().build_petri_net(loop)
        self.assertEqual(id_loop_list, [])


class TestIsCyclingExecRecognizesTaudIterations(unittest.TestCase):
    """
    Direct test of is_cycling_exec's actual regression: it must recognize
    a do/redo execution as "tau" (free, unproductive) via tau_ids when the
    net is ids-labelled -- the old `m.startswith("TAU_")` check never
    matches an opaque tree id, so on a real ids-labelled net it silently
    never detects the cycle at all (confirmed empirically: a full
    align_pn_all run against this exact fixture terminated in
    milliseconds whether or not the guard fired, for a tie this small --
    the reported production cost only shows up at the scale and structure
    of a real log's search space, see process-voids' report referenced
    above). This test isolates the piece that must be correct for the
    guard to have any chance of firing on a real, larger case: recognizing
    two consecutive tau-only do/redo executions.
    """

    def _do_redo_leafs_and_id_loop_list(self):
        tree = _tree_with_free_loop_cycle()
        net, im, fm, activity_to_id, tau_ids, id_loop_list = EbiOccurance().build_petri_net(tree)
        self.assertEqual(len(id_loop_list), 1)
        loop_id, loop_pm4py = id_loop_list[0]
        return loop_id, loop_pm4py, tau_ids

    def test_two_consecutive_tau_executions_are_recognized_as_cycling(self):
        loop_id, loop_pm4py, tau_ids = self._do_redo_leafs_and_id_loop_list()
        # the do-child is Xor(act, do_tau) -- pick its tau alternative
        # specifically, not whichever leaf get_leafs lists first
        do_tau_leaf = [l for l in get_leafs(loop_pm4py.children[0]) if l in tau_ids][0]
        redo_leaf = get_leafs(loop_pm4py.children[1])[0]  # redo-child is a bare Tau, no alternative

        # do(tau), redo(tau), do(tau), still running (no exit sentinel
        # yet) -- two full tau'd executions in a row, exactly the
        # unproductive repeat the guard exists to cut off
        agn = [
            (None, "TAU_entry_" + loop_id),
            ('>>', do_tau_leaf),
            ('>>', redo_leaf),
            ('>>', do_tau_leaf),
        ]
        self.assertTrue(is_cycling_exec(agn, loop_pm4py, loop_id, tau_ids=tau_ids))

    def test_same_sequence_is_missed_without_tau_ids(self):
        # this is the actual bug: on an ids-labelled net, an opaque id
        # like "2" never starts with "TAU_", so every do/redo execution
        # looks like a real (non-tau) activity to the old check --
        # the cycle is never recognized, and align_pn_all has nothing left
        # to bound the search on a real trace/model where this tie
        # genuinely doesn't resolve in milliseconds
        loop_id, loop_pm4py, tau_ids = self._do_redo_leafs_and_id_loop_list()
        # the do-child is Xor(act, do_tau) -- pick its tau alternative
        # specifically, not whichever leaf get_leafs lists first
        do_tau_leaf = [l for l in get_leafs(loop_pm4py.children[0]) if l in tau_ids][0]
        redo_leaf = get_leafs(loop_pm4py.children[1])[0]  # redo-child is a bare Tau, no alternative

        agn = [
            (None, "TAU_entry_" + loop_id),
            ('>>', do_tau_leaf),
            ('>>', redo_leaf),
            ('>>', do_tau_leaf),
        ]
        self.assertFalse(is_cycling_exec(agn, loop_pm4py, loop_id, tau_ids=None))


if __name__ == '__main__':
    unittest.main()
