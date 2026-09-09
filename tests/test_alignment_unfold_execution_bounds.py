"""
Unit test for a real crash reported by process-voids: build_execution_tree's
`assert tree_c is not None` (execution.py) firing on a real
inductive_noise20-discovered rtfm.xes tree, activity_gradual degradation,
trace variant ('Add penalty',) among others.

Root cause is in State.unfold() (alignment.py), not build_execution_tree
itself: when a Skip/TauPath at path position i is unfolded into a
different-length replacement, every *other* recorded Execution whose span
needs to grow to accommodate the length change should have its .stop
extended. The boundary check `start < i and stop > i` uses a strict `<`,
so an execution whose .start is exactly i (the skip is the very *first*
slot inside that execution's own span, not something before it) matches
neither this branch nor the `start > i` branch -- its .stop silently
stays stale (too narrow) whenever the replacement's length differs from
the original 1-slot skip. build_execution_tree then can't find that
execution's own children within its now-too-narrow bounds, hence the
assertion.

This only manifests when the length actually changes (a replacement of
length 1 hides the bug entirely, since stop += 0 is a no-op either way)
-- which is why one of two tied optimal alignments for this exact trace
crashes and the other doesn't (see the tree below: the same skip's two
tied resolutions replace a 1-slot Skip with either a 1-activity path
[Appeal to Judge] or a 2-activity path [Send Fine, Tau] -- only the
2-item replacement exposes the stale bound).

Run with:
    python -m unittest tests.test_alignment_unfold_execution_bounds -v
"""
import unittest

from skipalignments.alignment import Aligner
from skipalignments.execution import ExecutionManager
from skipalignments.processtree import Activity, And, Sequence, Tau, Xor

COST = 100000


def _rtfm_repro_tree():
    """
    Xor(Tau, And(
        Xor(Tau, Activity('Insert Date Appeal to Prefecture')),
        Sequence(
            Xor(
                Activity('Appeal to Judge'),
                Sequence(Activity('Send Fine'), Xor(Tau, Activity('Insert Fine Notification')))
            ),
            Xor(Tau, Activity('Add penalty'))
        )
    ))

    A real tree discovered by pm4py's inductive miner (noise_threshold=0.2)
    from rtfm.xes, as reported by process-voids.
    """
    def act(name):
        a = Activity(None, name, COST)
        a.id = name
        return a

    def tau(name):
        t = Tau(None, name, 0)
        t.id = name
        return t

    insert_date = act('Insert Date Appeal to Prefecture')
    tau1 = tau('tau1')
    xor1 = Xor(None, [tau1, insert_date])
    tau1.set_parent(xor1)
    insert_date.set_parent(xor1)
    xor1.id = 'xor1'

    appeal_judge = act('Appeal to Judge')
    send_fine = act('Send Fine')
    tau2 = tau('tau2')
    insert_notif = act('Insert Fine Notification')
    xor_notif = Xor(None, [tau2, insert_notif])
    tau2.set_parent(xor_notif)
    insert_notif.set_parent(xor_notif)
    xor_notif.id = 'xor_notif'

    seq_send_notif = Sequence(None, [send_fine, xor_notif])
    send_fine.set_parent(seq_send_notif)
    xor_notif.set_parent(seq_send_notif)
    seq_send_notif.id = 'seq_send_notif'

    xor_judge_or_seq = Xor(None, [appeal_judge, seq_send_notif])
    appeal_judge.set_parent(xor_judge_or_seq)
    seq_send_notif.set_parent(xor_judge_or_seq)
    xor_judge_or_seq.id = 'xor_judge_or_seq'

    tau3 = tau('tau3')
    add_penalty = act('Add penalty')
    xor_penalty = Xor(None, [tau3, add_penalty])
    tau3.set_parent(xor_penalty)
    add_penalty.set_parent(xor_penalty)
    xor_penalty.id = 'xor_penalty'

    seq_inner = Sequence(None, [xor_judge_or_seq, xor_penalty])
    xor_judge_or_seq.set_parent(seq_inner)
    xor_penalty.set_parent(seq_inner)
    seq_inner.id = 'seq_inner'

    and_node = And(None, [xor1, seq_inner])
    xor1.set_parent(and_node)
    seq_inner.set_parent(and_node)
    and_node.id = 'and_node'

    tau_root = tau('tau_root')
    root = Xor(None, [tau_root, and_node])
    tau_root.set_parent(root)
    and_node.set_parent(root)
    root.id = 'root'

    return root


class TestUnfoldExtendsStaleExecutionBoundsAtExactBoundary(unittest.TestCase):

    def test_coninciding_agns_does_not_crash_on_add_penalty_variant(self):
        tree = _rtfm_repro_tree()
        Aligner.set_level_incentive(0)
        aligner = Aligner(tree)
        variant = ['Add penalty']
        states, _ = aligner.align_normal_form(variant, [COST] * len(variant), all_optimal=True, timeout=30)
        self.assertEqual(len(states), 2, "fixture must reproduce exactly the reported tie")

        em = ExecutionManager()
        # must not raise
        em.coninciding_agns({','.join(variant): states})


if __name__ == '__main__':
    unittest.main()
