"""
Manual scaling probe (not a unit test -- deliberately not named test*.py
so `python -m unittest discover` skips it) for __search's tie-count
scaling, since the process-voids report this fix responds to was about
real-log scale (thousands of ties) that a fast unit test can't
demonstrate on its own. See tests/test_alignall_tie_explosion.py for a
fast, permanent regression test using the same fixture at a small n.

Confirms __reconstruct_alignment's O(depth) rebuild (was O(depth^2) via
repeated list-prepend) and closed's O(1) average set membership (was
O(n) linear scan over a list) together turn what was catastrophic
super-linear scaling into roughly linear-in-tie-count scaling:

    n=5  ties=120    elapsed=0.06s -> 0.03s
    n=6  ties=720    elapsed=0.99s -> 0.20s
    n=7  ties=5040   elapsed=43.02s -> 1.36s
    n=8  ties=40320  elapsed=(not measured pre-fix) -> 11.07s

(measured on this machine; expect variance elsewhere, but the *shape*
-- from clearly super-linear to roughly linear -- should reproduce).

Run with:
    skip/Scripts/python.exe tests/perf_search_scaling.py
"""
import time

from skipalignments.alignall import align_pn_all, set_max_tie_count
from skipalignments.probabilities import EbiOccurance
from skipalignments.processtree import Activity, And, Tau, Xor


def _tree_with_n_concurrent_optional_branches(n):
    """And of n independent, concurrent, individually-optional branches
    -- against a trace matching one branch's real activity, the cheapest
    resolution has the other n-1 branches take their Tau alternative, and
    all n branches' relative firing order is unconstrained -> n! tied
    optimal alignments. See test_alignall_tie_explosion.py."""
    branches = []
    for i in range(n):
        act = Activity(None, f'act{i}', 100000)
        act.id = f"a{i}"
        tau = Tau(None, f'tau{i}', 0)
        tau.id = f"t{i}"
        choice = Xor(None, [act, tau])
        act.set_parent(choice)
        tau.set_parent(choice)
        choice.id = f"c{i}"
        branches.append(choice)
    tree = And(None, branches)
    for b in branches:
        b.set_parent(tree)
    tree.id = "root"
    return tree


if __name__ == '__main__':
    set_max_tie_count(10_000_000)
    for n in (5, 6, 7, 8):
        tree = _tree_with_n_concurrent_optional_branches(n)
        net, im, fm, activity_to_id, tau_ids, id_loop_list = EbiOccurance().build_petri_net(tree)
        t0 = time.time()
        result = align_pn_all([activity_to_id['act0']], net, im, fm, id_loop_list, timeout=120, tau_ids=tau_ids)
        elapsed = time.time() - t0
        _, (opt_agns, timed_out, _) = result[0]
        print(f"n={n} ties={len(opt_agns)} timed_out={timed_out} elapsed={elapsed:.2f}s")
