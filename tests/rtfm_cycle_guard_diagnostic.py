"""
One-off diagnostic (not a unit test -- deliberately not named test*.py so
`python -m unittest discover` skips it) for the id_loop_list/tau_ids cycle
guard against a real log, per process-voids' report that the guard fix
(alignall.py's insert_cycle_checks/does_allow_tau_path/is_cycling_exec
wiring) had zero measurable effect on their rtfm.xes sweep (42/48
variants still clustering at ~100-104s, before and after).

Answers three questions:
  1. Is id_loop_list actually non-empty for rtfm's discovered tree?
  2. Does is_cycling ever fire (return True) during a real search on one
     of the slow variants?
  3. Does the discovered tree contain an And (parallel) node -- if so,
     the 4115-raw-ties number process-voids reported looks more like
     concurrent-branch interleaving explosion (a different, unguarded
     problem on this classical-alignment path) than tau-cycling.

Run with:
    skip/Scripts/python.exe tests/rtfm_cycle_guard_diagnostic.py
"""
import time

import pm4py

import skipalignments.alignall as alignall
from skipalignments.probabilities import EbiOccurance
from skipalignments.processtree import And, ProcessTree

LOG_PATH = "C:/working/data/rtfm.xes"


def has_and_node(tree):
    if isinstance(tree, And):
        return True
    return any(has_and_node(c) for c in tree.children)


def describe(tree, depth=0):
    label = type(tree).__name__
    name = getattr(tree, 'name', None)
    print("  " * depth + label + (f"({name})" if name else ""))
    for c in tree.children:
        describe(c, depth + 1)


def main():
    print(f"Reading log from {LOG_PATH} ...")
    t0 = time.time()
    log = pm4py.read_xes(LOG_PATH)
    print(f"  read in {time.time()-t0:.1f}s, {len(log)} events")

    print("Discovering process tree (inductive miner, noise_threshold=0.0) ...")
    t0 = time.time()
    pm4py_tree = pm4py.discover_process_tree_inductive(log, noise_threshold=0.0)
    print(f"  discovered in {time.time()-t0:.1f}s")

    tree = ProcessTree.from_pm4py(pm4py_tree, 100000, 0, 0)

    print("\n=== Q3: does the discovered tree contain an And node? ===")
    print("And node present:", has_and_node(tree))
    print("\nTree structure:")
    describe(tree)

    print("\n=== Q1: id_loop_list for this tree's build_petri_net ===")
    net, im, fm, activity_to_id, tau_ids, id_loop_list = EbiOccurance().build_petri_net(tree)
    print("len(id_loop_list):", len(id_loop_list))
    for loop_id, loop_pm4py in id_loop_list:
        print(" loop id:", loop_id, "operator:", loop_pm4py.operator)
    print("len(tau_ids):", len(tau_ids))

    print("\n=== Q2: does is_cycling ever fire during a real run? ===")
    variants = log.groupby('case:concept:name')['concept:name'].apply(list)
    distinct = {}
    for seq in variants:
        distinct.setdefault(tuple(seq), 0)
        distinct[tuple(seq)] += 1
    two_activity_variants = [v for v in distinct if len(v) == 2]
    print(f"{len(distinct)} distinct variants total, {len(two_activity_variants)} with exactly 2 activities")

    if not two_activity_variants:
        print("No 2-activity variant found -- picking the shortest available instead.")
        target = min(distinct, key=len)
    else:
        # match process-voids' variant 13 by weight if possible, else just
        # take the first 2-activity variant
        target = two_activity_variants[0]
    print("Target variant:", target, "count:", distinct[target])

    call_count = {"total": 0, "true": 0}
    orig_is_cycling = alignall.is_cycling

    def counting_is_cycling(*args, **kwargs):
        call_count["total"] += 1
        result = orig_is_cycling(*args, **kwargs)
        if result:
            call_count["true"] += 1
        return result

    alignall.is_cycling = counting_is_cycling

    activity_ids = [activity_to_id[a] for a in target]
    t0 = time.time()
    result = alignall.align_pn_all(activity_ids, net, im, fm, id_loop_list, timeout=30, tau_ids=tau_ids)
    elapsed = time.time() - t0

    alignall.is_cycling = orig_is_cycling

    _, (opt_agns, timed_out, _) = result[0]
    print(f"elapsed: {elapsed:.1f}s, timed_out flag: {timed_out}, num optimal alignments returned: {len(opt_agns)}")
    print(f"is_cycling calls: {call_count['total']}, returned True: {call_count['true']}")


if __name__ == '__main__':
    main()
