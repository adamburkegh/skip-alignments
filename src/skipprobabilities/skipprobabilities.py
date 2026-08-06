"""
Aggregator module that exposes the full skip-probability derivation
(process trees, alignments, skip alignments, and the derivation pipeline)
under a single import

Usage:
    from skipprobabilities import *
"""
import random
import pm4py
import pandas as pd

from processtree import *
from alignment import *
from execution import *
from probabilities import *
from skips import *
from alignall import *
from derivation import *
from logs import Logs

__all__ = [name for name in dir() if not name.startswith("_")]

def update_pair_taus(tree: "ProcessTree"):
    if isinstance(tree, Tau):
        if tree.parent is not None and len(tree.parent.children) == 2:
            other = tree.parent.children[0]
            if other == tree:
                other = tree.parent.children[1]
            if isinstance(other, Activity):
                tree.name = "TAU_" + other.name
            else:
                tree.name = "TAU_" + other.id
        else:
            tree.name = "TAU_" + str(tree.get_distance_to_root()) + str(random.random())
        return
    elif not isinstance(tree, Activity):
        for c in tree.children:
            update_pair_taus(c)
        return


def check_names(tree: "ProcessTree", names):
    """
    Asserts that every Activity leaf label occurring in `tree` is contained
    in `names` (typically the set of activity labels of the log the tree
    is meant to describe).
    """
    if isinstance(tree, Activity):
        assert tree.name in names
    elif isinstance(tree, Tau):
        pass
    else:
        for c in tree.children:
            check_names(c, names)


def get_variant_dict(log):
    """Returns {variant_tuple: trace_count}, sorted by descending count."""
    variants = dict()
    for k, v in pm4py.statistics.variants.log.get.get_variants_from_log_trace_idx(log).items():
        variants[k] = len(v)
    return dict(sorted(variants.items(), key=lambda x: -x[1]))


def get_activities(log):
    """Returns the list of distinct activity labels appearing in the log's variants."""
    variants = get_variant_dict(log)
    activities = []
    for var in variants.keys():
        for act in list(var):
            if act not in activities:
                activities.append(act)
    return activities


def generate_tree(activities, prob_sequence=0.25, prob_xor=0.25, prob_and=0.25,
                   prob_loop=0.25, prob_tau=0.4, max_children=4):
    operator_r = random.random()
    num_children = random.randint(2, max_children)
    if operator_r > prob_sequence + prob_xor + prob_and + prob_loop:
        # create a leaf node
        if random.random() < prob_tau or len(activities) == 0:
            return Tau(None, 'TAU', 0)
        else:
            return Activity(None, activities[random.randint(0, len(activities) - 1)], 100000)
    else:
        children = []
        while len(children) < num_children:
            c = generate_tree(
                [act for act in activities if act not in [x.name for x in children if isinstance(x, Activity)]],
                prob_sequence / 2, prob_xor / 2, prob_and / 2, prob_loop / 2, prob_tau, max_children,
            )
            if isinstance(c, Tau) and sum(isinstance(x, Tau) for x in children) > (
                    0 if operator_r < prob_sequence + prob_xor + prob_and else 1):
                # no two taus in non-loops
                continue
            if isinstance(c, Activity) and sum(isinstance(x, Activity) and x.name == c.name for x in children) > 0:
                # no duplicate label on same leafs
                continue
            children.append(c)
        if operator_r < prob_sequence:
            node = Sequence(None, children)
        elif operator_r < prob_sequence + prob_xor:
            node = Xor(None, children)
        elif operator_r < prob_sequence + prob_xor + prob_and:
            node = And(None, children)
        else:
            children = children[:2]
            node = Loop(None, children)
        for c in children:
            c.set_parent(node)
        return node
