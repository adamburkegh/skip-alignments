"""
Plain-Python equivalent of example.ipynb: computes skip probabilities for
the running example from the paper, without needing Jupyter.

Run with:
    pip install -e .
    python examples/quickstart.py
"""
import sys

import pandas as pd

from skipalignments import Activity, Sequence, Xor, Loop, DerivationPipeline

# print_blinded() below uses Unicode operator glyphs (-> x /\ etc.); some
# terminals (notably Windows' default cp1252 console) can't encode them
# without this.
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")


def build_example_tree():
    a = Activity(None, 'a', 100000)
    a.id = "4"
    c = Activity(None, 'c', 100000)
    c.id = "8"
    d = Activity(None, 'd', 100000)
    d.id = "9"
    e = Activity(None, 'e', 100000)
    e.id = "6"
    f = Activity(None, 'f', 100000)
    f.id = "7"

    choice = Xor(None, [c, d])
    c.set_parent(choice)
    d.set_parent(choice)
    choice.id = "5"

    sequence = Sequence(None, [a, choice])
    a.set_parent(sequence)
    choice.set_parent(sequence)
    sequence.id = "2"

    loop = Loop(None, [e, f])
    e.set_parent(loop)
    f.set_parent(loop)
    loop.id = "3"

    tree = Sequence(None, [sequence, loop])
    sequence.set_parent(tree)
    loop.set_parent(tree)
    tree.id = "1"
    return tree


def run_pipeline(output_dir="./var/quickstart"):
    tree = build_example_tree()

    log = pd.DataFrame({
        'case:concept:name': [1, 2, 2, 2, 3, 3, 3],
        'concept:name': ['b', 'a', 'f', 'e', 'a', 'c', 'e'],
        'time:timestamp': [pd.Timestamp(year=1000 + i, month=1, day=1) for i in range(7)],
    })

    # Model/log path distributions given directly, so this runs without
    # needing Ebi installed. See readme.md for the Ebi-backed alternative
    # (pn_log=/pn_method=), used when deriving these distributions from an
    # actual event log instead of specifying them by hand.
    model_dist = {
        ('4', '8', '6', '7', '6'): 0.1,
        ('4', '9', '6', '7', '6'): 0.1,
        ('4', '8', '6'): 0.3,
        ('4', '9', '6'): 0.3,
    }
    log_dist = {
        ('b',): 0.1,
        ('a', 'f', 'e'): 0.2,
        ('a', 'c', 'e'): 0.7,
    }

    derivation = DerivationPipeline(tree, log, pl=log_dist, pn_measure=model_dist)
    derivation.compute(output_dir)
    return derivation


if __name__ == "__main__":
    d = run_pipeline()
    print(d.print_blinded())
    d.stats()
