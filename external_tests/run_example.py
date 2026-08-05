"""
Standalone version of example.ipynb — run with:
    python run_example.py

Builds the process tree:  ->( ->(a, x(c,d)), loop(e,f) )
aligns it against an example event log and derives skip probabilities.
"""
import pandas as pd

from skipprobabilities import *

#Build the process tree 
a = Activity(None, 'a', 100000); a.id = "4"
c = Activity(None, 'c', 100000); c.id = "8"
d = Activity(None, 'd', 100000); d.id = "9"
e = Activity(None, 'e', 100000); e.id = "6"
f = Activity(None, 'f', 100000); f.id = "7"

choice = Xor(None, [c, d])
c.set_parent(choice); d.set_parent(choice); choice.id = "5"

sequence = Sequence(None, [a, choice])
a.set_parent(sequence); choice.set_parent(sequence); sequence.id = "2"

loop = Loop(None, [e, f])
e.set_parent(loop); f.set_parent(loop); loop.id = "3"

tree = Sequence(None, [sequence, loop])
sequence.set_parent(tree); loop.set_parent(tree); tree.id = "1"

#event log
log = pd.DataFrame({
    'case:concept:name': [1, 2, 2, 2, 3, 3, 3],
    'concept:name': ['b', 'a', 'f', 'e', 'a', 'c', 'e'],
    'time:timestamp': [pd.Timestamp(year=1000 + i, month=1, day=1) for i in range(7)],
})

#Probability distributions
#model_dist: probability per model path (keyed by leaf ids along that path)
model_dist = {
    ('4', '8', '6', '7', '6'): 0.1,
    ('4', '9', '6', '7', '6'): 0.1,
    ('4', '8', '6'): 0.3,
    ('4', '9', '6'): 0.3,
}
# log_dist: probability per trace variant
log_dist = {
    ('b',): 0.1,
    ('a', 'f', 'e'): 0.2,
    ('a', 'c', 'e'): 0.7,
}


derivation = DerivationPipeline(tree, log, pl=log_dist, pn_measure=model_dist)
derivation.compute("./example_out")

print()
print(derivation.print_blinded())
print()
derivation.stats()