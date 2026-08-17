# src/skipalignments/__init__.py
import os
import sys

_pkg_dir = os.path.dirname(os.path.abspath(__file__))
if _pkg_dir not in sys.path:
    sys.path.insert(0, _pkg_dir)


from derivation import DerivationPipeline, EbiWeights
from processtree import ProcessTree, LeafNode, Activity, Tau, Sequence, Xor, And, Loop
from alignment import Aligner, State, Mapper

__all__ = [
    "DerivationPipeline", "EbiWeights",
    "ProcessTree", "LeafNode", "Activity", "Tau", "Sequence", "Xor", "And", "Loop",
    "Aligner", "State", "Mapper",
]