"""
Test suite for the skipprob package's example pipeline.

Run with:
    pytest test_run_example.py -v

Assumes the package has been installed (e.g. `pip install -e .` from the
project root, per the pyproject.toml setup). If you haven't packaged it
yet and are still using the flat-folder layout, add the folder to
sys.path before running pytest (see the `conftest.py` note at the bottom
of this file).
"""
import os
import pickle
import shutil

import pytest

from skipprobabilities import run_pipeline, build_example_tree, DerivationPipeline
from skipprobabilities import Activity, Sequence, Xor, Loop, ProcessTree, LeafNode


OUTPUT_DIR = "./test_example_out"


@pytest.fixture(scope="module")
def derivation():
    """Runs the full example pipeline once and reuses the result across tests."""
    if os.path.exists(OUTPUT_DIR):
        shutil.rmtree(OUTPUT_DIR)
    d = run_pipeline(output_dir=OUTPUT_DIR)
    yield d
    shutil.rmtree(OUTPUT_DIR, ignore_errors=True)



# Tree construction
def test_build_example_tree_structure():
    tree = build_example_tree()
    assert isinstance(tree, Sequence)
    assert tree.id == "1"
    assert len(tree.children) == 2

    seq1, loop = tree.children
    assert isinstance(seq1, Sequence)
    assert isinstance(loop, Loop)

    a, choice = seq1.children
    assert isinstance(a, Activity) and a.name == "a"
    assert isinstance(choice, Xor)

    c, d = choice.children
    assert {c.name, d.name} == {"c", "d"}

    e, f = loop.children
    assert {e.name, f.name} == {"e", "f"}


def test_tree_leaf_labels():
    tree = build_example_tree()
    assert set(tree.get_leaf_labels()) == {"a", "c", "d", "e", "f"}



# Pipeline execution
def test_pipeline_runs_without_error(derivation):
    assert isinstance(derivation, DerivationPipeline)


def test_pipeline_produces_skip_probs(derivation):
    assert hasattr(derivation, "skip_probs")
    assert isinstance(derivation.skip_probs, dict)
    assert len(derivation.skip_probs) > 0


def test_skip_probs_cover_every_tree_node(derivation):
    def all_nodes(node):
        nodes = [node]
        for c in getattr(node, "children", []):
            nodes += all_nodes(c)
        return nodes

    nodes = all_nodes(derivation.tree)
    for node in nodes:
        assert node in derivation.skip_probs, f"Missing skip prob for node {node.id}"


def test_skip_probs_are_valid_probabilities(derivation):
    for node, prob in derivation.skip_probs.items():
        assert 0.0 <= prob <= 1.0, f"Node {node.id} has out-of-range prob {prob}"


def test_leaf_skip_probs_are_floats(derivation):
    for node, prob in derivation.skip_probs.items():
        if isinstance(node, LeafNode):
            assert isinstance(prob, float)


# Output artifacts written to disk
@pytest.mark.parametrize("filename", [
    "tree",
    "skip_dict",
    "trace_probs",
    "trace_counts",
    "skip_probs",
])
def test_output_files_written(derivation, filename):
    path = os.path.join(OUTPUT_DIR, filename)
    assert os.path.exists(path), f"Expected pickled output '{filename}' not found"
    with open(path, "rb") as fh:
        obj = pickle.load(fh)
    assert obj is not None


def test_pickled_skip_probs_matches_in_memory(derivation):
    path = os.path.join(OUTPUT_DIR, "skip_probs")
    with open(path, "rb") as fh:
        pickled = pickle.load(fh)
    assert len(pickled) == len(derivation.skip_probs)


# Sanity checks on reported stats/timings
def test_stats_runs_without_error(derivation, capsys):
    derivation.stats()
    captured = capsys.readouterr()
    assert "Skip alignment computation" in captured.out
    assert "Derivation skip probabilities" in captured.out


def test_print_blinded_runs_without_error(derivation):
    output = derivation.print_blinded()
    assert isinstance(output, str)
    assert len(output) > 0