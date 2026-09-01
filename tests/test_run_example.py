"""
Test suite for the skipalignments package's example pipeline.

Run with:
    python -m unittest tests.test_run_example -v

Assumes the package has been installed (e.g. `pip install -e .` from the
project root, per pyproject.toml).
"""
import contextlib
import io
import os
import pickle
import shutil
import tempfile
import unittest

import pandas as pd

from skipalignments import DerivationPipeline, Activity, Sequence, Xor, Loop, LeafNode


OUTPUT_DIR = "./test_example_out"


def build_example_tree():
    """The running-example process tree from examples/example.ipynb."""
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


def build_example_log():
    return pd.DataFrame({
        'case:concept:name': [1, 2, 2, 2, 3, 3, 3],
        'concept:name': ['b', 'a', 'f', 'e', 'a', 'c', 'e'],
        'time:timestamp': [pd.Timestamp(year=1000 + i, month=1, day=1) for i in range(7)],
    })


def run_pipeline(output_dir):
    tree = build_example_tree()
    log = build_example_log()
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


class TreeConstructionTests(unittest.TestCase):

    def test_build_example_tree_structure(self):
        tree = build_example_tree()
        self.assertIsInstance(tree, Sequence)
        self.assertEqual(tree.id, "1")
        self.assertEqual(len(tree.children), 2)

        seq1, loop = tree.children
        self.assertIsInstance(seq1, Sequence)
        self.assertIsInstance(loop, Loop)

        a, choice = seq1.children
        self.assertIsInstance(a, Activity)
        self.assertEqual(a.name, "a")
        self.assertIsInstance(choice, Xor)

        c, d = choice.children
        self.assertEqual({c.name, d.name}, {"c", "d"})

        e, f = loop.children
        self.assertEqual({e.name, f.name}, {"e", "f"})

    def test_tree_leaf_labels(self):
        tree = build_example_tree()
        self.assertEqual(set(tree.get_leaf_labels()), {"a", "c", "d", "e", "f"})


class PipelineTests(unittest.TestCase):
    """
    Runs the full example pipeline once for the class and reuses the result
    across tests, inside a throwaway temp directory: EbiOccurance writes
    model.pnml/log.xes relative to cwd, so this keeps those (and
    OUTPUT_DIR) out of the repo instead of leaking them into the working
    tree.
    """

    @classmethod
    def setUpClass(cls):
        cls._original_cwd = os.getcwd()
        cls._tmp_dir = tempfile.mkdtemp(prefix="skipalignments_test_")
        os.chdir(cls._tmp_dir)
        cls.derivation = run_pipeline(output_dir=OUTPUT_DIR)

    @classmethod
    def tearDownClass(cls):
        os.chdir(cls._original_cwd)
        shutil.rmtree(cls._tmp_dir, ignore_errors=True)

    def test_pipeline_runs_without_error(self):
        self.assertIsInstance(self.derivation, DerivationPipeline)

    def test_pipeline_produces_skip_probs(self):
        self.assertTrue(hasattr(self.derivation, "skip_probs"))
        self.assertIsInstance(self.derivation.skip_probs, dict)
        self.assertGreater(len(self.derivation.skip_probs), 0)

    def test_skip_probs_cover_every_tree_node(self):
        def all_nodes(node):
            nodes = [node]
            for c in getattr(node, "children", []):
                nodes += all_nodes(c)
            return nodes

        nodes = all_nodes(self.derivation.tree)
        for node in nodes:
            self.assertIn(node, self.derivation.skip_probs, f"Missing skip prob for node {node.id}")

    def test_skip_probs_are_valid_probabilities(self):
        # small float tolerance for accumulated summation error, not a
        # logic bound -- probabilities are genuinely in [0, 1]
        eps = 1e-9
        for node, prob in self.derivation.skip_probs.items():
            self.assertTrue(-eps <= prob <= 1.0 + eps, f"Node {node.id} has out-of-range prob {prob}")

    def test_leaf_skip_probs_are_floats(self):
        for node, prob in self.derivation.skip_probs.items():
            if isinstance(node, LeafNode):
                self.assertIsInstance(prob, float)

    def test_skip_prob_matches_paper_worked_example(self):
        # "Skip Probabilities for Subprocesses" (ICPM 2025), running
        # example: P_{N_5}(sk) = 0.22 for the choice node between the two
        # ways the appeal subprocess can go (here: 'c'/'d'), derived from
        # sigma_2 and sigma_3 each having their own optimal skip alignments
        # in normal form and coinciding alignments over this node.
        choice = self.derivation.tree.children[0].children[1]
        self.assertIsInstance(choice, Xor)
        self.assertAlmostEqual(self.derivation.skip_probs[choice], 0.22, places=2)

    def test_cond_prob_matches_paper_sigma2_skip_alignment_split(self):
        # Same paper, same example: sigma_2 has two optimal skip alignments
        # in normal form, sagn_21 (synchronizes 'f') and sagn_22 (performs
        # a log move on 'f' instead), with P_{sigma_2}(sagn_21) = 0.25 and
        # P_{sigma_2}(sagn_22) = 0.75. This checks the intermediate
        # conditional skip alignment probability, one layer below the final
        # skip_probs aggregate checked above.
        _, choice = self.derivation.tree.children[0].children
        states = self.derivation.skip_dict["a, f, e"]
        self.assertEqual(len(states), 2)
        seen_synced_f = False
        for state in states:
            synced_f = any(getattr(m, 'name', None) == 'f' for _, m in state.path)
            prob = self.derivation.cond_prob_[choice][state]
            if synced_f:
                seen_synced_f = True
                self.assertAlmostEqual(prob, 0.25, places=2)  # sagn_21
            else:
                self.assertAlmostEqual(prob, 0.75, places=2)  # sagn_22
        self.assertTrue(seen_synced_f)

    def test_output_files_written(self):
        for filename in ["tree", "skip_dict", "trace_probs", "trace_counts", "skip_probs"]:
            with self.subTest(filename=filename):
                path = os.path.join(OUTPUT_DIR, filename)
                self.assertTrue(os.path.exists(path), f"Expected pickled output '{filename}' not found")
                with open(path, "rb") as fh:
                    obj = pickle.load(fh)
                self.assertIsNotNone(obj)

    def test_pickled_skip_probs_matches_in_memory(self):
        path = os.path.join(OUTPUT_DIR, "skip_probs")
        with open(path, "rb") as fh:
            pickled = pickle.load(fh)
        self.assertEqual(len(pickled), len(self.derivation.skip_probs))

    def test_stats_runs_without_error(self):
        captured = io.StringIO()
        with contextlib.redirect_stdout(captured):
            self.derivation.stats()
        output = captured.getvalue()
        self.assertIn("Skip alignment computation", output)
        self.assertIn("Derivation skip probabilities", output)

    def test_print_blinded_runs_without_error(self):
        output = self.derivation.print_blinded()
        self.assertIsInstance(output, str)
        self.assertGreater(len(output), 0)


if __name__ == "__main__":
    unittest.main()
