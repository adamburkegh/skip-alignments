"""
Unit tests for DerivationPipeline's Toothpaste/PPT weight source
(DiscoverySource.TOOTHPASTE) -- wiring translate_ppt/compile_to_slpn (see
skip_alignments/ppt.py and skip_align_improvements.md) into the pipeline's
public constructor, the piece flagged as "not yet wired up" after the PPT
import work landed. EbiWeights was renamed to DiscoverySource in the same
change, since the enum selects a weight-*source*, not anything Ebi-specific
-- OCCURANCE and TOOTHPASTE both ultimately query Ebi at the same final
step, but neither the name nor (for TOOTHPASTE) any step before that has
anything to do with Ebi.

Run with:
    python -m unittest tests.test_derivation_toothpaste -v

The end-to-end test is skipped if a real Ebi binary isn't available (see
tests/test_slpn_weighting.py for the same convention).
"""
import os
import shutil
import tempfile
import unittest

import pandas as pd

import skip_alignments.probabilities as probabilities
from skip_alignments.derivation import DerivationPipeline, DiscoverySource
from skip_alignments.ppt import PPTNode, translate_ppt


def _find_ebi_executable():
    override = os.environ.get('EBI_EXECUTABLE')
    if override and shutil.which(override):
        return override
    hardlinked = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'skip', 'Scripts', 'ebi.exe')
    if os.path.isfile(hardlinked):
        return hardlinked
    return shutil.which('ebi')


_EBI_EXECUTABLE = _find_ebi_executable()


class TestDiscoverySourceRename(unittest.TestCase):

    def test_members_exist(self):
        self.assertTrue(hasattr(DiscoverySource, 'OCCURANCE'))
        self.assertTrue(hasattr(DiscoverySource, 'UNIFORM'))
        self.assertTrue(hasattr(DiscoverySource, 'TOOTHPASTE'))

    def test_importable_from_package_facade(self):
        import skip_alignments
        self.assertIs(skip_alignments.DiscoverySource, DiscoverySource)


class TestDerivationPipelineToothpasteConstruction(unittest.TestCase):
    """__init__ validation for the TOOTHPASTE source -- doesn't need Ebi."""

    def _tree_and_log(self):
        ppt = PPTNode('leaf', 1.0, name='a')
        tree, weights, loop_taus = translate_ppt(ppt)
        log = pd.DataFrame({
            'case:concept:name': [1],
            'concept:name': ['a'],
            'time:timestamp': [pd.Timestamp(year=2000, month=1, day=1)],
        })
        return tree, log, weights, loop_taus

    def test_toothpaste_requires_pn_ppt_weights(self):
        tree, log, weights, loop_taus = self._tree_and_log()
        with self.assertRaises(AssertionError):
            DerivationPipeline(tree, log, pn_method=DiscoverySource.TOOTHPASTE)

    def test_toothpaste_with_pn_ppt_weights_constructs(self):
        tree, log, weights, loop_taus = self._tree_and_log()
        derivation = DerivationPipeline(
            tree, log, pn_method=DiscoverySource.TOOTHPASTE, pn_ppt_weights=(weights, loop_taus),
        )
        self.assertEqual(derivation.pn_method, DiscoverySource.TOOTHPASTE)
        self.assertEqual(derivation.pn_ppt_weights, (weights, loop_taus))
        self.assertIsNone(derivation.pn_log)


@unittest.skipUnless(_EBI_EXECUTABLE, "ebi binary not available in this environment")
class TestDerivationPipelineToothpasteEndToEnd(unittest.TestCase):
    """
    Real end-to-end run: a translated PLoop tree, a small log with variants
    at 0/1/2 iterations, computed through DerivationPipeline.compute() using
    DiscoverySource.TOOTHPASTE -- no EbiOccurance estimation step anywhere
    in this path (see skip_align_improvements.md), just compile_to_slpn's
    weighted .slpn plus the generic Ebi trace_probs query.
    """

    @classmethod
    def setUpClass(cls):
        cls._original_cwd = os.getcwd()
        cls._tmp_dir = tempfile.mkdtemp(prefix="skip_alignments_toothpaste_derivation_test_")
        os.chdir(cls._tmp_dir)
        probabilities.EBI_EXECUTABLE = _EBI_EXECUTABLE

    @classmethod
    def tearDownClass(cls):
        os.chdir(cls._original_cwd)
        shutil.rmtree(cls._tmp_dir, ignore_errors=True)

    def test_pipeline_runs_and_produces_valid_skip_probs(self):
        w, rho = 1.0, 3.0
        ppt = PPTNode('ploop', w, rho=rho, children=[PPTNode('leaf', w, name='a')])
        tree, weights, loop_taus = translate_ppt(ppt)

        # variants at 1, 2, and 3 occurrences of 'a' -- exercises the loop
        log = pd.DataFrame({
            'case:concept:name': [1, 2, 2, 3, 3, 3],
            'concept:name':      ['a', 'a', 'a', 'a', 'a', 'a'],
            'time:timestamp': [pd.Timestamp(year=2000 + i, month=1, day=1) for i in range(6)],
        })

        derivation = DerivationPipeline(
            tree, log, pn_method=DiscoverySource.TOOTHPASTE, pn_ppt_weights=(weights, loop_taus),
        )
        derivation.compute(self._tmp_dir)

        skip_probs = derivation.results()
        self.assertEqual(set(skip_probs.keys()), _all_nodes(tree))
        for node, prob in skip_probs.items():
            self.assertGreaterEqual(prob, 0.0)
            self.assertLessEqual(prob, 1.0)


def _all_nodes(tree):
    nodes = {tree}
    for c in tree.children:
        nodes |= _all_nodes(c)
    return nodes


if __name__ == '__main__':
    unittest.main()
