"""
Regression coverage for DerivationPipeline.compute() writing model.pnml/
log.xes into whatever directory the caller's process happened to be
running from, rather than the pipeline's own caller-nominated output
directory (`path`) -- flagged by a process-voids report: fixed filenames
in cwd can collide across concurrent runs (their sweeps run 13 workers),
the same class of bug as the earlier smodel.slpn test leak in this repo.

write_tree_to_petri/write_log already accepted optional model/log_path
parameters (see CHANGELOG.md, 0.2.2) -- compute() itself just never
threaded them through, so the capability existed but wasn't used by the
pipeline's own default path.

Run with:
    python -m unittest tests.test_derivation_output_file_paths -v
"""
import os
import shutil
import tempfile
import unittest

from tests.test_run_example import build_example_tree, build_example_log
from skipalignments import DerivationPipeline


class TestComputeWritesModelPnmlUnderTheGivenOutputPath(unittest.TestCase):

    def setUp(self):
        self._original_cwd = os.getcwd()
        self._tmp_dir = tempfile.mkdtemp(prefix="skipalignments_test_")
        os.chdir(self._tmp_dir)

    def tearDown(self):
        os.chdir(self._original_cwd)
        shutil.rmtree(self._tmp_dir, ignore_errors=True)

    def test_model_pnml_lands_under_output_path_not_bare_cwd(self):
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
        output_dir = "example_out"
        derivation.compute(output_dir)

        self.assertTrue(
            os.path.exists(os.path.join(output_dir, "model.pnml")),
            "expected model.pnml under the pipeline's own output directory",
        )
        self.assertFalse(
            os.path.exists("model.pnml"),
            "model.pnml leaked into the bare process cwd instead of the output directory",
        )

    def test_explicit_model_path_override_is_still_honoured(self):
        tree = build_example_tree()
        log = build_example_log()
        model_dist = {
            ('4', '8', '6', '7', '6'): 0.1,
            ('4', '9', '6', '7', '6'): 0.1,
            ('4', '8', '6'): 0.3,
            ('4', '9', '6'): 0.3,
        }
        derivation = DerivationPipeline(tree, log, pn_measure=model_dist)
        custom_dir = "custom_models"
        os.mkdir(custom_dir)
        derivation.compute("example_out", model_path=os.path.join(custom_dir, "mymodel.pnml"))

        self.assertTrue(os.path.exists(os.path.join(custom_dir, "mymodel.pnml")))
        self.assertFalse(os.path.exists(os.path.join("example_out", "model.pnml")))


if __name__ == '__main__':
    unittest.main()
