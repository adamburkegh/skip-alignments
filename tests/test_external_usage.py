"""
Confirms the installed package works when imported the way a downstream
consumer (e.g. process-voids) would use it.

Run with:
    pip install -e .
    python -m unittest tests.test_external_usage -v
"""
import os
import shutil
import tempfile
import unittest

import pandas as pd

from skip_alignments import DerivationPipeline, Activity, Sequence, Xor


class ExternalUsageTests(unittest.TestCase):

    def test_import_works_from_outside_the_repo(self):
        self.assertIsNotNone(DerivationPipeline)
        self.assertIsNotNone(Activity)

    def test_pipeline_runs_from_outside_the_repo(self):
        # EbiOccurance writes model.pnml/log.xes relative to cwd, so run
        # from a throwaway temp directory instead of leaking them into the
        # repo.
        original_cwd = os.getcwd()
        tmp_dir = tempfile.mkdtemp(prefix="skip_alignments_test_")
        os.chdir(tmp_dir)
        try:
            a = Activity(None, 'a', 100000)
            a.id = "2"
            b = Activity(None, 'b', 100000)
            b.id = "3"
            tree = Sequence(None, [a, b])
            a.set_parent(tree)
            b.set_parent(tree)
            tree.id = "1"

            log = pd.DataFrame({
                'case:concept:name': [1],
                'concept:name': ['a'],
                'time:timestamp': [pd.Timestamp(year=2000, month=1, day=1)],
            })
            log_dist = {('a',): 1.0}
            model_dist = {('2', '3'): 1.0}

            derivation = DerivationPipeline(tree, log, pl=log_dist, pn_measure=model_dist)
            derivation.compute("./external_test_out")
            self.assertIsInstance(derivation, DerivationPipeline)
            self.assertGreater(len(derivation.skip_probs), 0)
        finally:
            os.chdir(original_cwd)
            shutil.rmtree(tmp_dir, ignore_errors=True)


if __name__ == "__main__":
    unittest.main()
