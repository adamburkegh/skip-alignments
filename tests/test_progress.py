"""
Unit tests for skipalignments.progress: the central, explicit on/off switch
for every tqdm progress bar in the package. Added after a lab run showed
progress bars still appearing even with per-module logger configuration --
disable_progress_bars() is meant to be a single, reliable way to silence
all of them regardless of any logger's level, called explicitly by the
caller rather than picked up implicitly from the environment.

Run with:
    python -m unittest tests.test_progress -v
"""
import unittest

import skipalignments.progress as progress


class TestProgressToggle(unittest.TestCase):

    def setUp(self):
        self._orig_disabled = progress.progress_bars_disabled()

    def tearDown(self):
        # module-global state -- must not leak between tests
        if self._orig_disabled:
            progress.disable_progress_bars()
        else:
            progress.enable_progress_bars()

    def test_enabled_by_default(self):
        progress.enable_progress_bars()
        self.assertFalse(progress.progress_bars_disabled())

    def test_disable_and_enable_round_trip(self):
        progress.disable_progress_bars()
        self.assertTrue(progress.progress_bars_disabled())
        progress.enable_progress_bars()
        self.assertFalse(progress.progress_bars_disabled())

    def test_disable_is_idempotent(self):
        progress.disable_progress_bars()
        progress.disable_progress_bars()
        self.assertTrue(progress.progress_bars_disabled())


if __name__ == '__main__':
    unittest.main()
