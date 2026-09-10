"""
Unit test for skipalignments.__version__, added on a process-voids
courtesy request: `import skipalignments; skipalignments.__version__`
raised AttributeError. Not load-bearing for process-voids itself (their
own run-log version banner already goes via importlib.metadata directly),
but the conventional attribute consumers reach for interactively.

Sourced from importlib.metadata.version("skipalignments") rather than a
hand-maintained string literal, so it can't drift from pyproject.toml's
own version the way a literal eventually would when someone bumps the
release and forgets the module.

Run with:
    python -m unittest tests.test_package_version -v
"""
import unittest
from importlib.metadata import version

import skipalignments


class TestPackageVersion(unittest.TestCase):

    def test_version_attribute_matches_installed_distribution_version(self):
        self.assertEqual(skipalignments.__version__, version("skipalignments"))


if __name__ == '__main__':
    unittest.main()
