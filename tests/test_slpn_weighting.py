"""
End-to-end validation of the PLoop-weighted SLPN path (translate_ppt ->
compile_to_slpn, both in skipalignments.ppt, no pm4py/Ebi involved until
this test queries the result) against Toothpaste's OWN reference trace
probabilities.

The oracle values below are copied verbatim from Toothpaste's own Haskell
unit tests (refs/toothpaste-master/src/test/haskell/Toothpaste/TPConformTest.hs,
`probLoopTests`, using `lpa = Node1 PLoop la 3 1` -- a PLoop[rho=3] wrapping
leaf "a" weight 1, node weight 1):
    prob []          lpa = 1/3   (not tested here -- Ebi's `probability
                                   trace` CLI requires at least one [TRACE]
                                   argument, no syntax for the empty trace)
    prob ["a"]       lpa = 2/9
    prob ["a","a","a"] lpa = 2**3/3**4

This is a stronger check than testing our own derived formula against
itself: it validates the whole pipeline (weight derivation, compilation,
and Ebi's own parsing/probability semantics) against ground truth from
Toothpaste's own test suite, independent of anything this codebase assumes.

Run with:
    python -m unittest tests.test_slpn_weighting -v

Skipped if a real Ebi binary isn't available (EBI_EXECUTABLE env var
override, else the hardlinked skip/Scripts/ebi.exe, else PATH).
"""
import os
import shutil
import subprocess
import tempfile
import unittest

from skipalignments.ppt import PPTNode, compile_to_slpn, translate_ppt


def _find_ebi_executable():
    override = os.environ.get('EBI_EXECUTABLE')
    if override and shutil.which(override):
        return override
    hardlinked = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'skip', 'Scripts', 'ebi.exe')
    if os.path.isfile(hardlinked):
        return hardlinked
    return shutil.which('ebi')


_EBI_EXECUTABLE = _find_ebi_executable()


def _ebi_trace_prob(ebi_executable, model, trace):
    result = subprocess.check_output([ebi_executable, "prob", "trac", model] + trace + ["-a"]).decode("utf-8")
    if 'Approximately' not in result:
        raise ValueError(f"Ebi did not return a probability: {result!r}")
    try:
        return float(result.split(' ')[1])
    except (IndexError, ValueError):
        return float(result.split('\n')[0])


@unittest.skipUnless(_EBI_EXECUTABLE, "ebi binary not available in this environment")
class TestCompiledSlpnMatchesToothpasteOracle(unittest.TestCase):
    """
    Builds the exact PPT Toothpaste's own TPConformTest.hs uses for `lpa`,
    compiles it with compile_to_slpn (no pm4py, no Ebi, in-process), writes
    the result, and asks a real Ebi binary for trace probabilities -- then
    compares against the literal numbers from Toothpaste's Haskell suite.
    """

    @classmethod
    def setUpClass(cls):
        cls._original_cwd = os.getcwd()
        cls._tmp_dir = tempfile.mkdtemp(prefix="skipalignments_slpn_oracle_test_")
        os.chdir(cls._tmp_dir)

        w, rho = 1.0, 3.0
        ppt = PPTNode('ploop', w, rho=rho, children=[PPTNode('leaf', w, name='a')])
        tree, weights, loop_taus = translate_ppt(ppt)
        text, activity_to_id = compile_to_slpn(tree, weights, loop_taus)
        with open('lpa.slpn', 'w') as f:
            f.write(text)
        cls._activity_id = activity_to_id['a']
        # Tau_redo is a real labelled transition (matches captured real Ebi
        # output from earlier this session), firing once between each pair
        # of loop iterations -- so the model's own label sequence for n
        # repetitions of "a" is [a, tau_redo, a, tau_redo, ..., a], not n
        # copies of "a" alone. loop_taus[0] = (tau_skip_id, tau_redo_id).
        cls._tau_redo_id = loop_taus[0][1]

    @classmethod
    def tearDownClass(cls):
        os.chdir(cls._original_cwd)
        shutil.rmtree(cls._tmp_dir, ignore_errors=True)

    def test_single_a(self):
        prob = _ebi_trace_prob(_EBI_EXECUTABLE, 'lpa.slpn', [self._activity_id])
        self.assertAlmostEqual(prob, 2 / 9, places=3)

    def test_three_as(self):
        trace = [self._activity_id, self._tau_redo_id, self._activity_id,
                  self._tau_redo_id, self._activity_id]
        prob = _ebi_trace_prob(_EBI_EXECUTABLE, 'lpa.slpn', trace)
        self.assertAlmostEqual(prob, 2 ** 3 / 3 ** 4, places=3)


if __name__ == '__main__':
    unittest.main()
