# TODO

Tracked improvements and technical debt.

## Ebi binary used instead of ebi Python package

A fairly broad Python wrapper is now available for most ebi features. The 
current requirement for ebi to be on the PATH or to override 
`probabilities.EBI_EXECUTABLE` can be obsoleted if this package will suffice.

This would also fix the CLI-shape fragility this code currently has: the
`disc occ`/`prob trac` subprocess calls in `probabilities.py` are written
against a specific Ebi CLI version's argument shape (already had to be
patched once — see CHANGELOG — after Ebi's `disc occ` gained a required
`stochastic-labelled-Petri-net` subcommand token), with no version pin or
check anywhere. A Python dependency would put Ebi under normal dependency
management instead.


## Add LICENSE and clarify LICENSE chain from upstream, particularly pm4py

As well as using as a library dependency, alignall.py is from pm4py.

`pm4py/algo/conformance/alignments/petri_net/variants/state_equation_a_star.py`


## `EbiOccurance` hardcodes output paths relative to cwd

`write_tree_to_petri` and `write_log` in `src/skipalignments/probabilities.py`
write `model.pnml` and `log.xes` as literal relative paths instead of
accepting them as parameters (unlike `ebi_slpn`, `validate_slpn`,
`update_slpn_weights`, and `ebi_trace_prob`, which already take a
`model`/`path` argument). Any caller — including tests — gets these files
dumped into whatever directory happened to be the process's cwd at the
time.

Currently worked around in `tests/test_run_example.py` and
`tests/test_external_usage.py` by running the pipeline inside a temp
directory (`tempfile.mkdtemp()` + `os.chdir()`, in `setUpClass`/
`tearDownClass` and a `try`/`finally` respectively) rather than fixing
the library.

Fix: give `write_tree_to_petri`/`write_log` the same `model`/`log`
parameters (with the current literals as defaults) that the other
`EbiOccurance` methods already have, and update `derivation.py`'s call
sites accordingly.

