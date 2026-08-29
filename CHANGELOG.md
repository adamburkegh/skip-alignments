# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

Cleanup pass bringing this package back up to parity with `process-voids`'s
independently-patched copy of the same engine code, and fixing the
regression introduced when it was originally extracted from that monorepo.
Thanks to Joshua Gong for the heavy lifting here.

### Fixed
- Restored the `EBI_EXECUTABLE` / `MISSING_ACTIVITY_WEIGHT` configuration
  knobs in `probabilities.py`, removing a hardcoded `ebi.exe` regression
  from the original extraction.
- Merged `process-voids`'s independent fixes back into this package: the
  zero-weight-activity fix when forcing SLPN weights (`update_visible_taus`
  renamed back to `update_slpn_weights`), more robust Ebi output parsing in
  `ebi_trace_prob`, and the `out=` kwarg on `ebi_slpn`.
- Fixed the package's internal imports: sibling modules now use real
  absolute imports (`from skipalignments.processtree import *`, etc.)
  instead of bare names that only worked via a `sys.path` hack in
  `__init__.py`. This was the root cause of the `NameError` that had
  stalled `process-voids` from depending on this package via `pip`.
- Rewrote `tests/test_run_example.py` and `tests/test_external_usage.py`
  against the current public API — both still referenced the pre-rename
  `skipprobabilities` package and a `run_pipeline` / `build_example_tree`
  API that no longer existed. Ported to the standard-library `unittest`
  rather than `pytest`, which the suite had assumed without ever declaring
  it as a dependency anywhere in the project. Run with `python -m
  unittest` from the project root.
- Added `tests/__init__.py` so `python -m unittest` discovers the suite
  from the project root; without it, bare discovery silently found zero
  tests instead of recursing into `tests/`.
- Test runs no longer leak `model.pnml`, `log.xes`, or output directories
  into the repository root; the pipeline now runs inside a temp directory
  for the duration of each test.
- `ebi_slpn()`'s `disc occ` call was missing a required subcommand token
  (`stochastic-labelled-Petri-net`) between `occ` and the log/model
  arguments — the current Ebi CLI rejects the file path in that slot with
  `Usage: ebi discover occurrence <COMMAND>`. Found via real-data testing
  in the parallel `process-voids` session (local Ebi build: 0.3.14).

### Added
- `__init__.py` is now a full facade over the package: `from skipalignments
  import *` exposes process trees, alignment, execution, probabilities,
  skips, and derivation, plus `update_pair_taus`, `check_names`,
  `get_variant_dict`, `get_activities`, and `generate_tree` — restoring the
  intent of the old `skipalignments.py` aggregator module.
- `probabilities` is exported as a submodule so callers can override its
  configuration directly, e.g.
  `skipalignments.probabilities.EBI_EXECUTABLE = "..."`.
- `examples/quickstart.py`, a plain-Python (no Jupyter) equivalent of
  `examples/example.ipynb`, runnable directly with `python
  examples/quickstart.py`. Writes its output to `var/quickstart/`
  (gitignored) rather than the repository root.

### Removed
- `main.py`, which referenced a non-existent `lib/` and
  `external_tests/run_example.py`.
- `src/skipalignments/skipalignments.py`, folded into `__init__.py`.

## [0.1.1] - 2026-08-17

### Changed
- Renamed the package to `skipalignments`.
- Regenerated `requirements.txt` from `pip list` after running
  `example.ipynb`.

## [0.1.0] - 2026-08-06

Initial extraction of the skip-probability engine into an installable
library, undoing the copy-paste fork relationship with `process-voids`.

### Added
- `requirements.txt`, and an initial `runexample.py` ported from the
  original Jupyter notebook, verified to match `example.ipynb`'s output.
- Restructured the project into a standard `src`-layout Python package
  installable via `pip install -e .` (`pyproject.toml`, setuptools
  backend).
- `main.py` demonstrating the package used as a library.
- A `pytest` suite covering example tree construction, end-to-end pipeline
  execution, skip-probability coverage and range checks, pickled output
  artifacts, and the reporting methods.

### Changed
- Refactored `run_example.py` into `build_example_tree()` and
  `run_pipeline(output_dir=...)` so the example could be invoked
  programmatically rather than only via script.
- Replaced assertion-based tests with user-defined test cases to better
  represent expected behaviour.
- Reorganised the repository: separated example notebooks, results, and
  package source into `examples/`, `results/`, and `src/`.
- Resolved `run_example.py` import/output issues caused by incorrect
  project structure (outputs had been relying on hard-coded values).
- Updated example notebooks to import via `from skipprobabilities import *`
  instead of directly from the original implementation, and verified
  outputs still matched the original notebooks.
