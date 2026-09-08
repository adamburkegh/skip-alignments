# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).


## [Unreleased]

### Added
- `alignall.TieExplosionError`, with a configurable ceiling
  (`get_max_tie_count`/`set_max_tie_count`, default 100,000) and an
  earlier warning threshold (`get_large_tie_count_threshold`/
  `set_large_tie_count_threshold`, default 10,000): `align_pn_all`/
  `align_pn_all_multi`'s classical A* search raises instead of continuing
  to enumerate once accumulated tied optimal alignments for a single
  variant exceed the ceiling. Mirrors `execution.py`'s
  `ShuffleExplosionError`/`MAX_SHUFFLE_COUNT`/`LARGE_SHUFFLE_COUNT_THRESHOLD`
  for the same underlying phenomenon on the other alignment path: a
  classical Petri net has no notion of "these two moves are concurrent,
  their order doesn't matter" the way the process-tree-native path's `And`
  node does, so the search counts every interleaving of concurrent
  branches as a genuinely distinct optimal alignment. Found via a
  process-voids report on rtfm.xes: a 2-activity variant against a
  heavily `And`/`Xor(Tau,_)`-nested discovered tree produced 4115+ tied
  optimal alignments well before the per-trace timeout, with no feedback
  and no way to bound or catch it — confirmed to be this and not the loop
  cycle-guard above (`id_loop_list` was empty for that tree; no `Loop`
  nodes at all). As with `MAX_SHUFFLE_COUNT`, there is no single correct
  default ceiling — it trades off how many ties an experiment actually
  needs (e.g. for a `|optimal alignments|` statistic) against how long a
  single variant may run; tune via the setter.

### Fixed
- `align_pn_all`/`align_pn_all_multi`'s `__search` scaled far worse than
  the number of tied optimal alignments should require, for two stacking
  reasons, both pure implementation inefficiencies with zero effect on
  output: `__reconstruct_alignment` rebuilt its alignment-so-far list via
  repeated `[parent.t] + alignment` while walking from a state up to the
  search root — O(depth) work per step, for O(depth) steps, so O(depth²)
  for a single reconstruction, called at least once per popped state (and
  more, since `all_optimal=True` keeps searching after the first
  solution); and `closed` (the visited-state set) was a plain list, so
  `is_closed`'s `search_tuple in closed` was an O(n) linear scan repeated
  n times. Together these dominated real run time once tie counts reached
  the thousands — exactly the process-voids rtfm.xes case above. Fixed by
  building the alignment list via `append` while walking up and reversing
  once at the end (O(depth) total), and by making `closed` a `set` of
  hashable `(tuple(alignment), marking)` pairs (pm4py's `Transition` and
  `Marking` are both hashable — id-based and frozenset-based respectively
  — so this is an identical equality check, just O(1) average instead of
  O(n)). Measured on a synthetic fixture (n independent concurrent
  optional branches, n! tied optimal alignments): n=7 (5040 ties) went
  from 43.0s to 1.4s, with an *identical* tie count before and after —
  see `tests/test_alignall_search_performance.py` (fast regression test)
  and `tests/perf_search_scaling.py` (full scaling numbers up to n=8,
  40320 ties, too slow for the normal suite).

### Changed
- Renamed `Aligner.align2` to `Aligner.align_normal_form`, to distinguish it
  from classical alignment (see `alignall.align_pn_all`/`align_pn_all_multi`)
  and to name what it actually does: compute optimal skip alignments in
  normal form, lumping an entirely-unwitnessed subtree into one Skip/TauPath
  move rather than one model move per missing leaf activity. `align2` is
  kept as a backward-compatible alias (`PendingDeprecationWarning`, not
  removed) for collaborators still depending on an ancestor project's API.
- `EbiOccurance.build_petri_net` now returns a 6-tuple
  (`net, im, fm, activity_to_id, tau_ids, id_loop_list`), adding `tau_ids`
  (the set of transition labels that are genuine `Tau` leaves) and
  `id_loop_list` (loops needing a cycle guard, from `insert_cycle_checks` —
  see below). **Breaking** for any direct caller of `build_petri_net` (not
  `write_tree_to_petri`, whose own return contract is unchanged).

### Fixed
- `align_pn_all`/`align_pn_all_multi`/`align_pn_all_for_one`/
  `align_pn_one_for_one`'s `net_model_move` priced a genuine `ProcessTree`
  `Tau` leaf's model move at the full activity cost (100000) instead of 0,
  whenever the net came from `EbiOccurance.build_petri_net` (which labels
  every leaf, `Tau` included, with its own opaque tree id via
  `to_pm4py(use_ids=True)` — the `"TAU"`-prefix check only matches
  `to_pm4py(use_ids=False)`'s labelling). This mispriced the model
  legitimately choosing its own free tau branch (e.g. `Xor(activity, Tau)`)
  as a deviation, which could bias the A* search away from a genuinely
  optimal alignment when comparing against equal-cost alternatives, not
  just mislabel a correct one after the fact. Found via a process-voids
  report on the running-example fixture (`payment_approval`'s
  `schedule_choice = Xor(s, Tau)`). Fixed by having `build_petri_net`
  additionally return `tau_ids`, and the four functions above take an
  optional `tau_ids` parameter to price those transitions at 0; omitting
  it preserves the old (mispriced) behaviour for hand-built nets that use
  `to_pm4py(use_ids=False)`'s `"TAU_"`-prefix convention instead.
- `align_pn_all`/`align_pn_all_multi` (classical alignment) never actually
  bounded the A* search on a loop that can execute entirely for free (both
  its do- and some redo-child allow a zero-cost path) — the existing
  `insert_cycle_checks`/`id_loop_list` cycle guard was never wired into
  `build_petri_net`, so `id_loop_list` was always empty. Even passing one
  in by hand wouldn't have helped: `does_allow_tau_path`
  (`ProcessTree.from_pm4py`) and `is_cycling_exec`'s do/redo tau-tracking
  both string-sniff a `"TAU_"` label prefix, the same
  `to_pm4py(use_ids=False)`-only convention behind the cost bug above — so
  on an ids-labelled net (`build_petri_net`'s only output), no loop was
  ever recognized as needing a guard, and no tau'd execution inside one
  was ever recognized as such during search. The result: a trace whose
  optimal alignment goes through such a loop has infinitely many
  equal-cost ties, and the search just burns the full per-trace timeout
  instead of cutting the tie off. Found via a process-voids report on a
  real 104k-case log (rtfm.xes): 42 of 48 trace variants clustered at the
  ~100s per-variant timeout regardless of variant size or weight — ~74 of
  a 76-minute run spent inside `align_pn_all`. Fixed by threading the same
  `tau_ids` through `ProcessTree.from_pm4py`/`does_allow_tau_path`/
  `insert_cycle_checks` (so the right loops get a guard inserted) and
  through `is_cycling`/`is_cycling_exec` (so a tau'd do/redo execution is
  recognized during search), and by having `build_petri_net` call
  `insert_cycle_checks` and return its `id_loop_list`. Also fixed a latent,
  previously-harmless bug this surfaced in `build_petri_net`'s label
  substitution loop: every transition label with no match in the original
  tree (pm4py's own invisible transitions, and now `insert_cycle_checks`'
  `"TAU_entry_"`/`"TAU_exit_"` sentinels too) was aliased through the same
  `None` dict key, silently collapsing distinct sentinel labels into
  whichever one was seen first — harmless while `None` was the only such
  label, not once there are others.

## [0.2.1] - 2026-09-05

### Changed
- Attempted to rename the PyPI distribution and importable package from
  `skipalignments` to `skip-alignments`/`skip_alignments`, to match the git
  repository's own name. Reverted: PyPI's anti-typosquatting check blocks
  registering `skip-alignments` as "too similar" to the (now-deleted)
  existing `skipalignments` project, and resolving that requires manual
  PyPI support intervention on a multi-week queue. The distribution and
  import name remain `skipalignments`, unchanged from `0.2.0`.

### Fixed
- Added defensiveness for divide by zero in 
-- `ExecutionManager.coninciding_agns` divided by zero (`ZeroDivisionError`)
  computing the "Compression in skip alignments" ratio whenever there were
  no variants to report on at all (an empty `skip_dict`, e.g. a
  degraded-to-nothing log) or a variant with zero computed states. Found
  via a real process-voids run. A zero-states variant's (empty) coinciding
  set is still recorded; it's just excluded from the ratio average rather
  than dividing by zero or skewing it with a fabricated value.
-- `DerivationPipeline.compute()`, same degraded-to-nothing root cause
-- `DerivationPipeline.stats()`, average sagns per variant, average time per
  variant among non-timed-out ones, and average agns per variant all
  divided by a count that's zero for a degraded-to-nothing result.

## [0.2.0] - 2026-09-05

First PyPI release.

### Added
- `skipalignments.ppt`: direct import of a Toothpaste Miner Probabilistic
  Process Tree (`.ptree` export) as a fixed-weight model, decoupled from
  log-driven occurrence weighting. `translate_ppt` maps PPT's
  PLoop/FLoop/Seq/Choice/Conc onto skip-alignments' own process-tree
  operators, including the loop-topology translation PPT's reusable
  loop-place and skip-alignments' 1-or-more redo-loop require; `compile_to_slpn`
  then renders the translated tree directly to Ebi's plaintext `.slpn`
  format, with every transition's weight attached at the moment it's
  created — no pm4py, no PNML file, and no Ebi subprocess call anywhere in
  the path. Validated against Toothpaste's own Haskell reference trace
  probabilities (`TPConformTest.hs`) via a real Ebi binary, not just this
  codebase's own derivation. See `skip_align_improvements.md`.
- `skipalignments.execution.ShuffleExplosionError`, with a configurable
  `MAX_SHUFFLE_COUNT` ceiling (`get_max_shuffle_count`/`set_max_shuffle_count`):
  a hard, catchable stop for the And-node interleaving search in
  `ExecutionTree.shuffle`, which is combinatorial in concurrent-branch
  length and could previously stall for tens of minutes with no feedback
  and no way to bound or catch it. Found via a real process-voids run
  against BPI2013 Incidents.
- A fast path for And-node interleaving (`_sync_only_merge`): when every
  element being interleaved is a synchronous log move, its relative order
  is already uniquely determined by the log's own order, so the
  combinatorial search is skipped in favour of a single sort.
- stdlib `logging` instrumentation in `execution.py` (predicted shuffle
  counts, generate/keep waste ratios, per-variant timing), silent by
  default — enable with
  `logging.getLogger("skipalignments.execution").setLevel(logging.DEBUG)`.
- `skipalignments.progress`: an explicit, package-wide
  `disable_progress_bars()`/`enable_progress_bars()` switch covering every
  `tqdm` progress bar in the library, since per-module logger configuration
  alone could leave bars reappearing unexpectedly (a child logger with no
  level of its own inherits its parent's).
- `DerivationPipeline` now has a real constructor path for the PPT/Toothpaste
  weight source: `DiscoverySource.TOOTHPASTE` plus a new `pn_ppt_weights`
  constructor argument (the `(weights, loop_taus)` pair from `translate_ppt`).
  `compute()`'s TOOTHPASTE branch calls `write_slpn` directly, with no
  `EbiOccurance.write_tree_to_petri`/pm4py/PNML step. Tested end to end,
  including a real Ebi query against the compiled `.slpn`.
- `LICENSE` (AGPL-3.0-only), declared in `pyproject.toml` (`license`/`license-files`)
  and bundled into both built distributions. Required by this project's
  dependence on, and adapted code from, pm4py (AGPL-3.0) — see the README's
  "Third Party Dependencies and licenses" section.
- `pyproject.toml` now carries full PyPI listing metadata: `readme`, `authors`,
  `keywords`, `classifiers`, and project `urls`. Added a `dev` extra
  (`build`, `twine`) for building/validating releases.

### Changed
- Renamed `EbiWeights` to `DiscoverySource` (`OCCURANCE`/`UNIFORM`, now also
  `TOOTHPASTE`). The old name was misleading the same way `EbiOccurance` is:
  Ebi is only ever the final query backend, the same for every source — the
  enum selects where model weights come from, not anything about Ebi
  specifically. **Breaking**: any code importing `EbiWeights` directly needs
  to update to `DiscoverySource`.

### Fixed
- `translate_ppt` used the same `model_move_cost` for both `Activity` and
  `Tau` nodes, violating the alignment engine's own invariant
  (`Aligner.align2`: `assert tau_cost < activity_cost`) — every real
  alignment run against a translated PPT tree failed that assertion. Found
  via the end-to-end `DerivationPipeline` test above, not by any of
  `ppt.py`'s own isolated tests. Fixed by giving `Tau` nodes their own
  `model_move_tau_cost` (default `0`, matching the codebase-wide convention
  used everywhere else, e.g. `ProcessTree.from_pm4py`'s callers).

## [0.1.2] - 2026-09-01

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
