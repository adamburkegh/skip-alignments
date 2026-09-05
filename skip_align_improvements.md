# skip-alignments improvements

Working notes on improvements/additions to skip-alignments. One section per
topic.

## PPT translation: importing Toothpaste Miner (PPT) models into skip-alignments

### Motivation

An experiment holds a discovered model fixed while degrading the log across
runs (missing information simulation). The model comes from Toothpaste
Miner, which performs *direct* stochastic discovery — it discovers its own
control-flow structure, not just weights for an existing one — and outputs a
Probabilistic Process Tree (PPT) with real per-transition weights, including
loop continuation/exit probabilities derived from a geometric-distribution
parameter ρ.

`DerivationPipeline` currently only supports two weight sources: log
occurrence counting via Ebi (`EbiWeights.OCCURANCE`, re-derives weights from
`pn_log` every run) and a fully precomputed path-probability dict
(`pn_measure`). Neither fits: re-deriving from `pn_log` would let the
model's weights drift with each run's log degradation, defeating the point
of holding the model fixed; `pn_measure` requires the caller to have already
done all trace-probability math themselves before ever touching
`DerivationPipeline`.

This section covers the third path: importing a PPT (structure + weights)
directly as the model, decoupled from whatever log is aligned against it.

### Why not use Toothpaste's own SLPN/PNML export directly

Toothpaste's own tooling exports a weighted PNML alongside the PPT. The
naive plan — hand that PNML to Ebi as `trace_probs`'s `model=` argument
unmodified — doesn't work, because it isn't just an ID-labelling mismatch.

`DerivationPipeline`'s alignment computation produces `model_path_ids` —
sequences of *this codebase's own* tree-node ids — and queries Ebi with
them. For that query to mean anything against a given `.slpn`, that file's
transition labels have to be exactly skip-alignments' own id scheme. That's
achievable for *visible* activities (rename to match). It is not achievable
for the *loop* construct, because the two systems encode "zero or more
repetitions" with genuinely different topologies:

- PPT's `\loopp{\rho}(x{:}w)` is a single reusable "loop place" with two
  outgoing edges every iteration (exit / continue-into-child-then-loop-back).
- skip-alignments' `Loop` node is the classic process-tree redo-loop
  (`do`; `(redo; do)*`), which is **1-or-more** by construction — confirmed
  both structurally (`alignment.py`'s `Loop` handling forces `children[0]`
  before any redo branch) and empirically (`Skip.unfold()`'s `Loop` branch
  in `processtree.py` always substitutes exactly one execution of
  `children[0]`, never zero).

So even after correct ID renaming, skip-alignments' own alignment machinery
can produce a `model_path` (zero occurrences of the loop's activity) that
provably cannot exist in PPT's SLPN's language, and vice versa. This isn't
fixable by renaming — the shapes disagree. The fix is to reconcile the
*structures*, not just the labels, and only then attach weights.

### The chosen approach

1. Import the PPT into skip-alignments' own `ProcessTree` classes
   (`Sequence`, `Xor`, `And`, `Loop`, `Activity`, `Tau`), translating loops
   per the scheme below, minting fresh unique ids as usual
   (`skipalignments/ppt.py`'s `translate_ppt`).
2. Compile the translated tree directly to Ebi's `.slpn` plaintext format
   (`compile_to_slpn`), with every transition's weight attached at the
   moment it's created — no pm4py, no PNML file, no Ebi subprocess call
   anywhere in this path, and no `EbiOccurance` involvement (that class is
   log-occurrence *estimation*, and Toothpaste already gives us exact
   weights — there's nothing to estimate). See "Implementation status"
   below for why this replaced an earlier pm4py-based design.
3. `ebi.trace_probs(var_C, model=<compiled slpn>)`, exactly as today.

### Loop translation

`\loopp{\rho}(x{:}w)` (single child, zero-or-more, geometric with exit
probability `1/ρ` per iteration) translates to:

```
Xor(Tau_skip, Loop(x, Tau_redo))
```

`Tau_skip` covers the zero-iteration case (skip-alignments' `Loop` can't
represent that on its own — see above). `Tau_redo` is the loop's redo child,
making each additional iteration silent from the model's own perspective
(only `x` is visible).

Compiling this fragment directly (`compile_to_slpn`'s block-structured
`Xor`/`Loop` cases) produces **5** transitions total: 3 labelled
(`Tau_skip`, `Tau_redo`, `x`) and 2 unlabelled structural ones — the `Xor`'s
own gate into the `Loop` branch (entering the loop the first time, distinct
from `Tau_redo`, which only fires when looping back from inside — merging
them would wrongly re-enable `Tau_skip` after every iteration) and the
loop's own exit transition. No separate "harmless wrapping" transitions
exist in this direct compiler (unlike the earlier pm4py-based design,
where a leaf `Xor` branch got a redundant extra gate transition it didn't
need).

### Weight derivation

Three of the four weighted transitions follow directly from PPT's own
stated conservation rules (Xor: children's weights sum to the parent's;
Loop: a child inherits the loop's own weight), applied to the translated
fragment as if it inherited the original PPT loop node's weight `w`:

- **Xor rule**, split to reproduce PPT's own exit/continue probabilities
  (`1/ρ`, `(ρ-1)/ρ`): `weight(Tau_skip) = w/ρ`; the `Loop` branch itself
  gets `w(ρ-1)/ρ` (recorded under the `Loop` node's own id, so
  `compile_to_slpn`'s `Xor` case can look up any branch's selection weight
  uniformly, leaf or compound).
- **Loop rule**: both of `Loop`'s children inherit the loop's own weight:
  `weight(x) = weight(Tau_redo) = w(ρ-1)/ρ`.

The fourth — the unlabelled "stop after this iteration" transition the
compiler inserts — has no PPT-native counterpart (PPT encodes this choice
once, on a single reusable place; skip-alignments' redo-loop encodes it as a
recurring decision inside a `do` cycle instead), so there's no conservation
rule to appeal to directly. It's pinned instead by memorylessness: the
geometric distribution requires *every* realization of "exit or continue"
to carry the same `1/ρ` : `(ρ-1)/ρ` ratio, or the iteration-count
distribution stops being geometric with parameter `ρ`. Since it competes
against `Tau_redo` (already pinned at `w(ρ-1)/ρ`) at the same place:

```
weight(stop-transition) = w/ρ
```

Summary, for a PPT loop node with weight `w` and parameter `ρ`:

| Transition                          | Weight        |
|--------------------------------------|---------------|
| `Tau_skip` (Xor branch)              | `w/ρ`         |
| loop-entry (unlabelled, structural)  | `w(ρ-1)/ρ`    |
| `x` (the loop's visible child)       | `w(ρ-1)/ρ`    |
| `Tau_redo`                           | `w(ρ-1)/ρ`    |
| loop-exit (unlabelled, structural)   | `w/ρ`         |

Non-loop PPT operators (`Sequence`, `Xor`, `And`, plain leaves) map onto
skip-alignments' equivalents structurally as-is; their weights are PPT's own
stated conservation values, carried straight across without translation.

### Toothpaste's actual export format (confirmed against source)

Toothpaste's own source (`refs/toothpaste-master`, downloaded from GitHub)
settles this rather than leaving it assumed. It exports three files per run
(`Main.hs`, `--pnetfile`/`--ptreefile`/`--traceprobfile`): a weighted PNML
(flat Petri net, no tree structure — see `ProcessFormats.hs`), a `.ptree`
file (the actual PPT, tree-shaped — what this importer targets), and
optionally trace probabilities. Only the `.ptree` file is ever read — the
PNML export is never touched by this importer, for two independent reasons:
Ebi's own PNML importer doesn't read stochastic weights from PNML at all
(confirmed against Ebi's own Rust source — every PNML-consuming Ebi command
assigns its own weights and ignores any `toolspecific` extension), and even
if it did, the PNML's topology is PPT's own loop-place shape, not
skip-alignments' redo-loop shape — see "Why not use Toothpaste's own
SLPN/PNML export directly" above.

The `.ptree` format (`ProbProcessTree.hs`, `formatPPTreeIndent`), confirmed
against real sample output (`results/2021_pn/teleclaims_k1.ptree`):

```
<indent>Seq:704.0
<indent>  "incoming claim":704.0
<indent>  PLoop[2.0]:344.0
<indent>    "some activity":344.0
<indent>  tau:47.0
<indent>  Conc:93.0
<indent>    ...
```

- Indentation is 2 spaces per level.
- A leaf line is `"<name>":<weight>`.
- A silent/tau leaf is `tau:<weight>`.
- `Node1` (unary) is `PLoop[<r>]:<weight>` or `FLoop[<r>]:<weight>`, followed
  by exactly one child line at indent+1. `PLoop` is the probabilistic loop;
  `FLoop` is the fixed-count loop.
- `NodeN` is `Seq:<weight>` / `Choice:<weight>` / `Conc:<weight>`, followed
  by 2+ child lines at indent+1.

Confirmed against the actual weight-derivation code too, not just the
paper: `TPConform.hs`'s `probPLoop`/`pathset` implementations use `1/r` for
the exit probability and `(r-1)/r` for the continuation scaling factor —
i.e. `r` **is** `ρ` exactly as used above, and the `w/ρ` / `w(ρ-1)/ρ`
formulas derived from the paper match Toothpaste's own implementation
directly, not just its description in the thesis chapter.

`FLoop[r]` (fixed loop) is simpler and needs no probability split: per
PPT's own stated rule, it's `Sequence`d `round(r)` copies of the child, each
inheriting the loop's own weight unchanged — no ρ, no synthetic taus.

### Implementation status

Implemented and tested, including an end-to-end check against a real Ebi
binary (0.3.14) and, independently, against Toothpaste's own Haskell
reference trace probabilities (`TPConformTest.hs`'s `probLoopTests`) — not
just this codebase's own derivation:

- `skipalignments/ppt.py`: `.ptree` parsing (`parse_ptree`), the PLoop
  translation and weight derivation (`translate_ppt`, recording a weight
  for every translated node, not just leaves, so `compile_to_slpn`'s `Xor`
  case can treat any branch — leaf or compound — uniformly), and the direct
  block-structured `ProcessTree`-to-`.slpn` compiler (`compile_to_slpn`,
  `write_slpn`).
- `tests/test_ppt.py`, `tests/test_slpn_weighting.py`: unit tests of the
  translation/compilation, plus an integration test comparing real Ebi
  output against Toothpaste's own `lpa = Node1 PLoop la 3 1` oracle values
  (`prob ["a"] lpa = 2/9`, `prob ["a","a","a"] lpa = 2**3/3**4`) — exact
  matches, validating the whole pipeline (weight derivation, compilation,
  and Ebi's own parsing) against ground truth independent of anything this
  codebase assumes.

**This replaced an earlier pm4py-based design** (build a Petri net via
pm4py, generate a uniform-weight-1 `.slpn` skeleton via `ebi discover
uniform`, then patch specific transitions' weights in by locating them
topologically) that turned out to be the wrong shape entirely, for reasons
surfaced during review:

1. `EbiOccurance` (where the pm4py path lived) is log-occurrence
   *estimation* machinery — Toothpaste already gives exact weights, so nesting
   this work inside an estimator class was a category error, not just a
   naming one.
2. pm4py's tree-to-net compiler is weight-oblivious, forcing weights to be
   re-attached *after* structure was built, by topologically finding
   pm4py's own anonymous helper transitions (`find_loop_structural_transitions`)
   — indirect, and only as robust as that topological search.
3. Compiling `ProcessTree` directly to `.slpn` text needs no external
   library at all: every operator (`Sequence`/`Xor`/`And`/`Loop`/`Activity`/
   `Tau`) has a standard block-structured Petri-net compilation, so a
   transition can be minted together with its weight in the same pass —
   structurally impossible to end up with an unweighted transition, rather
   than something verified after the fact.

`build_petri_net`/`find_loop_structural_transitions`/`build_pnml_id_to_weight`/
`discover_uniform_slpn`/`write_weighted_slpn` (the pm4py-based versions) no
longer exist; `compile_to_slpn` replaced all of them.

### Wired into DerivationPipeline

`DiscoverySource.TOOTHPASTE` (the enum was renamed from `EbiWeights` --
see "Naming: EbiWeights -> DiscoverySource" below) selects this path in
`compute()`: construct with

```python
tree, weights, loop_taus = translate_ppt(parse_ptree(open("model.ptree").read()))
derivation = DerivationPipeline(tree, aligned_log,
                                 pn_method=DiscoverySource.TOOTHPASTE,
                                 pn_ppt_weights=(weights, loop_taus))
derivation.compute(path)
```

`compute()`'s TOOTHPASTE branch calls `write_slpn` directly (no
`EbiOccurance.write_tree_to_petri`, no pm4py, no PNML) and then the same
generic `ebi.trace_probs`/`ebi.skip_agn_probs_traversal` every other source
already uses. Tested end to end in
`tests/test_derivation_toothpaste.py`, including a real
`ebi.trace_probs` query against a compiled `.slpn`.

A real bug surfaced by that end-to-end test, not by any of `ppt.py`'s own
isolated tests: `translate_ppt` used the same `model_move_cost` for both
`Activity` and `Tau` nodes, violating the alignment engine's own invariant
(`alignment.py`'s `Aligner.align2`: `assert tau_cost < activity_cost`) --
every real alignment run against a translated PPT tree failed that
assertion. Fixed by giving `Tau` nodes their own `model_move_tau_cost`
(default `0`, matching the codebase-wide convention used everywhere else,
e.g. `ProcessTree.from_pm4py`'s callers: activity `100000`, tau `0`),
threaded through `translate_ppt`/`_translate` as a new parameter alongside
`model_move_cost`.

### Naming: EbiWeights -> DiscoverySource

The weight-source enum (`OCCURANCE`/`UNIFORM`, now also `TOOTHPASTE`) was
named `EbiWeights`, which is misleading the same way `EbiOccurance` was
(see the PPT translation section above): Ebi is only ever the final query
backend for `trace_probs`, the same for every source -- OCCURANCE and
TOOTHPASTE both end up querying Ebi identically, but neither the enum's
name nor (for TOOTHPASTE) anything before that final step has anything to
do with Ebi specifically. Renamed to `DiscoverySource`, matching the
domain language this document already uses ("Toothpaste performs *direct*
stochastic discovery") -- occurrence-counting and PPT/Toothpaste both are
genuinely discovery methods in that sense; uniform is the trivial
baseline. `EbiWeights` no longer exists; consumers importing it directly
(rather than via `skipalignments.DerivationPipeline`'s own API) need to
update to `DiscoverySource`.

## Alignment computation performance: the And-node interleaving blowup

Surfaced using skip-alignments as a library from process-voids' dose-response
harness (`exp_disco_degrade`), which found a real combinatorial cost in the
existing `coninciding_agns`/`ExecutionTree.shuffle` machinery in
`execution.py`.

### The problem

`ExecutionTree.shuffle`'s `And` branch enumerates every order-preserving
interleaving of `k` concurrent branches via `all_order_preserving_shuffles`,
then `ExecutionManager.shuffle` filters that set down to the ones whose log
side matches the trace's actual observed order. For branches of lengths
`l1..lk` that's `multinomial(sum(l_i); l1,...,lk)` candidates generated, most
of them thrown away.

A real run against BPI2013 Incidents (mined via the inductive miner) stalled
in this exact code path. Standard-library `logging` instrumentation (see
below) pinned it down precisely:

- A single `lengths=[15,11]` interleaving (`predicted_count=7,726,160`)
  recurred for 24+ minutes with no completion and no way to bound or catch
  it.
- A `lengths=[13,9]` interleaving (`predicted_count=497,420`), called ~90
  times for the same variant, generated **41.78 million** candidates to keep
  **659** -- a `waste_ratio` of 1.000, i.e. essentially all generated work
  was discarded.

Both stalls occurred at the *undegraded* baseline (dose level 0.0), so
neither is related to log degradation, and neither is addressed by caching
alignment results across dose-response levels -- this is a stage-2
(`coninciding_agns`) cost, independent of stage 1
(`compute_skipalignments`).

### Logging instrumentation

`execution.py` uses stdlib `logging` (`logging.getLogger(__name__)`), silent
by default:

- `all_order_preserving_shuffles` calls in the `And` branch log their
  predicted output size (computed analytically via `_predicted_shuffle_count`,
  no enumeration needed) -- `debug` normally, `warning` (visible without
  opting in) above `LARGE_SHUFFLE_COUNT_THRESHOLD`.
- `ExecutionManager.shuffle` logs `generated`/`passed_sync_filter`/`yielded`
  counts per call, at `debug` -- this is what surfaced the 41.78M/659 waste
  ratio.
- `coninciding_agns` logs per-variant elapsed time at `debug`, and its
  `tqdm` progress bar has its own child logger
  (`skipalignments.execution.progress`) so it can be silenced independently
  of debug logging -- useful in an automated sweep where a wall of bars
  per dose-response point is noise, not signal.

Enable with `logging.getLogger("skipalignments.execution").setLevel(logging.DEBUG)`.

### Two mitigations

**A hard ceiling.** Above `MAX_SHUFFLE_COUNT`, `shuffle()`'s `And` branch
raises `ShuffleExplosionError` (carrying `lengths`/`predicted_count`)
instead of ever calling `all_order_preserving_shuffles`. Turns the
7,726,160 case into an immediate, structured, catchable failure instead of
an unbounded stall.

**Skip the search when there's nothing to search for.** A *synchronous*
move's log label is uniquely position-suffixed by `remove_log_moves` (e.g.
`"Accepted3"`), so its relative order across all branches is already fixed
by the log itself -- the only genuine freedom is where *model-only* (`'>>'`)
moves land, since they have no log anchor. `_sync_only_merge(paths,
sync_rank)` returns the single correct interleaving directly (a sort, not a
search) when every element across all branches is synchronous, falling back
to the unchanged `all_order_preserving_shuffles` the instant a model-only
move is found:

```python
def _sync_only_merge(paths, sync_rank: Dict[str, int]):
    merged = []
    for p in paths:
        for e in p:
            if e[0] not in sync_rank:
                return None          # a model-only move -- fall back
            merged.append(e)
    merged.sort(key=lambda e: sync_rank[e[0]])
    return [merged]
```

`ExecutionTree.shuffle` takes an optional `sync_rank: Dict[str, int] = None`
parameter (threaded through recursion, since it's per-call context derived
from the specific state/log being processed -- unlike the thresholds below,
it can't live as shared config). `ExecutionManager.shuffle` builds it from
the `sync_moves_log` it already computes and passes it in; `None` preserves
the old unconstrained behaviour exactly for any caller that doesn't opt in.

Real branches in the BPI log are dominated by synchronous moves (repeated
`Accepted`/`Queued`/etc.), so this collapses most of the observed waste to
O(1) per call; the ceiling catches whatever remains genuinely explosive.

### Configuring the thresholds

`LARGE_SHUFFLE_COUNT_THRESHOLD`/`MAX_SHUFFLE_COUNT` are process-wide
configuration, not per-call context -- `shuffle()`'s `And` branch already
reads them as bare module-global names, so a setter mutates that same
global and takes effect for every subsequent call at any recursion depth,
with nothing threaded through any call signature (contrast with
`sync_rank` above, which genuinely has to be threaded):

```python
from skipalignments.execution import set_max_shuffle_count, set_large_shuffle_count_threshold

set_large_shuffle_count_threshold(500_000)  # warn earlier
set_max_shuffle_count(2_000_000)            # abort earlier
```

`get_max_shuffle_count()`/`get_large_shuffle_count_threshold()` read the
current values; both setters reject non-positive input.

### Tests

`tests/test_execution.py`: `TestSyncOnlyMerge` and `TestConstrainedAndShuffle`
prove the fast path on/off equivalence -- not just that the fast answer is
*a* valid candidate, but that it equals exactly what filtering the
unconstrained output the same way `ExecutionManager.shuffle` does would
keep, and `TestExecutionManagerShuffleWithFastPath` proves it through the
real call site (not just `ExecutionTree.shuffle` directly).
`TestShuffleExplosionCeiling` and `TestConfigurableThresholds` cover the
ceiling and its config API, including that the fast path bypasses the
ceiling entirely (it never computes a candidate count at all) and that
lowering `MAX_SHUFFLE_COUNT` takes effect with no threading.

### Not yet done

Process-voids observed the *same* stall independently, back to back, for
both the `activity` and `trace` degradation dimensions at dose level 0.0 --
structurally identical computations (no degradation applied means the same
tree and log either way). That's a harness-side deduplication opportunity
(reuse the 0.0-level result across dimensions), not a library change.
