import logging
import math
import time
from typing import Dict, Set

from tqdm import tqdm
from skip_alignments.processtree import *
from skip_alignments.alignment import State
from skip_alignments.progress import progress_bars_disabled

logger = logging.getLogger(__name__)
# A separate child logger for progress-bar visibility, so it can be silenced
# independently of the debug instrumentation above -- e.g. debug logging on
# but bars off: logging.getLogger("skip_alignments.execution").setLevel(logging.DEBUG);
# logging.getLogger("skip_alignments.execution.progress").setLevel(logging.WARNING)
progress_logger = logger.getChild("progress")


def _bars_disabled() -> bool:
    """
    True if this module's progress bar(s) should be suppressed -- either
    because the package-wide switch (skip_alignments.progress) says so, which
    always wins, or because progress_logger isn't at DEBUG (the per-module
    default: off unless explicitly opted into).
    """
    return progress_bars_disabled() or not progress_logger.isEnabledFor(logging.DEBUG)

# Above this predicted order-preserving-shuffle count, log a warning (visible
# by default) instead of just a debug line -- a concrete early signal for the
# combinatorial blowup that can otherwise stall for hours with no feedback
# (see the BPI2013 Incidents stall investigated with process-voids).
LARGE_SHUFFLE_COUNT_THRESHOLD = 1_000_000

# Hard stop: above this predicted count, refuse to enumerate at all and
# raise instead of just warning. A BPI2013 Incidents run hit a single
# lengths=[15,11] (predicted_count=7,726,160) interleaving that recurred for
# 24+ minutes with no completion and no way to bound or catch it -- this
# turns that into an immediate, structured, catchable failure instead.
# Above LARGE_SHUFFLE_COUNT_THRESHOLD but below this: still attempted, just
# warned about.
MAX_SHUFFLE_COUNT = 5_000_000


def get_large_shuffle_count_threshold() -> int:
    return LARGE_SHUFFLE_COUNT_THRESHOLD


def set_large_shuffle_count_threshold(n: int) -> None:
    """Above this predicted count, shuffle()'s And-branch logs a warning
    instead of a debug line. Takes effect immediately for every subsequent
    call, at any recursion depth -- this is process-wide configuration, not
    per-call context, so there's nothing to thread through call signatures."""
    global LARGE_SHUFFLE_COUNT_THRESHOLD
    if n <= 0:
        raise ValueError(f"LARGE_SHUFFLE_COUNT_THRESHOLD must be positive, got {n}")
    LARGE_SHUFFLE_COUNT_THRESHOLD = n


def get_max_shuffle_count() -> int:
    return MAX_SHUFFLE_COUNT


def set_max_shuffle_count(n: int) -> None:
    """Above this predicted count, shuffle()'s And-branch raises
    ShuffleExplosionError instead of enumerating. Takes effect immediately
    for every subsequent call, at any recursion depth -- see
    set_large_shuffle_count_threshold."""
    global MAX_SHUFFLE_COUNT
    if n <= 0:
        raise ValueError(f"MAX_SHUFFLE_COUNT must be positive, got {n}")
    MAX_SHUFFLE_COUNT = n


class ShuffleExplosionError(Exception):
    """
    Raised instead of enumerating when an And-branch interleaving's
    predicted output size (see _predicted_shuffle_count) exceeds
    MAX_SHUFFLE_COUNT. Carries the branch lengths and predicted count so a
    caller can log or report which variant/state triggered it, rather than
    the run hanging silently or a Python exception with no context.
    """
    def __init__(self, lengths, predicted_count: int):
        self.lengths = lengths
        self.predicted_count = predicted_count
        super().__init__(
            f"all_order_preserving_shuffles: predicted_count={predicted_count} for "
            f"lengths={lengths} exceeds MAX_SHUFFLE_COUNT={MAX_SHUFFLE_COUNT}"
        )


def _sync_only_merge(paths, sync_rank: Dict[str, int]):
    """
    If every (log, model) element across all paths is a synchronous move --
    i.e. every element's log side is a key in sync_rank, meaning none is a
    model-only '>>' move -- its final relative order is already uniquely
    determined by the log's own order (sync_rank, built from sync_moves_log;
    each label is uniquely position-suffixed by remove_log_moves). In that
    case there is exactly one order-preserving interleaving consistent with
    the log, and it can be computed directly by sorting on rank, with no
    combinatorial search at all.

    Returns [merged_list] (matching all_order_preserving_shuffles' return
    shape: a list of one candidate) when every element is synchronous, or
    None the moment a model-only element is found -- signalling the caller
    to fall back to the general (unconstrained) algorithm, since a
    model-only move's placement genuinely isn't determined by the log.
    """
    merged = []
    for p in paths:
        for e in p:
            if e[0] not in sync_rank:
                return None
            merged.append(e)
    merged.sort(key=lambda e: sync_rank[e[0]])
    return [merged]


def _predicted_shuffle_count(paths) -> int:
    """
    The exact size of all_order_preserving_shuffles(*paths)'s output, without
    generating it: interleaving sequences of lengths l1..lk has
    multinomial(sum(l_i); l1,...,lk) = (sum)! / (l1!*...*lk!) orderings,
    computed here via successive math.comb (choosing each list's positions
    among the slots not yet claimed by earlier lists).
    """
    total = 1
    consumed = 0
    for p in paths:
        total *= math.comb(consumed + len(p), len(p))
        consumed += len(p)
    return total

class ExecutionTree(object):
    parent:Optional["ExecutionTree"] = None
    children:List["ExecutionTree"] = []
    execution:Execution = None

    def __init__(self, parent:Optional["ExecutionTree"], children:List["ExecutionTree"], execution:Execution):
        self.parent = parent
        self.children = children
        self.execution = execution
    
    def set_parent(self, parent:"ExecutionTree"):
        self.parent = parent
    
    def get_distance_to_root(self):
        """
        Returns the level of the current subtree.
        """
        if self.parent is None:
            return 0
        return self.parent.get_distance_to_root() + 1
    
    def __str__(self) -> str:
        child_string = "\n".join([c.__str__() for c in self.children])
        res = (" " * self.get_distance_to_root()*2)
        if isinstance(self.execution.node, LeafNode):
            res += "ACT"
        elif isinstance(self.execution.node, Sequence):
            res += "→"
        elif isinstance(self.execution.node, Xor):
            res += "×"
        elif isinstance(self.execution.node, And):
            res += "∧"
        elif isinstance(self.execution.node, Loop):
            res += "↺"
        else:
            raise ValueError("Execution is covering an unexpected type of node:", self.execution.node)
        res += "(" + self.execution.__str__().replace('\n', " ") + ")"
        res += "\n" + child_string
        return res
    
    def __repr__(self):
        return self.__str__()

    def all_order_preserving_shuffles(self, *paths):
        inter = []
        for comb in itertools.combinations(range(len(paths[0])+len(paths[1])), len(paths[0])):
            inter.append(self.assign(paths[1], list(zip(comb, paths[0]))))
        if len(paths) == 2:
            return inter
        res = []
        for new_path in inter:
            res += self.all_order_preserving_shuffles(new_path, *paths[2:])
        return res

    def assign(self, fillup_list, insert_list):
        l = []
        for i in range(len(insert_list)+len(fillup_list)):
            if len(insert_list) > 0 and insert_list[0][0] == i:
                l.append(insert_list[0][1])
                insert_list = insert_list[1:]
            else:
                l.append(fillup_list[0])
                fillup_list = fillup_list[1:]
        return l
    
    def shuffle(self, state:State, sync_rank: Optional[Dict[str, int]] = None) -> List[List[tuple]]:
        # computes the reshuffling according to the model semantics
        # output: List of lists where the inner list is the new partial alignment
        # sync_rank (optional): label -> position in the log's own true sync
        # order (see ExecutionManager.shuffle, which builds and passes it).
        # Used only by the And case, to skip combinatorial search entirely
        # when every element being interleaved is a synchronous move -- see
        # _sync_only_merge. None preserves the old unconstrained behaviour
        # exactly, for any caller that doesn't pass it.
        if isinstance(self.execution.node, LeafNode):
            # no change, return the move
            return [state.path[self.execution.start:self.execution.stop]]
        elif isinstance(self.execution.node, Sequence) or isinstance(self.execution.node, Loop):
            # no change, fix order
            res_c = []
            for c in self.children:
                res_c.append(c.shuffle(state, sync_rank))
            res = []
            for path in list(itertools.product(*res_c)):
                # path = [ [agn1 for c1], [agn1 for c2], ... ]
                agn = []
                for sub_path in path:
                    agn += sub_path
                res.append(agn)
            return res
        elif isinstance(self.execution.node, Xor):
            # no change, fix selection
            return self.children[0].shuffle(state, sync_rank)
        elif isinstance(self.execution.node, And):
            # changes order
            res_c = []
            for c in self.children:
                res_c.append(c.shuffle(state, sync_rank))
            # res_c is [ [ [agn1 for c1], [agn2 for c1], ... ], [ [agn1 for c2], [agn2 for c2], ... ], ...]
            #            child1 c1                              child2 c2
            res = []
            for path in list(itertools.product(*res_c)):
                # path = [ [agn1 for c1], [agn1 for c2], ... ]
                if sync_rank is not None:
                    fast = _sync_only_merge(path, sync_rank)
                    if fast is not None:
                        res += fast
                        continue
                # cheap (just math.comb), so computed unconditionally -- the
                # warning/error below must stay visible even when debug
                # logging for this module is off, since they're the signal
                # for a blowup that can otherwise stall for hours (or never
                # finish at all) with no feedback
                predicted = _predicted_shuffle_count(path)
                if predicted > MAX_SHUFFLE_COUNT:
                    raise ShuffleExplosionError([len(p) for p in path], predicted)
                if predicted > LARGE_SHUFFLE_COUNT_THRESHOLD:
                    logger.warning("all_order_preserving_shuffles: lengths=%s predicted_count=%d (large)",
                                    [len(p) for p in path], predicted)
                elif logger.isEnabledFor(logging.DEBUG):
                    logger.debug("all_order_preserving_shuffles: lengths=%s predicted_count=%d",
                                 [len(p) for p in path], predicted)
                res += self.all_order_preserving_shuffles(*path)
            return res
        else:
            raise ValueError("Execution is covering an unexpected type of node:", self.execution.node)

class ExecutionManager(object):
    def __init__(self):
        pass

    def correct_to_narrowest_moves(self, state:State):
        # we might need to narrow in the execution indices
        for i in range(len(state.executions)):
            while state.executions[i].stop > state.executions[i].start:
                # log move on start
                if state.path[state.executions[i].start][1] == '>>':
                    state.executions[i].start += 1
                # log move on stop
                elif state.path[state.executions[i].stop-1][1] == '>>':
                    state.executions[i].stop -= 1
                # not contained node, Skip or Taupath
                elif (isinstance(state.path[state.executions[i].start][1], Skip) or isinstance(state.path[state.executions[i].start][1], TauPath)) and not state.executions[i].node.contains_tree(state.path[state.executions[i].start][1].node):
                    state.executions[i].start += 1
                # not contained node, Skip or Taupath
                elif (isinstance(state.path[state.executions[i].stop-1][1], Skip) or isinstance(state.path[state.executions[i].stop-1][1], TauPath)) and not state.executions[i].node.contains_tree(state.path[state.executions[i].stop-1][1].node):
                    state.executions[i].stop -= 1
                # not contained node
                elif not isinstance(state.path[state.executions[i].start][1], Skip) and not isinstance(state.path[state.executions[i].start][1], TauPath) and not state.executions[i].node.contains_tree(state.path[state.executions[i].start][1]):
                    state.executions[i].start += 1
                # not contained node
                elif not isinstance(state.path[state.executions[i].stop-1][1], Skip) and not isinstance(state.path[state.executions[i].stop-1][1], TauPath) and not state.executions[i].node.contains_tree(state.path[state.executions[i].stop-1][1]):
                    state.executions[i].stop -= 1
                else:
                    break
    
    def remove_log_moves(self, state:State):
        # expects correct_to_narrowest_moves before
        i = 0
        log_moves = []
        log_path = []
        while i < len(state.path):
            if state.path[i][1] == '>>':
                log_moves.append(state.path[i][0] + str(len(log_path)))
                log_path.append(state.path[i][0] + str(len(log_path)))
                for exec in state.executions:
                    if exec.start < i and exec.stop > i:
                        # ........ logmove .......
                        # start... stop    ....... -> is okay, as only part before logmove contains execution
                        # start............stop... -> deduce end position
                        exec.stop -= 1
                    elif exec.start > i:
                        exec.start -= 1
                        exec.stop -= 1
                    # i is never a start position after applying correct_to_narrowest_moves
                state.path.pop(i)
            else:
                if state.path[i][0] != '>>':
                    # sync move
                    state.path[i] = (state.path[i][0] + str(len(log_path)), state.path[i][1])
                    log_path.append(state.path[i][0])
                i += 1
        return log_moves, log_path
    
    def build_execution_tree(self, state:State, node:Execution=None) -> ExecutionTree:
        # expects remove_log_moves before
        state.executions = sorted(state.executions, key=lambda x:(x.start, -x.stop, x.node.get_distance_to_root()))

        if node is None:
            node = state.executions[0]
            state.executions.remove(node)
        
        tree = ExecutionTree(None, [], node)

        if isinstance(tree.execution.node, LeafNode):
            return tree
        if isinstance(tree.execution.node, Sequence) or isinstance(tree.execution.node, And):
            # search for the children in bounds of tree.execution
            for c in tree.execution.node.children:
                tree_c = None
                for i in range(len(state.executions)):
                    if state.executions[i].node == c and tree.execution.start <= state.executions[i].start and tree.execution.stop >= state.executions[i].stop:
                        exec_c = state.executions.pop(i)
                        tree_c = self.build_execution_tree(state, exec_c)
                        break
                assert tree_c is not None
                tree.children.append(tree_c)
                tree_c.set_parent(tree)
            return tree
        if isinstance(tree.execution.node, Xor):
            # search for the executed child in bounds of tree.execution
            tree_c = None
            for i in range(len(state.executions)):
                if state.executions[i].node.parent == tree.execution.node and tree.execution.start <= state.executions[i].start and tree.execution.stop >= state.executions[i].stop:
                    exec_c = state.executions.pop(i)
                    tree_c = self.build_execution_tree(state, exec_c)
                    break
            assert tree_c is not None
            tree.children.append(tree_c)
            tree_c.set_parent(tree)
            return tree
        if isinstance(tree.execution.node, Loop):
            # search for all executed children in bounds of tree.execution
            # note that they already appear in appropriate order in the execution
            i = 0
            tree_c = []
            while i < len(state.executions) and state.executions[i].start < tree.execution.stop:
                if state.executions[i].node.parent == tree.execution.node and tree.execution.start <= state.executions[i].start and tree.execution.stop >= state.executions[i].stop:
                    exec_c = state.executions.pop(i)
                    tree_c.append(self.build_execution_tree(state, exec_c))
                else:
                    i += 1
            assert len(tree_c) > 0
            tree.children = tree_c
            for c in tree_c:
                c.set_parent(tree)
            return tree
        raise ValueError("Unexpected type of execution found:", tree.execution)
    
    def shuffle(self, state:State, log_moves:List[str], log_path:List[str], execution_tree:ExecutionTree):
        # state is representing an alignment, not a skip alignment
        # state is not modified but used to copy the data
        # output: generator that yields all the coinciding optimal alignments for an optiaml skip alignment
        sync_moves_log = [l for l in log_path if l not in log_moves]
        # each sync move's log label is uniquely position-suffixed by
        # remove_log_moves, so this rank is exact -- see _sync_only_merge
        sync_rank = {l: i for i, l in enumerate(sync_moves_log)}
        coinciding_agns = execution_tree.shuffle(state, sync_rank)
        generated = len(coinciding_agns)
        passed_sync_filter = 0
        yielded = 0
        for agn in coinciding_agns:
            # agn is a list of pairs, not a state
            if [l for l,_ in agn if l != '>>'] == sync_moves_log:
                passed_sync_filter += 1
                # merge in log moves
                base_agn = [[(l, m) for l,m in agn]]
                new_base = []
                for l in log_moves:
                    for a in base_agn:
                        merge_start = 0
                        merge_end = len(a)
                        log_path_copy = [str(s) for s in log_path]
                        #print(log_moves, log_path)
                        while merge_start < len(a):
                            #print("merge_start =", merge_start, a[merge_start][0])
                            if a[merge_start][0] == '>>' or a[merge_start][0] == log_path_copy[0]:
                                if a[merge_start][0] == log_path_copy[0]:
                                    log_path_copy = log_path_copy[1:]
                                merge_start += 1
                            else:
                                break
                        # go back to the last '>>'
                        while merge_start > 0 and a[merge_start-1][0] == '>>':
                            merge_start -= 1
                        # possibly there is no move with '>>' to sink in, so start at merge_start
                        merge_end = merge_start
                        while merge_end < len(a):
                            if a[merge_end][0] == '>>':
                                merge_end += 1
                            else:
                                break
                        #print("Merging", l, "at", merge_start, merge_end, "into", a)
                        for sink_in in execution_tree.all_order_preserving_shuffles(a[merge_start:merge_end], [(l, '>>')]):
                            new_base.append(a[:merge_start] + sink_in + a[merge_end:])
                    base_agn = [[(m1,m2) for m1,m2 in na] for na in new_base]
                    new_base = []
                for a in base_agn:
                    if [l for l,_ in a if l != '>>'] == log_path:
                        yielded += 1
                        yield a
        if logger.isEnabledFor(logging.DEBUG):
            waste_ratio = 1 - (passed_sync_filter / generated) if generated else 0.0
            logger.debug("shuffle: generated=%d passed_sync_filter=%d yielded=%d waste_ratio=%.3f",
                          generated, passed_sync_filter, yielded, waste_ratio)
        return
    
    def coninciding_agns(self, skip_dict:Dict[str, List[State]]):
        C = {} # state -> set(agns); state can be a skip alignment or an unfolded skip alignment
        global_C = {} # var -> list(set(agns)); list of sets of coinciding agns for the sagns in nf per variant
        ratio_per_var = []
        debug = logger.isEnabledFor(logging.DEBUG)
        for var, states in tqdm(skip_dict.items(), disable=_bars_disabled()):
            var_start = time.process_time_ns()
            already_found = []
            for state in states:
                self.correct_to_narrowest_moves(state)
                log_moves, log_path = self.remove_log_moves(state)
                for s in state.unfold():
                    exec_tree = self.build_execution_tree(s.copy())
                    C_state = self.shuffle(s, log_moves, log_path, exec_tree)
                    C[s] = frozenset([tuple(c) for c in C_state])
                    if state not in C:
                        C[state] = C[s]
                    else:
                        C[state] = frozenset(set(C[state]).union(set(C[s])))
                    already_found.append(C[s])
            if debug:
                logger.debug("coninciding_agns: variant=%r states=%d elapsed_ns=%d",
                              var, len(states), time.process_time_ns() - var_start)
            # a variant with zero states (e.g. a timed-out/degraded-away
            # variant) has no compression ratio to report -- record its
            # (empty) coinciding set but don't let it divide by zero, and
            # don't let it skew the average with a fabricated value
            if states:
                ratio_per_var.append(sum(len(x) for x in already_found)/len(states))
            global_C[var] = already_found
        if ratio_per_var:
            print("Compression in skip alignments was 1 :", sum(ratio_per_var)/len(ratio_per_var))
        else:
            print("Compression in skip alignments was 1 : N/A (no variants with computed alignments)")
        return C, global_C
    
    def validate(self, global_C:Dict[str, Set[List[List|tuple]]]):
        for var, agns in global_C.items():
            for i in range(len(agns)):
                for j in range(i+1, len(agns)):
                    if len(set(agns[i]).intersection(set(agns[j]))) > 0:
                        raise ValueError("Overlap in coinciding agns violating theorem:", var)
        return
    
    def coinciding_agns_var(self, global_C:Dict[str, List[Set[List|tuple]]]):
        var_C = {} # var -> set(agns)
        for k, v in global_C.items():
            var_C[k] = set()
            for fs in v:
                var_C[k] = var_C[k].union(set(fs))
            var_C[k] = list(var_C[k])
        return var_C