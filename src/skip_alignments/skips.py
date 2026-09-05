from typing import Dict
from skip_alignments.alignment import State
from skip_alignments.processtree import *


class Skipper(object):
    def __init__(self):
        pass

    def get_tree_node_by_id(self, tree:ProcessTree, id:str):
        if tree.id == id:
            return tree
        for c in tree.children:
            res = self.get_tree_node_by_id(c, id)
            if res is not None:
                return res
        return None

    def fix_tree_references(self, tree:ProcessTree, agn:State):
        for i in range(len(agn.path)):
            if agn.path[i][1] != '>>':
                if isinstance(agn.path[i][1], Skip):
                    # Skip
                    agn.path[i] = (agn.path[i][0], Skip(self.get_tree_node_by_id(tree, agn.path[i][1].node.id), agn.path[i][1].skip_cost))
                elif isinstance(agn.path[i][1], TauPath):
                    # Tau Skip
                    agn.path[i] = (agn.path[i][0], TauPath(self.get_tree_node_by_id(tree, agn.path[i][1].node.id)))
                else:
                    # Activity
                    agn.path[i] = (agn.path[i][0], self.get_tree_node_by_id(tree, agn.path[i][1].id))

    def fix_sagns(self, tree:ProcessTree, skip_dict:Dict[str, List[State]]):
        for agns in skip_dict.values():
            for agn in agns:
                self.fix_tree_references(tree, agn)

    def count_skip_executions(self, tree:ProcessTree, state:State, number_of_executions:int):
        alignment = state.path
        model_trace = [p[1] for p in alignment if p[1] != '>>']
        counts = 0
        for m in model_trace:
            if (isinstance(m, Skip) or isinstance(m, TauPath)) and m.node == tree:
                counts += 1
        return counts * number_of_executions

    def count_non_skip_executions(self, tree:ProcessTree, state:State, number_of_executions:int):
        alignment = state.path
        model_trace = [p[1] for p in alignment if p[1] != '>>']
        if isinstance(tree, LeafNode):
            return model_trace.count(tree) * number_of_executions
        elif isinstance(tree, Sequence):
            child_executions = []
            for c in tree.children:
                skips = self.count_skip_executions(c, state, number_of_executions)
                non_skips = self.count_non_skip_executions(c, state, number_of_executions)
                child_executions.append(skips+non_skips)
            assert len(child_executions) == child_executions.count(child_executions[0]) # all children executed equally often
            return child_executions[0]
        elif isinstance(tree, Xor):
            child_executions = []
            for c in tree.children:
                skips = self.count_skip_executions(c, state, number_of_executions)
                non_skips = self.count_non_skip_executions(c, state, number_of_executions)
                child_executions.append(skips+non_skips)
            return sum(child_executions)
        elif isinstance(tree, And):
            child_executions = []
            for c in tree.children:
                skips = self.count_skip_executions(c, state, number_of_executions)
                non_skips = self.count_non_skip_executions(c, state, number_of_executions)
                child_executions.append(skips+non_skips)
            assert len(child_executions) == child_executions.count(child_executions[0]) # all children executed equally often
            return child_executions[0]
        elif isinstance(tree, Loop):
            child_executions = []
            for c in tree.children:
                skips = self.count_skip_executions(c, state, number_of_executions)
                non_skips = self.count_non_skip_executions(c, state, number_of_executions)
                child_executions.append(skips+non_skips)
            assert child_executions[0] - sum(child_executions[1:]) >= 0 # no more redo parts than do parts
            return (child_executions[0] - sum(child_executions[1:])) # number of isolated executions
    
    def node_reached(self, node:ProcessTree, state:State) -> bool:
        """
        Whether `node` was reached at all by `state` -- either directly
        (a synchronous move or its own Skip/TauPath), or because it's
        contained within a coarser Skip/TauPath placed on one of its
        ancestors, when the whole ancestor subtree went unwitnessed as one
        block. count_skip_executions/count_non_skip_executions alone only
        catch the former (they use exact-identity matching, which is also
        relied on by Sequence/Xor/And/Loop's own child-count bookkeeping,
        so they can't be made containment-aware themselves without
        double-counting there).

        Containment alone isn't enough, though: attributing a lump skip to
        one specific descendant is only meaningful if something else in
        the same alignment actually synchronized -- otherwise the skip is
        maximally coarse precisely because *nothing* in the trace informs
        which descendant it should be pinned to (e.g. a trace with no
        relation to the model at all skips the whole tree at the root, and
        that doesn't mean every leaf in the tree individually "executed as
        a skip"). So containment only counts as reached when paired with
        at least one genuine synchronous move elsewhere in the alignment.
        """
        skip_cnt = self.count_skip_executions(node, state, 1)
        nskip_cnt = self.count_non_skip_executions(node, state, 1)
        if skip_cnt+nskip_cnt > 0:
            return True
        model_trace = [m for _, m in state.path if m != '>>']
        # id-based, not m.node.contains_tree(node): states computed via the
        # ProcessPoolExecutor in alignall.align_sk_all carry deserialized
        # copies of the tree in their path, not the same objects as `node`
        # (which comes from the caller's own tree) -- same .id, different
        # identity, and ProcessTree has no __eq__ override, so identity-based
        # containment silently fails across that boundary. Comparing by id
        # is what the rest of this codebase already does to bridge it (see
        # get_tree_node_by_id/fix_tree_references above).
        has_containing_skip = any(
            (isinstance(m, Skip) or isinstance(m, TauPath)) and self._contains_id(m.node, node.id)
            for m in model_trace
        )
        if not has_containing_skip:
            return False
        return any(isinstance(m, LeafNode) for m in model_trace)

    def _contains_id(self, container:ProcessTree, target_id:str) -> bool:
        if container.id == target_id:
            return True
        return any(self._contains_id(c, target_id) for c in container.children)

    def _conditional_skip_prob(self, node:ProcessTree, state:State):
        skip_cnt = self.count_skip_executions(node, state, 1)
        nskip_cnt = self.count_non_skip_executions(node, state, 1)
        if skip_cnt+nskip_cnt == 0:
            # ambiguous: either this node's branch genuinely wasn't part of
            # this execution (e.g. the untaken side of an Xor), or it's
            # masked inside a coarser Skip/TauPath on an ancestor. Only the
            # latter should count as skipped.
            return 1 if self.node_reached(node, state) else 0
        return skip_cnt/(skip_cnt+nskip_cnt)
    
    def _traverse_tree(self, tree:ProcessTree):
        if isinstance(tree, LeafNode):
            return [tree]
        nodes = []
        for c in tree.children:
            nodes += self._traverse_tree(c)
        return nodes + [tree]
        

    def conditional_skip_prob(self, tree:ProcessTree, skip_dict:Dict[str, List[State]]):
        # computes P(skip n|sagn)
        # output: sagn_to_node_to_prob: state -> node -> P(skip node | state)
        sagn_to_node_to_prob = {}
        for _, states in skip_dict.items():
            for state in states:
                sagn_to_node_to_prob[state] = {}
                for node in self._traverse_tree(tree):
                    sagn_to_node_to_prob[state][node] = self._conditional_skip_prob(node, state)
        return sagn_to_node_to_prob