from typing import Dict
from alignment import State
from processtree import *


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
    
    def _conditional_skip_prob(self, node:ProcessTree, state:State):
        skip_cnt = self.count_skip_executions(node, state, 1)
        nskip_cnt = self.count_non_skip_executions(node, state, 1)
        if skip_cnt+nskip_cnt == 0:
            # we did never come to the subtree
            return 0
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