"""
Import Toothpaste Miner's .ptree export format (Probabilistic Process
Trees) into skip-alignments' own ProcessTree classes, with SLPN weights
derived from the PPT's own w/rho values.

Format and weight-derivation formulas confirmed against Toothpaste's own
Haskell source (ProbProcessTree.hs's formatPPTreeIndent for the .ptree
grammar, TPConform.hs's probPLoop/pathset for the weight formulas) -- see
ppt_translation.md for the full derivation.
"""
import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

from skipalignments.processtree import Activity, And, LeafNode, Loop, ProcessTree, Sequence, Tau, Xor

DEFAULT_MODEL_MOVE_COST = 100000


@dataclass
class PPTNode:
    kind: str  # 'leaf', 'tau', 'ploop', 'floop', 'seq', 'choice', 'conc'
    weight: float
    name: Optional[str] = None
    rho: Optional[float] = None  # the 'r' value, for ploop/floop
    children: List['PPTNode'] = field(default_factory=list)


_LEAF_RE = re.compile(r'^"(?P<name>.*)":(?P<weight>[-\d.eE+]+)\s*$')
_TAU_RE = re.compile(r'^tau:(?P<weight>[-\d.eE+]+)\s*$')
_NODE1_RE = re.compile(r'^(?P<op>PLoop|FLoop)\[(?P<rho>[-\d.eE+]+)\]:(?P<weight>[-\d.eE+]+)\s*$')
_NODEN_RE = re.compile(r'^(?P<op>Seq|Choice|Conc):(?P<weight>[-\d.eE+]+)\s*$')

_NODEN_KIND = {'Seq': 'seq', 'Choice': 'choice', 'Conc': 'conc'}
_NODE1_KIND = {'PLoop': 'ploop', 'FLoop': 'floop'}


def parse_ptree(text: str) -> PPTNode:
    """Parses Toothpaste's .ptree text export format into a PPTNode tree."""
    lines = [l for l in text.split('\n') if l.strip() != '']
    if not lines:
        raise ValueError("Empty .ptree input")
    node, next_index = _parse_lines(lines, 0, 0)
    if next_index != len(lines):
        raise ValueError(f"Unexpected trailing content starting at line {next_index}: {lines[next_index]!r}")
    return node


def _indent_of(line: str) -> int:
    stripped = line.lstrip(' ')
    spaces = len(line) - len(stripped)
    if spaces % 2 != 0:
        raise ValueError(f"Odd indentation (expected 2 spaces per level): {line!r}")
    return spaces // 2


def _parse_lines(lines: List[str], index: int, indent: int) -> Tuple[PPTNode, int]:
    line = lines[index]
    line_indent = _indent_of(line)
    if line_indent != indent:
        raise ValueError(f"Expected indent {indent}, got {line_indent} at line {index}: {line!r}")
    content = line.strip()

    m = _LEAF_RE.match(content)
    if m:
        return PPTNode('leaf', float(m.group('weight')), name=m.group('name')), index + 1

    m = _TAU_RE.match(content)
    if m:
        return PPTNode('tau', float(m.group('weight'))), index + 1

    m = _NODE1_RE.match(content)
    if m:
        child, next_index = _parse_lines(lines, index + 1, indent + 1)
        node = PPTNode(_NODE1_KIND[m.group('op')], float(m.group('weight')),
                        rho=float(m.group('rho')), children=[child])
        return node, next_index

    m = _NODEN_RE.match(content)
    if m:
        children = []
        next_index = index + 1
        while next_index < len(lines) and _indent_of(lines[next_index]) > indent:
            child, next_index = _parse_lines(lines, next_index, indent + 1)
            children.append(child)
        if not children:
            raise ValueError(f"{m.group('op')} node at line {index} has no children")
        node = PPTNode(_NODEN_KIND[m.group('op')], float(m.group('weight')), children=children)
        return node, next_index

    raise ValueError(f"Could not parse .ptree line: {content!r}")


class _IdGen:
    def __init__(self):
        self._n = 0

    def next(self) -> str:
        self._n += 1
        return f"ppt{self._n}"


def translate_ppt(ppt: PPTNode, model_move_cost: int = DEFAULT_MODEL_MOVE_COST) -> Tuple[ProcessTree, Dict[str, float], List[Tuple[str, str]]]:
    """
    Translates a parsed PPTNode tree into a skip-alignments ProcessTree,
    returning (tree, weights, loop_taus):
      - weights maps each leaf/tau id -- including the synthetic
        Tau_skip/Tau_redo introduced by the PLoop translation -- to its
        derived SLPN weight.
      - loop_taus is a list of (tau_skip_id, tau_redo_id) pairs, one per
        translated PLoop, needed by find_loop_structural_transitions to
        locate and weight the two unlabelled Petri-net transitions each
        translated loop additionally requires.
    See ppt_translation.md for the derivation.
    """
    id_gen = _IdGen()
    weights: Dict[str, float] = {}
    loop_taus: List[Tuple[str, str]] = []
    tree = _translate(ppt, None, id_gen, model_move_cost, weights, loop_taus)
    return tree, weights, loop_taus


def _translate(ppt: PPTNode, parent: Optional[ProcessTree], id_gen: _IdGen,
               model_move_cost: int, weights: Dict[str, float],
               loop_taus: List[Tuple[str, str]]) -> ProcessTree:
    if ppt.kind == 'leaf':
        node = Activity(parent, ppt.name, model_move_cost)
        node.id = id_gen.next()
        weights[node.id] = ppt.weight
        return node

    if ppt.kind == 'tau':
        node = Tau(parent, 'tau', model_move_cost)
        node.id = id_gen.next()
        weights[node.id] = ppt.weight
        return node

    if ppt.kind in ('seq', 'choice', 'conc'):
        cls = {'seq': Sequence, 'choice': Xor, 'conc': And}[ppt.kind]
        node = cls(parent, [])
        node.id = id_gen.next()
        # Recorded on every node, not just leaves, so a compound node can be
        # weighted as an Xor branch the same way a leaf is -- see
        # compile_to_slpn's Xor case. PPT's own weight is defined on every
        # node kind, this just stops translate_ppt from discarding it for
        # anything but leaves/taus.
        weights[node.id] = ppt.weight
        node.children = [_translate(c, node, id_gen, model_move_cost, weights, loop_taus) for c in ppt.children]
        return node

    if ppt.kind == 'floop':
        # Fixed loop: syntactic shorthand for a Sequence repeating the child
        # round(rho) times, each copy inheriting the loop's own weight
        # unchanged -- PPT's own stated rule for fixed loops, no ratio
        # split (unlike PLoop, which needs the geometric-distribution split).
        if len(ppt.children) != 1:
            raise ValueError("FLoop must have exactly one child")
        count = round(ppt.rho)
        node = Sequence(parent, [])
        node.id = id_gen.next()
        weights[node.id] = ppt.weight
        node.children = [_translate(ppt.children[0], node, id_gen, model_move_cost, weights, loop_taus)
                          for _ in range(count)]
        return node

    if ppt.kind == 'ploop':
        if len(ppt.children) != 1:
            raise ValueError("PLoop must have exactly one child")
        w = ppt.weight
        rho = ppt.rho
        continue_factor = (rho - 1) / rho

        xor = Xor(parent, [])
        xor.id = id_gen.next()
        weights[xor.id] = w

        tau_skip = Tau(xor, 'tau_skip', model_move_cost)
        tau_skip.id = id_gen.next()
        weights[tau_skip.id] = w / rho

        loop = Loop(xor, [])
        loop.id = id_gen.next()
        # Matches Tau_redo's weight below (both are the Loop rule's "child
        # inherits the loop's own weight") -- recorded under the Loop
        # node's own id so compile_to_slpn's Xor case can look up "the
        # weight of choosing this branch" the same way it does for any
        # other Xor child, without special-casing PLoop-translated Loops.
        weights[loop.id] = w * continue_factor

        x = _translate(ppt.children[0], loop, id_gen, model_move_cost, weights, loop_taus)
        # The child inherits the *translated* loop's own weight (w(rho-1)/rho,
        # not PPT's own w), per the Loop rule -- but if x is itself a compound
        # subtree, its descendants' weights were computed relative to PPT's
        # original (unscaled) weight, so the whole subtree needs rescaling by
        # the same factor to stay internally consistent, not just its root.
        # (Mirrors Toothpaste's own `scale` operation in ProbProcessTree.hs.)
        _scale_weights(x, weights, continue_factor)

        tau_redo = Tau(loop, 'tau_redo', model_move_cost)
        tau_redo.id = id_gen.next()
        weights[tau_redo.id] = w * continue_factor

        loop.children = [x, tau_redo]
        xor.children = [tau_skip, loop]
        loop_taus.append((tau_skip.id, tau_redo.id))
        return xor

    raise ValueError(f"Unknown PPT node kind: {ppt.kind!r}")


def _scale_weights(node: ProcessTree, weights: Dict[str, float], factor: float) -> None:
    if node.id in weights:
        weights[node.id] *= factor
    for c in node.children:
        _scale_weights(c, weights, factor)


def compile_to_slpn(tree: ProcessTree, weights: Dict[str, float],
                     loop_taus: List[Tuple[str, str]]) -> Tuple[str, Dict[str, str]]:
    """
    Compiles tree directly into Ebi's plaintext .slpn format ('stochastic
    labelled Petri net', confirmed against real `ebi discover uniform`
    output -- see ppt_translation.md), with every transition's weight
    attached at the moment it's created. No Petri-net library and no Ebi
    subprocess call is involved: every operator (Sequence/Xor/And/Loop/
    Activity/Tau) has a standard block-structured Petri-net compilation
    (single entry place, single exit place per block), so this walks `tree`
    directly instead of delegating structure to pm4py and reattaching
    weights afterward by searching for them. A transition can't end up
    without a weight -- it's a parameter to whatever code creates the
    transition, not a separate lookup step.

    `weights` (from translate_ppt) must carry an entry for every leaf/tau
    id and every node that's a direct child of an Xor (its own selection
    weight, including translated Loop nodes -- see translate_ppt's ploop
    branch). `loop_taus` ties each translated Loop node to its
    (tau_skip_id, tau_redo_id) pair, used here to find the loop's own exit
    weight (the sibling Tau_skip's weight, per ppt_translation.md's table)
    via the Loop's own redo child.

    Returns (slpn_text, activity_to_id): activity_to_id maps each Activity's
    name to the tree id chosen to represent it (first-seen wins for
    repeated activity names).
    """
    places: List[int] = []
    transitions: List[dict] = []
    activity_to_id: Dict[str, str] = {}
    redo_to_skip = {redo_id: skip_id for skip_id, redo_id in loop_taus}

    def new_place() -> int:
        p = len(places)
        places.append(p)
        return p

    def add_transition(label: Optional[str], weight: float, in_places: List[int], out_places: List[int]) -> None:
        transitions.append({'label': label, 'weight': weight, 'in': in_places, 'out': out_places})

    def compile_node(node: ProcessTree, entry: int, exit: Optional[int] = None) -> int:
        ex = exit if exit is not None else new_place()

        if isinstance(node, Activity):
            activity_to_id.setdefault(node.name, node.id)
            add_transition(node.id, weights[node.id], [entry], [ex])
            return ex

        if isinstance(node, Tau):
            add_transition(node.id, weights[node.id], [entry], [ex])
            return ex

        if isinstance(node, Sequence):
            if not node.children:
                raise ValueError("Sequence with no children")
            cur = entry
            for i, c in enumerate(node.children):
                cur = compile_node(c, cur, ex if i == len(node.children) - 1 else None)
            return cur

        if isinstance(node, Xor):
            if not node.children:
                raise ValueError("Xor with no children")
            for c in node.children:
                if isinstance(c, LeafNode):
                    # the leaf's own transition already carries its
                    # selection weight -- no separate gate needed
                    compile_node(c, entry, ex)
                else:
                    gate_exit = new_place()
                    add_transition(None, weights[c.id], [entry], [gate_exit])
                    compile_node(c, gate_exit, ex)
            return ex

        if isinstance(node, And):
            split_targets = [new_place() for _ in node.children]
            add_transition(None, 1, [entry], split_targets)
            join_sources = [compile_node(c, p) for c, p in zip(node.children, split_targets)]
            add_transition(None, 1, join_sources, [ex])
            return ex

        if isinstance(node, Loop):
            if len(node.children) != 2:
                raise ValueError("Loop must have exactly two children (do, redo)")
            do, redo = node.children
            mid = new_place()
            compile_node(do, entry, mid)
            compile_node(redo, mid, entry)
            tau_skip_id = redo_to_skip[redo.id]
            add_transition(None, weights[tau_skip_id], [mid], [ex])
            return ex

        raise ValueError(f"Unsupported ProcessTree node type: {type(node)!r}")

    entry = new_place()
    compile_node(tree, entry)

    lines = ['stochastic labelled Petri net', '# number of places', str(len(places)), '# initial marking']
    lines.extend('1' if p == entry else '0' for p in places)
    lines.append('# number of transitions')
    lines.append(str(len(transitions)))
    for i, t in enumerate(transitions):
        lines.append(f'# transition {i}')
        lines.append(f"label {t['label']}" if t['label'] is not None else 'silent')
        lines.append('# weight')
        lines.append(str(t['weight']))
        lines.append('# number of input places')
        lines.append(str(len(t['in'])))
        lines.extend(str(p) for p in t['in'])
        lines.append('# number of output places')
        lines.append(str(len(t['out'])))
        lines.extend(str(p) for p in t['out'])
    return '\n'.join(lines), activity_to_id


def write_slpn(tree: ProcessTree, weights: Dict[str, float], loop_taus: List[Tuple[str, str]],
                out: str = 'smodel.slpn') -> Dict[str, str]:
    """
    Compiles tree to .slpn via compile_to_slpn and writes it to `out`. The
    only file-system/Ebi-adjacent step in the whole PPT-import path --
    everything upstream (parse_ptree, translate_ppt, compile_to_slpn) is
    pure. Returns activity_to_id, as compile_to_slpn does.
    """
    text, activity_to_id = compile_to_slpn(tree, weights, loop_taus)
    with open(out, 'w') as f:
        f.write(text)
    return activity_to_id
