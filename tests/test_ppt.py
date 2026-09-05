"""
Unit tests for skip_alignments.ppt: parsing Toothpaste Miner's .ptree export
format and translating it into skip-alignments' own ProcessTree classes with
derived SLPN weights.

Format and weight formulas confirmed against refs/toothpaste-master
(Toothpaste's own Haskell source) -- see ppt_translation.md for the
derivation.

Run with:
    python -m unittest tests.test_ppt -v
"""
import unittest

from skip_alignments.processtree import Activity, And, Loop, Sequence, Tau, Xor
from skip_alignments.ppt import PPTNode, compile_to_slpn, parse_ptree, translate_ppt


class TestParsePTree(unittest.TestCase):

    def test_parse_single_leaf(self):
        ppt = parse_ptree('"a":5.0\n')
        self.assertEqual(ppt.kind, 'leaf')
        self.assertEqual(ppt.name, 'a')
        self.assertEqual(ppt.weight, 5.0)
        self.assertEqual(ppt.children, [])

    def test_parse_tau(self):
        ppt = parse_ptree('tau:3.0\n')
        self.assertEqual(ppt.kind, 'tau')
        self.assertEqual(ppt.weight, 3.0)

    def test_parse_sequence_of_leaves(self):
        text = 'Seq:10.0\n  "a":10.0\n  "b":10.0\n'
        ppt = parse_ptree(text)
        self.assertEqual(ppt.kind, 'seq')
        self.assertEqual(ppt.weight, 10.0)
        self.assertEqual([c.kind for c in ppt.children], ['leaf', 'leaf'])
        self.assertEqual([c.name for c in ppt.children], ['a', 'b'])

    def test_parse_ploop_wrapping_single_leaf(self):
        # the shape actually seen in real Toothpaste output
        # (results/2021_pn/teleclaims_k1.ptree): PLoop directly wrapping one leaf
        text = 'PLoop[2.0]:344.0\n  "some activity":344.0\n'
        ppt = parse_ptree(text)
        self.assertEqual(ppt.kind, 'ploop')
        self.assertEqual(ppt.weight, 344.0)
        self.assertEqual(ppt.rho, 2.0)
        self.assertEqual(len(ppt.children), 1)
        self.assertEqual(ppt.children[0].name, 'some activity')

    def test_parse_floop(self):
        text = 'FLoop[3.0]:7.0\n  "a":7.0\n'
        ppt = parse_ptree(text)
        self.assertEqual(ppt.kind, 'floop')
        self.assertEqual(ppt.rho, 3.0)

    def test_parse_choice_and_conc(self):
        text = (
            'Choice:8.0\n'
            '  "a":5.0\n'
            '  "b":3.0\n'
        )
        ppt = parse_ptree(text)
        self.assertEqual(ppt.kind, 'choice')
        self.assertEqual(len(ppt.children), 2)

        text2 = 'Conc:8.0\n  "a":4.0\n  "b":4.0\n'
        ppt2 = parse_ptree(text2)
        self.assertEqual(ppt2.kind, 'conc')

    def test_parse_nested_real_shape(self):
        # a trimmed, representative excerpt of the real sample's shape:
        # Seq -> [leaf, Choice -> [PLoop -> leaf, tau]]
        text = (
            'Seq:704.0\n'
            '  "incoming claim":704.0\n'
            '  Choice:704.0\n'
            '    PLoop[2.0]:344.0\n'
            '      "check info":344.0\n'
            '    tau:360.0\n'
        )
        ppt = parse_ptree(text)
        self.assertEqual(ppt.kind, 'seq')
        self.assertEqual(ppt.children[0].kind, 'leaf')
        choice = ppt.children[1]
        self.assertEqual(choice.kind, 'choice')
        self.assertEqual(choice.children[0].kind, 'ploop')
        self.assertEqual(choice.children[0].children[0].name, 'check info')
        self.assertEqual(choice.children[1].kind, 'tau')


class TestTranslatePPTStructure(unittest.TestCase):

    def test_leaf_translates_to_activity(self):
        ppt = PPTNode('leaf', 5.0, name='a')
        tree, weights, loop_taus = translate_ppt(ppt)
        self.assertIsInstance(tree, Activity)
        self.assertEqual(tree.name, 'a')
        self.assertAlmostEqual(weights[tree.id], 5.0)

    def test_tau_translates_to_tau(self):
        ppt = PPTNode('tau', 3.0)
        tree, weights, loop_taus = translate_ppt(ppt)
        self.assertIsInstance(tree, Tau)
        self.assertAlmostEqual(weights[tree.id], 3.0)

    def test_seq_choice_conc_map_to_sequence_xor_and(self):
        for kind, cls in [('seq', Sequence), ('choice', Xor), ('conc', And)]:
            ppt = PPTNode(kind, 10.0, children=[
                PPTNode('leaf', 10.0, name='a'),
                PPTNode('leaf', 10.0, name='b'),
            ])
            tree, weights, loop_taus = translate_ppt(ppt)
            self.assertIsInstance(tree, cls)
            self.assertEqual(len(tree.children), 2)
            self.assertIs(tree.children[0].parent, tree)

    def test_floop_translates_to_sequence_of_repeated_child(self):
        ppt = PPTNode('floop', 7.0, rho=3.0, children=[PPTNode('leaf', 7.0, name='a')])
        tree, weights, loop_taus = translate_ppt(ppt)
        self.assertIsInstance(tree, Sequence)
        self.assertEqual(len(tree.children), 3)
        for c in tree.children:
            self.assertIsInstance(c, Activity)
            self.assertEqual(c.name, 'a')
            self.assertAlmostEqual(weights[c.id], 7.0)
        # each copy must have a distinct id
        self.assertEqual(len({c.id for c in tree.children}), 3)

    def test_all_ids_unique(self):
        ppt = PPTNode('seq', 10.0, children=[
            PPTNode('leaf', 10.0, name='a'),
            PPTNode('ploop', 10.0, rho=2.0, children=[PPTNode('leaf', 10.0, name='b')]),
        ])
        tree, weights, loop_taus = translate_ppt(ppt)
        ids = _all_ids(tree)
        self.assertEqual(len(ids), len(set(ids)))


class TestModelMoveCosts(unittest.TestCase):
    """
    The alignment engine's own invariant (alignment.py's Aligner.align2:
    `assert tau_cost < activity_cost`, matching the codebase-wide convention
    of model_move_activity_cost=100000 / model_move_tau_cost=0 used
    everywhere else, e.g. ProcessTree.from_pm4py's callers) was never
    honoured by translate_ppt -- it reused the same model_move_cost for
    both Activity and Tau nodes, so every real alignment run against a
    translated PPT tree hit that assertion. Caught by an end-to-end
    DerivationPipeline test, not by translate_ppt's own isolated tests.
    """

    def test_activity_leaf_gets_default_activity_cost(self):
        ppt = PPTNode('leaf', 5.0, name='a')
        tree, weights, loop_taus = translate_ppt(ppt)
        self.assertEqual(tree.skip_cost, 100000)

    def test_tau_leaf_gets_zero_cost_by_default(self):
        ppt = PPTNode('tau', 3.0)
        tree, weights, loop_taus = translate_ppt(ppt)
        self.assertEqual(tree.skip_cost, 0)

    def test_ploop_synthetic_taus_get_zero_cost_by_default(self):
        w, rho = 344.0, 2.0
        ppt = PPTNode('ploop', w, rho=rho, children=[PPTNode('leaf', w, name='check info')])
        tree, weights, loop_taus = translate_ppt(ppt)
        tau_skip, loop = tree.children
        x, tau_redo = loop.children
        self.assertEqual(tau_skip.skip_cost, 0)
        self.assertEqual(tau_redo.skip_cost, 0)
        self.assertEqual(x.skip_cost, 100000)

    def test_costs_are_overridable_and_stay_ordered(self):
        ppt = PPTNode('ploop', 10.0, rho=2.0, children=[PPTNode('leaf', 10.0, name='a')])
        tree, weights, loop_taus = translate_ppt(ppt, model_move_cost=500, model_move_tau_cost=7)
        tau_skip, loop = tree.children
        x, tau_redo = loop.children
        self.assertEqual(x.skip_cost, 500)
        self.assertEqual(tau_skip.skip_cost, 7)
        self.assertEqual(tau_redo.skip_cost, 7)
        self.assertLess(tau_skip.skip_cost, x.skip_cost)


class TestTranslatePLoopWeights(unittest.TestCase):
    """
    See ppt_translation.md for the derivation. For a PLoop[rho] node with
    weight w wrapping child x, the translation Xor(Tau_skip, Loop(x, Tau_redo))
    should carry:
        weight(Tau_skip) = w/rho
        weight(x)         = w(rho-1)/rho
        weight(Tau_redo)  = w(rho-1)/rho
    weight(loop.id) is also set to w(rho-1)/rho (the Loop rule again,
    recorded under the Loop node's own id so it can be looked up uniformly
    as an Xor branch weight by compile_to_slpn -- see TestCompileToSlpn
    below for the two additional unlabelled structural transitions
    compile_to_slpn derives from this).
    """

    def test_ploop_wrapping_leaf(self):
        w, rho = 344.0, 2.0
        ppt = PPTNode('ploop', w, rho=rho, children=[PPTNode('leaf', w, name='check info')])
        tree, weights, loop_taus = translate_ppt(ppt)

        self.assertIsInstance(tree, Xor)
        self.assertEqual(len(tree.children), 2)
        tau_skip, loop = tree.children
        self.assertIsInstance(tau_skip, Tau)
        self.assertIsInstance(loop, Loop)
        self.assertEqual(len(loop.children), 2)
        x, tau_redo = loop.children
        self.assertIsInstance(x, Activity)
        self.assertEqual(x.name, 'check info')
        self.assertIsInstance(tau_redo, Tau)

        self.assertAlmostEqual(weights[tau_skip.id], w / rho)
        expected_continue = w * (rho - 1) / rho
        self.assertAlmostEqual(weights[x.id], expected_continue)
        self.assertAlmostEqual(weights[tau_redo.id], expected_continue)
        self.assertAlmostEqual(weights[loop.id], expected_continue)

        self.assertEqual(loop_taus, [(tau_skip.id, tau_redo.id)])

    def test_ploop_weights_conserve_total_mass(self):
        # Xor rule: children's weights sum to the parent's -- check this
        # holds for the outer Xor(Tau_skip, Loop-branch) split specifically,
        # using the loop-branch's own weight (w(rho-1)/rho, carried by its
        # children per the Loop rule) as the second term.
        w, rho = 100.0, 4.0
        ppt = PPTNode('ploop', w, rho=rho, children=[PPTNode('leaf', w, name='a')])
        tree, weights, loop_taus = translate_ppt(ppt)
        tau_skip, loop = tree.children
        x, tau_redo = loop.children
        self.assertAlmostEqual(weights[tau_skip.id] + weights[x.id], w)

    def test_ploop_different_rho_values(self):
        for rho in (1.5, 2.0, 5.0, 10.0):
            w = 50.0
            ppt = PPTNode('ploop', w, rho=rho, children=[PPTNode('leaf', w, name='a')])
            tree, weights, loop_taus = translate_ppt(ppt)
            tau_skip, loop = tree.children
            x, tau_redo = loop.children
            with self.subTest(rho=rho):
                self.assertAlmostEqual(weights[tau_skip.id], w / rho)
                self.assertAlmostEqual(weights[x.id], w * (rho - 1) / rho)
                self.assertAlmostEqual(weights[tau_redo.id], w * (rho - 1) / rho)

    def test_nested_ploop_inside_sequence(self):
        # matches the real sample's shape: Seq[a, PLoop[rho](b)]
        w = 704.0
        rho = 2.0
        ppt = PPTNode('seq', w, children=[
            PPTNode('leaf', w, name='a'),
            PPTNode('ploop', w, rho=rho, children=[PPTNode('leaf', w, name='b')]),
        ])
        tree, weights, loop_taus = translate_ppt(ppt)
        self.assertIsInstance(tree, Sequence)
        a, xor = tree.children
        self.assertIsInstance(a, Activity)
        self.assertIsInstance(xor, Xor)
        tau_skip, loop = xor.children
        self.assertAlmostEqual(weights[tau_skip.id], w / rho)


def _parse_slpn(text):
    """Parses compile_to_slpn's output back into a list of
    (label_or_None, weight, in_places, out_places) per transition, for
    assertions -- independent of the writer's own internals."""
    lines = text.split('\n')
    assert lines[0] == 'stochastic labelled Petri net'
    num_places = int(lines[2])
    marking = [int(x) for x in lines[4:4 + num_places]]
    i = 4 + num_places
    assert lines[i] == '# number of transitions'
    num_transitions = int(lines[i + 1])
    i += 2
    transitions = []
    for _ in range(num_transitions):
        assert lines[i].startswith('# transition')
        label = None if lines[i + 1] == 'silent' else lines[i + 1][len('label '):]
        assert lines[i + 2] == '# weight'
        weight = float(lines[i + 3])
        n_in = int(lines[i + 5])
        in_places = [int(x) for x in lines[i + 6:i + 6 + n_in]]
        j = i + 6 + n_in
        n_out = int(lines[j + 1])
        out_places = [int(x) for x in lines[j + 2:j + 2 + n_out]]
        transitions.append((label, weight, in_places, out_places))
        i = j + 2 + n_out
    return num_places, marking, transitions


class TestCompileToSlpn(unittest.TestCase):
    """
    compile_to_slpn attaches every transition's weight at the moment it's
    created (see ppt.py), so a translated PLoop's 4 weight-bearing
    transitions (2 labelled Tau_skip/Tau_redo, covered above, plus 2
    unlabelled structural ones for Loop/Xor's own control-flow routing --
    see ppt_translation.md) should all come out correctly weighted with no
    separate patch-in step.
    """

    def test_single_activity(self):
        ppt = PPTNode('leaf', 5.0, name='a')
        tree, weights, loop_taus = translate_ppt(ppt)
        text, activity_to_id = compile_to_slpn(tree, weights, loop_taus)
        num_places, marking, transitions = _parse_slpn(text)
        self.assertEqual(len(transitions), 1)
        label, weight, in_places, out_places = transitions[0]
        self.assertEqual(label, tree.id)
        self.assertEqual(weight, 5.0)
        self.assertEqual(activity_to_id, {'a': tree.id})

    def test_ploop_has_exactly_five_transitions_two_unlabelled(self):
        w, rho = 344.0, 2.0
        ppt = PPTNode('ploop', w, rho=rho, children=[PPTNode('leaf', w, name='check info')])
        tree, weights, loop_taus = translate_ppt(ppt)
        text, activity_to_id = compile_to_slpn(tree, weights, loop_taus)
        num_places, marking, transitions = _parse_slpn(text)

        # 3 labelled (tau_skip, tau_redo, the leaf) + 2 unlabelled: the
        # Xor's own gate into the Loop branch (entering the loop the first
        # time -- genuinely distinct from Tau_redo, which only fires when
        # looping back from inside; reusing one place for both would
        # wrongly re-enable Tau_skip after every iteration) and the loop's
        # own exit transition. No separate "harmless wrapping" transitions
        # exist here (unlike the old pm4py-based path), since a leaf Xor
        # branch (Tau_skip) uses its own transition as the branch gate
        # directly instead of getting a redundant extra one.
        self.assertEqual(len(transitions), 5)
        labels = [label for label, _, _, _ in transitions]
        self.assertEqual(labels.count(None), 2)
        self.assertEqual(set(l for l in labels if l is not None), set(loop_taus[0]) | {activity_to_id['check info']})

    def test_ploop_weights_match_formula_table(self):
        w, rho = 344.0, 2.0
        ppt = PPTNode('ploop', w, rho=rho, children=[PPTNode('leaf', w, name='check info')])
        tree, weights, loop_taus = translate_ppt(ppt)
        tau_skip_id, tau_redo_id = loop_taus[0]
        text, activity_to_id = compile_to_slpn(tree, weights, loop_taus)
        num_places, marking, transitions = _parse_slpn(text)

        stop_weight = w / rho
        continue_weight = w * (rho - 1) / rho
        leaf_id = activity_to_id['check info']

        by_label = {label: weight for label, weight, _, _ in transitions if label is not None}
        self.assertAlmostEqual(by_label[tau_skip_id], stop_weight)
        self.assertAlmostEqual(by_label[tau_redo_id], continue_weight)
        self.assertAlmostEqual(by_label[leaf_id], continue_weight)

        # 2 unlabelled: the Xor's own gate into the Loop branch (weighted
        # like Tau_redo -- both are the Loop rule's "continue" weight) and
        # the loop's own exit transition (weighted like Tau_skip -- both
        # are the "stop" weight).
        unlabelled_weights = sorted(w_ for label, w_, _, _ in transitions if label is None)
        self.assertEqual(len(unlabelled_weights), 2)
        self.assertAlmostEqual(unlabelled_weights[0], stop_weight)
        self.assertAlmostEqual(unlabelled_weights[1], continue_weight)

    def test_two_independent_loops_resolve_separately(self):
        w1, rho1 = 10.0, 2.0
        w2, rho2 = 20.0, 5.0
        ppt = PPTNode('seq', w1 + w2, children=[
            PPTNode('ploop', w1, rho=rho1, children=[PPTNode('leaf', w1, name='a')]),
            PPTNode('ploop', w2, rho=rho2, children=[PPTNode('leaf', w2, name='b')]),
        ])
        tree, weights, loop_taus = translate_ppt(ppt)
        self.assertEqual(len(loop_taus), 2)
        text, activity_to_id = compile_to_slpn(tree, weights, loop_taus)
        num_places, marking, transitions = _parse_slpn(text)

        # 2 unlabelled transitions per loop (gate + exit -- see
        # test_ploop_has_exactly_five_transitions_two_unlabelled), 4 total,
        # each pair carrying that loop's own stop/continue weights.
        unlabelled_weights = sorted(w_ for label, w_, _, _ in transitions if label is None)
        expected = sorted([w1 / rho1, w1 * (rho1 - 1) / rho1, w2 / rho2, w2 * (rho2 - 1) / rho2])
        self.assertEqual(unlabelled_weights, expected)

    def test_arc_indices_reference_valid_places_and_single_initial_token(self):
        w, rho = 344.0, 2.0
        ppt = PPTNode('ploop', w, rho=rho, children=[PPTNode('leaf', w, name='check info')])
        tree, weights, loop_taus = translate_ppt(ppt)
        text, activity_to_id = compile_to_slpn(tree, weights, loop_taus)
        num_places, marking, transitions = _parse_slpn(text)

        self.assertEqual(sum(marking), 1)
        for label, weight, in_places, out_places in transitions:
            for idx in in_places + out_places:
                self.assertTrue(0 <= idx < num_places)


def _all_ids(tree):
    ids = [tree.id]
    for c in tree.children:
        ids += _all_ids(c)
    return ids


if __name__ == '__main__':
    unittest.main()
