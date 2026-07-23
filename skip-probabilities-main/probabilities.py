from typing import Dict
from processtree import *
from alignment import *
import random
import pandas as pd
import subprocess
from tqdm import tqdm
    
class EbiOccurance(object):
    def __init__(self):
        pass

    def _id_to_activity(self, tree:ProcessTree, id:str):
        if tree.id == id:
            assert isinstance(tree, LeafNode)
            return tree.name
        for c in tree.children:
            res = self._id_to_activity(c, id)
            if res is not None:
                return res
        return None

    def write_tree_to_petri(self, tree:ProcessTree) -> Dict[str, str]:
        # transitions:
        # activity -> activity_id
        # tau -> tau_id
        # helper-tau -> silent
        # NOTE: activities with the same label are assigned the same activity_id
        # returns a dict activity -> selected activity_id
        tree_pm4py = tree.to_pm4py(True)
        net, im, fm = pm4py.convert.convert_to_petri_net(tree_pm4py)
        activity_to_id = {}
        for t in net.transitions:
            # t.label is the id
            if self._id_to_activity(tree, t.label) in activity_to_id:
                t.label = activity_to_id[self._id_to_activity(tree, t.label)]
            else:
                activity_to_id[self._id_to_activity(tree, t.label)] = t.label
        pm4py.write_pnml(net, im, fm, 'model.pnml')
        return activity_to_id
    
    def write_log(self, log:pd.DataFrame, rename_dict:Dict[str, str]) -> pd.DataFrame:
        # log is a pm4py event log
        # rename_dict is the output of write_tree_to_petri used to rename the acitvities in the log with a uniform id
        log = log.copy()
        log['concept:name'] = log['concept:name'].apply(lambda x: rename_dict[x] if x in rename_dict else x)
        pm4py.write_xes(log, 'log.xes')
        return log
    
    def ebi_slpn(self, model='model.pnml', log='log.xes', out='smodel.slpn'):
        subprocess.check_call([r"ebi.exe", "disc", "occ", log, model, "-o", out])
    
    def validate_slpn(self, tree:ProcessTree, path='smodel.slpn'):
        file = open(path,"r")
        lines = file.readlines()
        file.close()
        seen_ids = []
        fix_needed = False
        for i in range(len(lines)-1):
            if lines[i].startswith("# transition") and lines[i+1].startswith("label"):
                id = lines[i+1][6:]
                if id in seen_ids:
                    print("There is at least one duplicate label that needs to be fixed. Use the graphical models to correct the names. Missing:", id)
                    self.reconstruct_petri_net(tree, lines)
                    fix_needed = True
                    break
                    raise ValueError("Duplicate transition ids not cleared:", id)
                else:
                    seen_ids.append(id)
        if fix_needed:
            self.question_duplicates(tree, path)
        return

    def question_duplicates(self, tree:ProcessTree, path):
        input("Please fix duplicate transitions. If you think cou are done, press enter.")
        self.validate_slpn(tree, path)
    
    def reconstruct_petri_net(self, tree:ProcessTree, lines:List[str]):
        num_places = int(lines[2])
        num_transitions = int(lines[2+num_places+3])
        places = [pm4py.objects.petri_net.obj.PetriNet.Place(str(i), [], []) for i in range(num_places)]
        transitions = [pm4py.objects.petri_net.obj.PetriNet.Transition(str(i), None, [], []) for i in range(num_transitions)]
        arcs = []
        transition_counter = 0
        for i in range(len(lines)-1):
            if lines[i].startswith("# transition"):
                if lines[i+1].startswith("label"):
                    transitions[transition_counter].label = lines[i].strip()# + "\n" + lines[i+1].strip().split("label ")[1]# + " l=" + str(i+1)
                # places in (ignore first line)
                j=1
                while i+5+j < len(lines) and not lines[i+5+j].startswith("#"):
                    place_idx = int(lines[i+5+j])
                    arc = pm4py.objects.petri_net.obj.PetriNet.Arc(places[place_idx], transitions[transition_counter])
                    arcs.append(arc)
                    places[place_idx].out_arcs.append(arc)
                    transitions[transition_counter].in_arcs.append(arc)
                    j += 1
                # places in (ignore first line)
                j += 2
                while i+5+j < len(lines) and not lines[i+5+j].startswith("#"):
                    place_idx = int(lines[i+5+j])
                    arc = pm4py.objects.petri_net.obj.PetriNet.Arc(transitions[transition_counter], places[place_idx])
                    arcs.append(arc)
                    places[place_idx].in_arcs.append(arc)
                    transitions[transition_counter].out_arcs.append(arc)
                    j += 1
                transition_counter += 1
        pn = pm4py.objects.petri_net.obj.PetriNet("Net", places, transitions, arcs)
        pm4py.view_petri_net(pn, format='png')

        tree_pm4py = tree.to_pm4py(True)
        net, im, fm = pm4py.convert.convert_to_petri_net(tree_pm4py)
        pm4py.view_petri_net(net, format='png')
                    
    
    def update_visible_taus(self, tree:ProcessTree, path='smodel.slpn'):
        self.validate_slpn(tree, path)
        file = open(path,"r")
        lines = file.readlines()
        file.close()
        for i in range(len(lines)-1):
            if lines[i].startswith("# transition") and lines[i+1].startswith("label TAU_"):
                lines[i+3] = "1" + lines[i+3][1:]
        file = open(path, "w")
        file.writelines(lines)
        file.close()
        
    def ebi_trace_prob(self, trace:List[str], model='smodel.slpn'):
        time_start = time.process_time_ns()
        res = subprocess.check_output([r"ebi.exe", "prob", "trac", model] + trace + ["-a"]).decode("utf-8")
        time_end = time.process_time_ns()
        if not 'Approximately' in res:
            raise ValueError("Ebi did not return a probability")
        return float(res.split('\n')[0]), (time_end-time_start)
    
    def trace_probs(self, agns:Dict[str, List[List|tuple]], measure:Optional[Dict[List[str], float]]=None, model='smodel.slpn'):
        # input: agns:    var -> List[agn]
        #        measure: var id list -> model prob
        # output: pi2(agn) IDS -> float, var -> pi2(agn) -> int
        trace_counts = {}
        trace_probs_d = {}
        for var, ass in tqdm(agns.items()):
            if var not in trace_counts:
                trace_counts[var] = {}
            for agn in ass:
                model_path = tuple([m for m in list(zip(*agn))[1] if m != '>>'])
                if model_path not in trace_counts[var]:
                    trace_counts[var][model_path] = 1
                else:
                    trace_counts[var][model_path] = trace_counts[var][model_path] + 1
        if measure is not None:
            for var, ass in tqdm(agns.items()):
                for agn in ass:
                    model_path = tuple([m for m in list(zip(*agn))[1] if m != '>>'])
                    model_path_ids = tuple([n.id for n in model_path])
                    if not model_path_ids in trace_probs_d:
                        trace_probs_d[model_path_ids] = measure[model_path_ids]
            return trace_probs_d, trace_counts, -1
        
        from concurrent.futures import ProcessPoolExecutor, as_completed
        with ProcessPoolExecutor(max_workers=14) as executor:
            futures = []
            checked_ids = []
            for var, ass in tqdm(agns.items()):
                for agn in ass:
                    model_path = tuple([m for m in list(zip(*agn))[1] if m != '>>'])
                    model_path_ids = tuple([n.id for n in model_path])
                    if not model_path_ids in checked_ids:
                        futures.append(executor.submit(self.ebi_trace_prob, list(model_path_ids), model))
                        checked_ids.append(model_path_ids)
            print("Futures created")
            progress = tqdm(total=len(futures))
            for future in as_completed(futures):
                progress.update()
            total_time = 0
            for index, _ in enumerate(futures):
                trace_prob, needed_time = futures[index].result()
                trace_probs_d[checked_ids[index]] = trace_prob
                total_time += needed_time
        return trace_probs_d, trace_counts, total_time
    
    def _skip_agn_probs(self, state:State, node:ProcessTree, var:str, agns:Set[List|tuple], trace_probs:Dict[List[LeafNode], float], trace_counts:Dict[List[LeafNode], int]):
        prob = 0 # skip agn prob
        if not node.id in [e.node.id for e in state.executions]+[m.node.id for l,m in state.path if isinstance(m,Skip)]+[m.node.id for l,m in state.path if isinstance(m,TauPath)]:
            return prob
        for agn in agns:
            model_path = tuple([m for m in list(zip(*agn))[1] if m != '>>'])
            model_path_ids = tuple([n.id for n in model_path])
            prob += trace_probs[model_path_ids]/trace_counts[var][model_path]
        return prob
    
    def skip_agn_probs_per_node(self, node:ProcessTree, skip_dict:Dict[str,List[State]], C:Dict[State,Set[List|tuple]], trace_probs:Dict[tuple,float], trace_counts:Dict[tuple,int]):
        # input: node:         node to calculate for
        #        skip_dict:    var -> List[sagn]
        #        C:            sagn -> Set[agn] coinciding agns
        #        trace_probs:  pi2(agn) -> float
        #        trace_counts: pi2(agn) -> int
        # output: score_sagn:  sagn -> sum(partial prob of coinciding agns)
        #         cond_prob:   sagn -> P(sagn|var)
        score_sagn = {}
        cond_prob = {}
        for var, states in skip_dict.items():
            for state in states:
                score_sagn[state] = self._skip_agn_probs(state, node, var, C[state], trace_probs, trace_counts)
            for state in states:
                cond_prob[state] = score_sagn[state]/sum(score_sagn[s] for s in states) if sum(score_sagn[s] for s in states) > 0 else 0.0
        return score_sagn, cond_prob
    
    def skip_agn_probs_traversal(self, tree:ProcessTree, skip_dict:Dict[str,List[State]], C:Dict[State,Set[List|tuple]], trace_probs:Dict[tuple,float], trace_counts:Dict[tuple,int]):
        # input: tree:         process tree whose nodes are alanyzed
        #        skip_dict:    var -> List[sagn]
        #        C:            sagn -> Set[agn] coinciding agns
        #        trace_probs:  pi2(agn) -> float
        #        trace_counts: pi2(agn) -> int
        # output: score_sagn:  sagn -> sum(partial prob of coinciding agns)
        #         cond_prob:   sagn -> P(sagn|var)
        score_sagn = {}
        cond_prob = {}
        if isinstance(tree, LeafNode):
            score, cprob = self.skip_agn_probs_per_node(tree, skip_dict, C, trace_probs, trace_counts)
            score_sagn[tree] = score
            cond_prob[tree] = cprob
            return score_sagn, cond_prob
        else:
            score, cprob = self.skip_agn_probs_per_node(tree, skip_dict, C, trace_probs, trace_counts)
            score_sagn[tree] = score
            cond_prob[tree] = cprob
            for c in tree.children:
                nscore, ncprob = self.skip_agn_probs_traversal(c, skip_dict, C, trace_probs, trace_counts)
                for k,v in nscore.items():
                    score_sagn[k] = v
                for k,v in ncprob.items():
                    cond_prob[k] = v
            return score_sagn, cond_prob
        
    
    # def skip_agn_probs(self, tree:ProcessTree, skip_dict:Dict[str,List[State]], C:Dict[State,Set[List|tuple]], trace_probs:Dict[tuple,float], trace_counts:Dict[tuple,int]):
    #     # input: skip_dict:    var -> List[sagn]
    #     #        C:            sagn -> Set[agn] coinciding agns
    #     #        trace_probs:  pi2(agn) -> float
    #     #        trace_counts: pi2(agn) -> int
    #     # output: score_sagn:  sagn -> sum(partial prob of coinciding agns)
    #     #         cond_prob:   sagn -> P(sagn|var)
    #     score_sagn = {}
    #     cond_prob = {}
    #     for var, states in tqdm(skip_dict.items()):
    #         for state in states:
    #             score_sagn[node][state] = EbiOccurance()._skip_agn_probs(state, node, var, C[state], trace_probs, trace_counts)
    #         for state in states:
    #             cond_prob[node][state] = score_sagn[node][state]/sum(score_sagn[node][s] for s in states) if sum(score_sagn[node][s] for s in states) > 0 else 0.0
    #     return score_sagn, cond_prob