
from pathlib import Path
import pickle
from typing import Tuple
import pm4py
from skipalignments.alignment import *
from skipalignments.processtree import *
import pandas as pd
from tqdm import tqdm
from skipalignments.alignall import *
import statistics
import random
from tqdm import tqdm
from skipalignments.execution import *
from skipalignments.probabilities import *
from skipalignments.skips import Skipper

class EbiWeights(Enum):
    OCCURANCE = 1
    UNIFORM = 2

class DerivationPipeline(object):
    
    def __init__(self, tree:ProcessTree, aligned_log:Any, pl:Dict[Tuple[str], float]=None, pn_log:Any=None, pn_method:str=None, pn_measure=None, sagn_timeout=600):
        """
        Derivation pipeline for skip probabilities.

        tree: The process tree whose skip probabilities are computed
        aligned_log: The event log for which the probabilities are derived
        pl: Probability distribution on the trace variants in aligned_log; default: frequency distribution
        pn_log: The event log used to derive the model distribution
        pn_method: The stochastic method used to derive the model distribution; only used if pn_log is not None
        pn_measure: Probability distribution on the model paths
        sagn_timeout: Timeout for the skip alignment computation in s; default 10 min
        """
        self.tree = tree
        self.aligned_log = aligned_log
        if pl is not None:
            self.pl = pl
        else:
            self.pl = self.variant_prob_dist(self.aligned_log)
        if pn_log is not None and pn_method is not None:
            self.pn_log = pn_log
            self.pn_method = pn_method
        else:
            assert pn_measure is not None
            self.pn_measure = pn_measure
            self.pn_log = None
            self.pn_method = None
        self.sagn_timeout = sagn_timeout
        
        self.variants = self.get_variant_dict(self.aligned_log)
    
    def write(self, obj:Any, path:str, name:str):
        file = open(path + "/" + name,"wb")
        pickle.dump(obj, file)
        file.close()
    
    def compute(self, path:str, sagns_precomputed:Dict[str, List[State]]=None, slpn_path='smodel.slpn'):
        # outputs
        if path.endswith('/'):
            path = path[:-1]
        Path(path).mkdir(parents=True, exist_ok=True)

        self.write(self.tree, path, "tree")
        
        # compute skip alignments
        print("1/4\tComputing all optimal skip alignments in normal form.")
        if sagns_precomputed is None:
            skip_dict, skip_times = self.compute_skip_alignments(self.sagn_timeout)
            self.write(skip_dict, path, "skip_dict")
        else:
            skip_dict = sagns_precomputed
            skip_times = {k:-1 for k in skip_dict.keys()}
        self.skip_times = skip_times
        self.skip_dict = skip_dict
        self.rename_sagns(skip_dict)
        skip_dict_backup = {k:[s.copy() for s in v] for k,v in skip_dict.items()}
        self.skip_dict_backup = skip_dict_backup

        # compute skip alignment probabilities
        print("2/4\tComputing optimal skip alignments in normal form probabilities.")
        em = ExecutionManager()
        if self.pn_method == EbiWeights.OCCURANCE or self.pn_method is None:
            ebi = EbiOccurance()
        else:
            raise ValueError("Not implemented Wbi weights method used:", self.pn_method)
        time_start_ns_agns = time.process_time_ns()
        C, global_C = em.coninciding_agns(skip_dict)
        time_stop_ns_agns = time.process_time_ns()
        self.agn_time = ((time_stop_ns_agns-time_start_ns_agns)/sum(len(v) for k,v in skip_dict.items()), time_stop_ns_agns-time_start_ns_agns) # avg, total
        em.validate(global_C)
        self.C = C
        self.global_C = global_C
        var_C = em.coinciding_agns_var(global_C)
        self.var_C = var_C
        activity_to_id = ebi.write_tree_to_petri(self.tree)
        print(activity_to_id)
        if self.pn_log is not None:
            id_filtered_log_rf = ebi.write_log(self.pn_log, activity_to_id)
            ebi.ebi_slpn(out=slpn_path)
            ebi.validate_slpn(self.tree, slpn_path)
            ebi.update_slpn_weights(self.tree, slpn_path)
            time_start_ns_model_paths = time.process_time_ns()
            trace_probs, trace_counts, prob_time = ebi.trace_probs(var_C, model=slpn_path)
            time_stop_ns_model_paths = time.process_time_ns()
            self.trace_prob_time = (prob_time/len(trace_probs), prob_time) # avg, total
        else:
            trace_probs, trace_counts, prob_time = ebi.trace_probs(var_C, self.pn_measure)
            self.trace_prob_time = (-1, -1)
        self.write(trace_probs, path, "trace_probs")
        self.write(trace_counts, path, "trace_counts")
        time_start_ns_derivation1 = time.process_time_ns()
        score_sagn, cond_prob_ = ebi.skip_agn_probs_traversal(self.tree, skip_dict, C, trace_probs, trace_counts)
        self.cond_prob_ = cond_prob_
        time_stop_ns_derivation1 = time.process_time_ns()

        # skips
        print("3/4\tComputing skips from all optimal skip alignments in normal form.")
        s = Skipper()
        s.fix_sagns(self.tree, skip_dict_backup)
        time_start_ns_derivation2 = time.process_time_ns()
        skips = s.conditional_skip_prob(self.tree, skip_dict_backup)
        self.skips = skips
        time_stop_ns_derivation2 = time.process_time_ns()

        # skip probabilities
        print("4/4\tComputing skip probabilities.")
        id_to_skip_dict = {s.name:s for _,sl in skip_dict.items() for s in sl}
        self.id_to_skip_dict = id_to_skip_dict
        id_to_skip_dict_backup = {s.name:s for _,sl in skip_dict_backup.items() for s in sl}
        self.id_to_skip_dict_backup = id_to_skip_dict_backup
        time_start_ns_derivation3 = time.process_time_ns()
        self.skip_probs = self.recur_node_prob(self.tree, self.variants, skip_dict, skips, cond_prob_, self.pl, id_to_skip_dict, id_to_skip_dict_backup, skip_dict_backup)
        time_stop_ns_derivation3 = time.process_time_ns()
        self.derivation_time = (((time_stop_ns_derivation1-time_start_ns_derivation1)+(time_stop_ns_derivation2-time_start_ns_derivation2)+(time_stop_ns_derivation3-time_start_ns_derivation3))/len(self.variants.keys()), (time_stop_ns_derivation1-time_start_ns_derivation1)+(time_stop_ns_derivation2-time_start_ns_derivation2)+(time_stop_ns_derivation3-time_start_ns_derivation3)) # avg, total
        self.write(self.skip_probs, path, "skip_probs")

        print("Done.")
    
    def results(self):
        return self.skip_probs
    
    def stats(self):
        print("---=== Skip alignment computation (timeout = " + str(self.sagn_timeout) + "s) ===---")
        print("Avg. number of sagns per trace variant [incl timeouts]:", sum(len(v) for v in self.skip_dict.values())/len(self.variants))
        print("Total number of sagns [incl timeouts]:", sum(len(v) for v in self.skip_dict.values()))
        print("Avg. time per log trace variant (ns) [no timeouts]:", sum(v for v in self.skip_times.values() if v != -1)/sum(v!=-1 for v in self.skip_times.values()))
        print("Total time all sagns (ns) [no timeouts]:", sum(v for v in self.skip_times.values() if v != -1))

        print("---=== Unfolding skip alignments ===---")
        print("Avg. number of agns for a trace variant [incl timeouts] (ns):", sum(len(v) for v in self.var_C.values())/len(self.var_C))
        print("Total number of agns [incl timeouts] (ns):", sum(len(v) for v in self.var_C.values()))
        print("Avg. time per set of agns for a skip alignment (ns):", self.agn_time[0])
        print("Total time for all agns (ns):", self.agn_time[1])

        print("---=== Ebi calls for model paths ===---")
        print("Avg. time per model path (ns):", self.trace_prob_time[0])
        print("Total time for all model paths (ns):", self.trace_prob_time[1])

        print("---=== Derivation skip probabilities ===---")
        print("Avg. time per log trace variant (ns):", self.derivation_time[0])
        print("Total time all log trace variant (ns):", self.derivation_time[1])
    
    def print_blinded(self, tree:ProcessTree=None):
        if tree is None:
            return self.print_blinded(self.tree)
        if isinstance(tree, LeafNode):
            return tree.__str__() + ", " + ("[ " if tree.get_cheapest_execution(0)[1] else "") + str(self.skip_probs[tree]) + (" ]" if tree.get_cheapest_execution(0)[1] else "")
        else:
            if isinstance(tree, Sequence):
                operator = "→"
            elif isinstance(tree, Xor):
                operator = "×"
            elif isinstance(tree, And):
                operator = "∧"
            elif isinstance(tree, Loop):
                operator = "↺"
            else:
                operator = "UNKNOWN"
            child_string = "\n".join([self.print_blinded(c) for c in tree.children])
            return (" " * tree.get_distance_to_root()*2) + operator + ", " + ("[ " if tree.get_cheapest_execution(0)[1] else "") + str(self.skip_probs[tree]) + (" ]" if tree.get_cheapest_execution(0)[1] else "") + "\n" + child_string

    # setup
    def get_variant_dict(self, log):
        variants = dict()
        for k,v in pm4py.statistics.variants.log.get.get_variants_from_log_trace_idx(log).items():
            variants[k] = len(v)
        variants = dict(sorted(variants.items(), key=lambda x: -x[1]))
        return variants

    def variant_prob(self, variant:tuple, variants:Dict):
        if variant in variants:
            return variants[variant]/sum(variants.values())
        else: return 0.
    
    def variant_prob_dist(self, log):
        var = self.get_variant_dict(log)
        variant_probs = {k:self.variant_prob(k, var) for k,_ in var.items()}
        return variant_probs

    # skip alignments
    def compute_skip_alignments(self, timeout=600):
        skip_dict = {}
        skip_times = {}
        Aligner.set_level_incentive(0)
        variant_strings = [list(var) for var in self.variants.keys()]
        futures = align_sk_all(variant_strings, self.tree, timeout=timeout)
        for index, variant in enumerate(futures):
            agns, t = futures[index].result()
            skip_dict[", ".join(variant_strings[index])] = agns
            skip_times[", ".join(variant_strings[index])] = t
        print("Number of computed optimal skip alignments in normal form:", sum(len(v) for k,v in skip_dict.items()))
        print("Timeouts:", sum(v==-1 for k,v in skip_times.items()))
        return skip_dict, skip_times
    
    def rename_sagns(self, states:Dict[str,List[State]]):
        counter = 0
        for var,statel in states.items():
            for s in statel:
                s.name = counter
                counter += 1
    
    # skip probabilities
    def prob_per_state_and_node(self, node:ProcessTree, state_id:int, variant:tuple[str], skip_prob:Dict[State,Dict[ProcessTree,float]], sagn_prob:Dict[State,float], variant_prob:Dict[tuple[str],float], id_to_skip_dict:Dict[int,State], id_to_skip_dict_backup:Dict[int,State]):
        f1 = skip_prob[id_to_skip_dict_backup[state_id]][node]
        f2 = sagn_prob[node][id_to_skip_dict[state_id]]
        f3 = variant_prob[variant]
        return f1*f2*f3
    
    def prob_per_variant_and_node(self, node:ProcessTree, variant:tuple[str], skip_dict:Dict[str,List[State]], skip_prob:Dict[State,Dict[ProcessTree,float]], sagn_prob:Dict[State,float], variant_prob:Dict[tuple[str],float], id_to_skip_dict:Dict[int,State], id_to_skip_dict_backup:Dict[int,State], skip_dict_backup:Dict[str, List[State]]):
        prob = 0.
        for s in skip_dict[", ".join(variant)]:
            prob += self.prob_per_state_and_node(node, s.name, variant, skip_prob, sagn_prob, variant_prob, id_to_skip_dict, id_to_skip_dict_backup)
        # we need to fix the conditional skip alignment probability for those sagns that never reach the node
        prob_of_sagns_reaching_n = 0.
        for s in skip_dict_backup[", ".join(variant)]:
            prob_of_sagns_reaching_n += Skipper().node_reached(node, s) * sagn_prob[node][id_to_skip_dict[s.name]]
        if prob_of_sagns_reaching_n == 0:
            return 0.
        return prob * 1/prob_of_sagns_reaching_n
    
    def prob_per_node(self, node:ProcessTree, variants:Dict[tuple[str],int], skip_dict:Dict[str,List[State]], skip_prob:Dict[State,Dict[ProcessTree,float]], sagn_prob:Dict[State,float], variant_prob:Dict[tuple[str],float], id_to_skip_dict:Dict[int,State], id_to_skip_dict_backup:Dict[int,State], skip_dict_backup:Dict[str, List[State]]):
        prob = 0.
        for v in variants.keys():
            prob += self.prob_per_variant_and_node(node, v, skip_dict, skip_prob, sagn_prob, variant_prob, id_to_skip_dict, id_to_skip_dict_backup, skip_dict_backup)
        # we need to fix the variant probability for those variants that never reach the node (i.e., in no skip alignment)
        prob_of_traces_reaching_n = 0.
        for v in variants.keys():
            reached = any(Skipper().node_reached(node, s) for s in skip_dict_backup[", ".join(v)])
            prob_of_traces_reaching_n += reached * variant_prob[v]
        if prob_of_traces_reaching_n == 0:
            return 0.
        return prob * 1/prob_of_traces_reaching_n
    
    # skip probabilities
    def recur_node_prob(self, tree:ProcessTree, variants:Dict[tuple[str],int], skip_dict:Dict[str,List[State]], skip_prob:Dict[State,Dict[ProcessTree,float]], sagn_prob:Dict[State,float], variant_prob:Dict[tuple[str],float], id_to_skip_dict:Dict[int,State], id_to_skip_dict_backup:Dict[int,State], skip_dict_backup:Dict[str, List[State]]):
        result = {}
        if isinstance(tree, LeafNode):
            result[tree] = self.prob_per_node(tree, variants, skip_dict, skip_prob, sagn_prob, variant_prob, id_to_skip_dict, id_to_skip_dict_backup, skip_dict_backup)
            return result
        else:
            for c in tree.children:
                res_tmp = self.recur_node_prob(c, variants, skip_dict, skip_prob, sagn_prob, variant_prob, id_to_skip_dict, id_to_skip_dict_backup, skip_dict_backup)
                for k,v in res_tmp.items():
                    result[k] = v
            result[tree] = self.prob_per_node(tree, variants, skip_dict, skip_prob, sagn_prob, variant_prob, id_to_skip_dict, id_to_skip_dict_backup, skip_dict_backup)
            return result