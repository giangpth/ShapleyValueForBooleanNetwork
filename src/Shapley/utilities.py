import argparse as ap
from typing import Dict, Any
import hashlib
import json
import math
import copy 
import random 
import numpy as np
from sklearn.metrics import ndcg_score
from math import isfinite
import random
import math

from Shapley.exceptions import InforError
from Shapley.rankingEvaluation import compare_rankings, print_rank_comparison
from Shapley.visualization import showNetwork

def tie_buckets_abs(values, atol=1e-4, rtol=1e-4, rep="max"):
    """
    values: dict[key] -> float (will use absolute value)
    returns: dict[key] -> (bucket_id, bucket_rep)
    """
    items = [(k, abs(v)) for k, v in values.items() if isfinite(v)]
    items.sort(key=lambda kv: kv[1], reverse=True)

    def close(a, b):
        return abs(a-b) <= max(atol, rtol*max(abs(a), abs(b)))

    buckets, current = [], []
    for k, val in items:
        if not current:
            current = [(k, val)]
        else:
            # compare to representative of current bucket
            ref = (max(x[1] for x in current) if rep=="max"
                   else (sum(x[1] for x in current)/len(current)))
            if close(val, ref):
                current.append((k, val))
            else:
                buckets.append(current)
                current = [(k, val)]
    if current:
        buckets.append(current)

    # assign bucket ids and representatives
    out = {}
    for bid, bucket in enumerate(buckets, start=1):
        rep_val = (max(x[1] for x in bucket) if rep=="max"
                   else (sum(x[1] for x in bucket)/len(bucket)))
        for k, _ in bucket:
            out[k] = (bid, rep_val)
    return out

def dense_ranks_from_buckets(bucket_map):
    # bucket_id already encodes order (1 = top)
    return {k: bucket_map[k][0] for k in bucket_map}

# Example: build ranks for baseline and propagated, then use your tie-aware metrics
# base_ranks = dense_ranks_from_buckets(tie_buckets_abs(baseline, atol=1e-6, rtol=1e-6))
# prop_ranks = dense_ranks_from_buckets(tie_buckets_abs(propagated, atol=1e-6, rtol=1e-6))


def rank_dict_values(d, method="dense", tol=None, largest=True):
    """
    Rank by absolute value.
      - largest=True -> DESC (largest magnitude rank 1)
      - largest=False -> ASC  (smallest magnitude rank 1)
    method: 'dense', 'min' (competition), 'max' (modified competition),
            'average' (fractional), 'ordinal'
    tol: floats within tol are treated as equal
    """
    if not d:
        return {}

    items = [(k, abs(v)) for k, v in d.items()]
    items.sort(key=lambda kv: kv[1], reverse=largest)  # DESC if largest=True

    groups = []
    gid = []
    prev = None
    g = -1
    for _, val in items:
        same = (val == prev) if tol is None else (prev is not None and abs(val - prev) <= tol)
        if not same:
            g += 1
            prev = val
        gid.append(g)

    first_pos, last_pos = {}, {}
    for i, g in enumerate(gid, start=1):
        first_pos.setdefault(g, i)
        last_pos[g] = i

    out = {}
    for i, (k, _) in enumerate(items, start=1):
        g = gid[i-1]
        if method == "dense":
            r = g + 1
        elif method == "min":        # all ties get first position in the tie block
            r = first_pos[g]
        elif method == "max":        # all ties get last position in the tie block
            r = last_pos[g]
        elif method == "average":    # mean of first..last positions
            r = (first_pos[g] + last_pos[g]) / 2
        elif method == "ordinal":    # break ties by order in 'items'
            r = i
        else:
            raise ValueError("Unknown method")
        out[k] = r
    return out


def dict_hash(dictionary: Dict[str, Any]) -> str:
    """MD5 hash of a dictionary."""
    dhash = hashlib.md5()
    # We need to sort arguments so {'a': 1, 'b': 2} is
    # the same as {'b': 2, 'a': 1}
    encoded = json.dumps(dictionary, sort_keys=True).encode()
    dhash.update(encoded)
    return dhash.hexdigest() 

def print_ranked_keys(d, orispecies=None):
    if orispecies:
        subd = {k: v for k, v in d.items() if k in orispecies}
    else:
        subd = d

    # sort by value descending
    sorted_items = sorted(subd.items(), key=lambda x: abs(x[1]), reverse=True)

    ranks = []
    current_rank = 1
    last_value = None

    # build rank groups
    for key, value in sorted_items:
        if value != last_value:
            ranks.append((current_rank, [key]))
            current_rank += 1
            last_value = value
        else:
            ranks[-1][1].append(key)

    # print with max 4 keys per row
    MAX_PER_ROW = 6

    for rank, keys in ranks:
        # slice keys into groups of size 4
        for idx in range(0, len(keys), MAX_PER_ROW):
            chunk = keys[idx:idx + MAX_PER_ROW]
            rank_label = str(rank) if idx == 0 else ""   # only show rank on first line
            print(f"{rank_label}\t{', '.join(chunk)}")





def parseArgumentsAnalysis():
    """
    Parse command line arguments.

    Returns:
        argparse.ArgumentParser: The argument parser with defined options. 
    """
    
    parser = ap.ArgumentParser(description='Input the expression file and ')
    parser.add_argument('-n', '--name', help="Name of the system", required=True)
    parser.add_argument('-e', '--expression', type=str, help="path to expression file", required=True)
    parser.add_argument('-o', '--output', type=str, help='list of interested components to examine', required=True)
    parser.add_argument('-r', '--sampling', help='Sampling strategy for input states: "all" (default), "uniform", "antithetic", or "stratified"', type=str, default='all')
    parser.add_argument('-p', '--samplePer', help='Percentage of combinations to sample if using a sampling strategy different from "all". Default is 60.0', type=float, default=60)
    parser.add_argument('-c', '--simPer', help='Percentage of nodes to be simulated. Default is 0.0', type=float, default=0)
    parser.add_argument('-k', '--scheme', help='Scheme of calculation, can be "Shapley" or "Uniform", default is "Shapley" ', type=str, default='Shapley')
    parser.add_argument('-m', '--mode', help='Mode of calculation, can be "Test" or "Propagation" or "Simulation", default is "Propagation" ', type=str, default='Propagation')
    parser.add_argument('-s', '--simprop', help='Probability to simulate an original node during the propagation', type = float, default=0.3)
    parser.add_argument('-v', '--verbose', help="Print log or not", action='store_true')
    parser.add_argument('-d', '--decimal', help="Log decimal mappings for hardcore debug mode", action='store_true')
    parser.add_argument('-t', '--test', help="Perform full network perturbation to test with the obtained result", action='store_true')
    
    
    return parser 


def decimal_pairs(sortedinput, target):
    n = len(sortedinput)
    k = n - 1 - sortedinput.index(target)   # target's bit position
    bit = 1 << k
    for r in range(1 << (n - 1)):
        low = r & (bit - 1)                 # bits below the target
        high = (r >> k) << (k + 1)          # bits above the target, shifted up
        off = high | low
        yield off | bit, off                # (target=1, target=0)

def genRowids2simulate(numinputs, errthreshold, mode='Shapley', debug=False):
    """
    Generate a set of row IDs to simulate based on the total number of rows and an error threshold.
    Args:
        numinput (int): The number of input nodes.
        errthreshold (float): The error threshold (between 0 and 1).
        mode (str): The mode of operation, can be 'Shapley' or 'Uniform', default is 'Shapley'.
        debug (bool): If True, print debug information.
    Returns:
        set: A set of row IDs to simulate.
    """
    if mode == 'Uniform':
        totalrows = pow(2, numinputs)
        simrows = math.ceil((1 - errthreshold) * totalrows)
        if debug:
            print(f"Total rows: {totalrows}, need to simulate at least {simrows} rows according to the error threshold {errthreshold}")
        rowidstosimulate = set()
        while len(rowidstosimulate) < simrows:
            randint = random.randint(0, totalrows - 1) 
            rowidstosimulate.add(randint)
    else: # Shapley mode
        maxval = numinputs + 1
        maxerr = maxval*(1.0-errthreshold)
        # simulate row with highest prop first, stop adding rows when the error is below threshold
        rowidstosimulate = set()
        cumerr = 0.0
        for size in range(0, int(numinputs/2) + 1): 
            numrowswithsize = math.comb(numinputs, size)
            perrowerr = 1.0 / numrowswithsize
            for i in range(numrowswithsize):
                bitstr = '1'*size + '0'*(numinputs - size)
                bitlist = list(bitstr)
                # permute the bitlist to get the ith combination
                from itertools import permutations
                perm = set(permutations(bitlist))
                perm = sorted(list(perm))
                selectedbitstr = ''.join(perm[i])
                rowid = int(selectedbitstr, 2)
                rowidstosimulate.add(rowid)
                cumerr += perrowerr
                if cumerr >= maxerr:
                    break
            if cumerr >= maxerr/2:
                break
        for size in range(numinputs, int(numinputs/2), -1): 
            numrowswithsize = math.comb(numinputs, size)
            perrowerr = 1.0 / numrowswithsize
            for i in range(numrowswithsize):
                bitstr = '1'*size + '0'*(numinputs - size)
                bitlist = list(bitstr)
                # permute the bitlist to get the ith combination
                from itertools import permutations
                perm = set(permutations(bitlist))
                perm = sorted(list(perm))
                selectedbitstr = ''.join(perm[i])
                rowid = int(selectedbitstr, 2)
                rowidstosimulate.add(rowid)
                cumerr += perrowerr
                if cumerr >= maxerr:
                    break
            if cumerr >= maxerr:
                break
    if debug:
        print(f"Row IDs to simulate: {sorted(list(rowidstosimulate))}") 
    return rowidstosimulate


def readfile(path, debug=False):
    """
    Read a file and return its lines as a list.
    Args:
        path (str): The path to the file.
        debug (bool): If True, print each line as it is read.
    Returns:
        list: A list of lines from the file, stripped of leading/trailing whitespace.
    Raises:
        IOError: If the file cannot be opened.
    """
    try:
        f = open(path, 'r')
        lines = []
        for line in f:
            line = line.strip()
            if debug:
                print(line)
            lines.append(line)
        return lines
    except IOError:
        print("Error opening file")

def top_rank_keys(d, top_n=10):
    # get unique ranks sorted ascending (best rank first)
    ranks = sorted(set(d.values()))
    # take up to top_n distinct ranks
    cutoff_ranks = set(ranks[:top_n])
    # return all keys whose rank is in these cutoff ranks
    return {k for k, v in d.items() if v in cutoff_ranks}

def top_matches(d1, d2, top_n=10):
    top1 = top_rank_keys(d1, top_n)
    top2 = top_rank_keys(d2, top_n)
    print(top1.difference(top2))
    print(top2.difference(top1))
    return len(top1 & top2), top1 & top2, len(top1)



def getOrderedList(inputnames, internames, debug=False):
    """
    Get sorted lists of input and intermediate species names.
    Args:
        inputnames (set): A set of input species names.
        internames (set): A set of intermediate species names.
        debug (bool): If True, print the sorted lists and their lengths.
    Returns:
        tuple: A tuple containing two lists - sorted input names and sorted intermediate names.
    """
    sortedinput = sorted(list(inputnames))
    if debug:
        print('Sorted input names')
        print(sortedinput)
        print("There are {} inputs with the range of input: {} {}".format(len(sortedinput),0,2**(len(sortedinput)) - 1))
    sortedinter = sorted(list(internames))
    if debug:
        print('Sorted internames:')
        print(sortedinter)
        print("There are {} intermediate nodes with the range of inter: {} {}".format(len(sortedinter), 0,2**(len(sortedinter)) - 1))
    return sortedinput, sortedinter 


def toDecimal(state, sortedinput, sortedinter=[], debug=False):
    def bits(names):
        return ''.join('1' if state[n] is True else '0' for n in names)
    inbits, interbits = bits(sortedinput), bits(sortedinter)
    if debug:
        print(inbits, interbits)
    return int(inbits, 2), (int(interbits, 2) if interbits else 0)


def merge2states(state1, state2, debug=False):
    """
    Merge two states by performing a logical OR operation on their values.
    For final state of periodict attractor, if a node is True in any of the states, it will be True in the merged state.
    Args:
        state1 (dict): The first state dictionary.
        state2 (dict): The second state dictionary.
        debug (bool): If True, print debug information.
    Returns:
        dict: The merged state dictionary.
    """
    blinking = set()
    for item in state1.keys():
        # print(state1[item], state2[item])
        state1[item] = state1[item] or state2[item] 
        # state1[item] = state1[item] and state2[item]
        if state1[item] ^ state2[item]:
            blinking.add(item)
    return state1, blinking 


def subsets(s):  
    """ Return all subsets of a set  
    Args:
        s (list): A list of elements.
    Returns:
        list: A list of lists, each representing a subset of the input list.
    """
    if len(s) == 0:  
        return [[]]  
    x = subsets(s[:-1])  
    return x + [[s[-1]] + y for y in x] 

def processPrecomputSimulations(intactbitable, output, rowids, KOres, KIres, debug=False):
    """
    Process precomputed simulation results for intact, knockout, and knockin scenarios.
    Args:
        intactbitable (dict): A dictionary mapping input states to output states for the intact network.
        output (str): The output species name to consider.
        rowids (set): A set of row IDs to consider.
        KOres (dict): A dictionary mapping intermediate nodes to their knockout simulation results.
        KIres (dict): A dictionary mapping intermediate nodes to their knockin simulation results.
        debug (bool): If True, print debug information.
    Returns:
        dict: A dictionary mapping rows accounted for intermediate nodes.
    """
    # print(intactbitable)
    # print(KOres)
    # print(KIres)
    precountedrow = dict()
    for internode, koresults in KOres.items():
        # print (koresults)
        # print (KIres[internode])

        precountedrow[internode] = set()
        if internode not in KIres:
            print("Cannot find knockin results for node {}".format(internode))

        for rowid in rowids:
            intactresult = intactbitable[rowid]
            if rowid in koresults:
                koresult = koresults[rowid]
            else:
                print("Cannot find knockout result for node {} at row {}".format(internode, rowid))
                continue
            if rowid in KIres[internode]:
                kiresult = KIres[internode][rowid]
            else:
                print("Cannot find knockin result for node {} at row {}".format(internode, rowid))
                continue

            # print("Node {}: intact output {}, KO output {}, KI output {} at row {}".format(internode, intactresult[output], koresult[output], kiresult[output], rowid))
            if intactresult[output] != koresult[output] or intactresult[output] != kiresult[output]:
                precountedrow[internode].add(rowid)
                
    return precountedrow

def compareOutputs(orioutput, outputnode, inputnames, kooutput = None, kioutput = None, debug=False):
    korow, kirow = set(), set() 
    # kovalue, kivalue = 0.0, 0.0 
    sortedinputs = sorted(list(inputnames))
    if kooutput or kioutput:
        if kooutput:
            assert len(orioutput) == len(kooutput)
            for rowid, oristate in enumerate(orioutput):
                bitlist = [0 if not oristate[input] else 1 for input in sortedinputs]
                inputid = int(''.join(map(str, bitlist)), 2)

                kostate = kooutput[rowid]
                # print(kostate)
                if kostate[outputnode] != oristate[outputnode]:
                    # if kostate[pernode] == oristate[pernode]:
                    #     print("SUSPICOUS situation: no perturbation but output changes")
                    # else:
                    korow.add(inputid)
                else:
                    if debug:
                        print(inputid)
        if kioutput:
            assert len(orioutput) == len(kioutput)
            for rowid, oristate in enumerate(orioutput):
                bitlist = [0 if not oristate[input] else 1 for input in sortedinputs]
                inputid = int(''.join(map(str, bitlist)), 2)
                kistate = kioutput[rowid]
                if kistate[outputnode] != oristate[outputnode]:
                    # if kistate[pernode] == oristate[pernode]:
                        
                    #     print("SUSPICOUS situation: no perturbation but output changes")
                    # else:
                    kirow.add(inputid)
    else:
        print("No perturbed table provided")

    return korow, kirow
    

def reportResults(kovalues, kivalues, inputnames, internames, outputnames, korows = None, kirows = None, debug = False):
    print(f"-------RESULT for output {outputnames}--------")
    if korows and kirows:
        print("Ordered input nodes:")
        print(sorted(inputnames))
        for internode in inputnames + internames: 
            print(f"----{internode}-----")
            print("KO\t{}: {}".format(len(korows[internode]), sorted(list(korows[internode]))))
            print("KI\t{}: {}".format(len(kirows[internode]), sorted(list(kirows[internode]))))
            print("All\t{}: {}".format(len(set(korows[internode]).union(set((kirows[internode])))), \
                                      sorted(list(set(korows[internode]).union(set((kirows[internode])))))))

    for internode in inputnames + internames:
            kovalues[internode] = round(kovalues[internode], 4)
            kivalues[internode] = round(kivalues[internode], 4)
            print("{:20} \t\t\t KO: {:10} | KI: {:10}".format(internode, kovalues[internode], kivalues[internode]))

    orispecies = kovalues.keys()
    print("-----KO RANKING-----")
    print_ranked_keys(kovalues, orispecies)
    print("-----KI RANKING-----")
    print_ranked_keys(kivalues, orispecies)
        

def genTableFromOutput(simoutputs, inputnames, blinkings, sortedinter, outputname, scheme='Shapley', debug=False):
    """
    Generate a table from simulation outputs, indexing rows by species states.
    Args:
        simoutputs (list): A list of dictionaries, each representing an output state with species names as keys and their boolean states as values.
        inputnames (set): A set of input species names.
        sortedinput (list): A sorted list of input species names.
        sortedinter (list): A sorted list of intermediate species names.
        outputnames (set): A set of output species names.
        scheme (str): The scheme of calculation, can be 'Shapley' or 'Uniform', default is 'Shapley'.
        debug (bool): If True, print debug information.
    Returns:
        dict: A dictionary mapping row IDs to dictionaries of species states.
        dict: A dictionary mapping species names to sets of row IDs where they are True.
        dict: A dictionary mapping species names to sets of row IDs where they are False.
    """
    print("----Generate table from output----")
    # if debug:
    #     print("Blinkings:")
    #     for inputid, bset in enumerate(blinkings):
    #         print(f"Row {inputid}: blinking nodes: {bset}")
    index = dict() # for each node, save the IDs of row that the node is TRUE
    aindex = dict() # for each node, save the IDs of row that the node is TRUE
    bindex = dict() # for each node, save the IDs of row that the node is BLINKING

    dictresult = dict()
    sortedintput = sorted(list(inputnames))
    for inter in sortedinter:
        index[inter] = set()
        aindex[inter] = set()
        bindex[inter] = set()
    for input in inputnames:
        index[input] = set()
        aindex[input] = set()
        bindex[input] = set()

    index[outputname] = set() 
    aindex[outputname] = set() 
    bindex[outputname] = set()
    
    # print("Integer verson of output:")
    for rowid, line in enumerate(simoutputs):
        bitlist = [0 if not line[input] else 1 for input in sortedintput]
        inputid = int(''.join(map(str, bitlist)), 2)
        # inp, inter = toDecimal(line, sortedinput, sortedinter)
        # print(inp, inter)
        # print(line)
        # each line is a dictionary with keys are name of species and value is true or false
        
        # for test 
        # if inputid in [0, 1, 3, 4, 5, 6, 7]:
        #     print(f"-------TESTING LINE: {inputid}-------")
        #     print(line)
        # print(line)
        size = 0
        for input in inputnames:
            if line[input]:
                size += 1
        line['SIZE'] = size
        if scheme == 'Uniform':
            line['PROP'] = round(1.0 / pow(2, len(inputnames)),4)
            # line['PROP'] = round(math.factorial(size)*math.factorial(len(inputnames) - size)/math.factorial(len(inputnames)),4)
        else: # Shapley mode
            line['PROP'] = round(math.factorial(size)*math.factorial(len(inputnames) - size)/math.factorial(len(inputnames)),4)
        blinking = blinkings[rowid] # this is the set of nodes that blink in this row 
        for inter in sortedinter:
            if line[inter]:
                if inter in blinking:
                    bindex[inter].add(inputid)
                index[inter].add(inputid)
            else:
                aindex[inter].add(inputid)
        for input in inputnames:
            if line[input]:
                if input in blinking:
                    bindex[input].add(inputid)
                index[input].add(inputid)
            else:
                aindex[input].add(inputid)
       
        if line[outputname]:
            if outputname in blinking:
                bindex[outputname].add(inputid)
            index[outputname].add(inputid)
        else:
            aindex[outputname].add(inputid)

        dictresult[inputid] = line
    return dictresult, index, aindex, bindex

def filterrows(op, nodetofilter, childnode, index, aindex, extranodes, allrows):
    """
    Filter rows based on the state of a specific node and a logical operator.
    Args:
        op (str): The logical operator ('AND' or 'OR').
        nodetofilter (str): The name of the node to filter on.
        onenode (str): The name of the node being evaluated.
        index (dict): A dictionary mapping node names to sets of row IDs where they are True.
        aindex (dict): A dictionary mapping node names to sets of row IDs where they are False.
        extranodes (list): A list of extra nodes to be considered in the filtering.
        allrows (set): A set of all row IDs to be filtered.
    Returns:
        set: A set of filtered row IDs.
    """
    print("FILTER")
    if op == 'OR':
        
        if nodetofilter not in extranodes:
            # print("NOT EXTRANODE")
            childrows = aindex[nodetofilter]
            # countedchild[onenode] = childrows
            print("Counting the rows that {} is FALSE for {} in operator {}".format(nodetofilter, childnode, op))
            allrows = allrows.intersection(childrows)
            print("\tRows after filter are {}".format(sorted(list(allrows))))
        else:
            # print("EXTRANODE")
            coms = nodetofilter.split("_to_")
            if len(coms) == 2 and coms[0] == coms[1]:
                print("{} is selfloop, do not count it".format(nodetofilter))
            else:
                childrows = aindex[nodetofilter]
                # countedchild[onenode] = childrows
                print("Counting the rows that {} is FALSE for {} in operator {}".format(nodetofilter, childnode, op))
                allrows = allrows.intersection(childrows)
                print("\tRows after filter are {}".format(sorted(list(allrows))))
    elif op == 'AND':
        if nodetofilter not in extranodes:
            # print("NOT EXTRANODE")
            childrows = index[nodetofilter]
            # countedchild[onenode] = childrows
            print("Counting the rows that {} is TRUE for {} in operator {}".format(nodetofilter, childnode, op))
            allrows = allrows.intersection(childrows)
            print("\tRows after filter are {}".format(sorted(list(allrows))))
        else:
            # print("EXTRANODE")
            coms = nodetofilter.split("_to_")
            if len(coms) == 2 and coms[0] == coms[1]:
                print("{} is selfloop, do not count it".format(nodetofilter))
            else:
                childrows = index[nodetofilter]
                # countedchild[onenode] = childrows
                print("Counting the rows that {} is TRUE for {} in operator {}".format(nodetofilter, childnode, op))
                allrows = allrows.intersection(childrows)
                print("Rows after filter are {}".format(sorted(list(allrows))))
        # pass
    else:
        print("Do not support {} operator".format(op))
    
    return allrows


def random_percentage_selection(elements, percentage):
    """
    Randomly selects approximately 'percentage'% of elements from the given iterable.
    
    Args:
        elements (list or set): The collection of elements to choose from.
        percentage (float): The percentage to select (e.g., 30 for 30%).
    
    Returns:
        list: A list of randomly selected elements.
    """
    if not elements:
        return []
    
    # Convert set to list if needed (random.sample requires a sequence)
    elem_list = list(elements)
    
    # Calculate the number of elements to select
    # Use math.ceil for rounding up, or int() for truncation, or round() for nearest
    num_to_select = max(1, math.ceil(len(elem_list) * (percentage / 100.0)))
    
    # Randomly sample without replacement
    selected = random.sample(elem_list, num_to_select)
    
    return selected


def genAllInputForSHAP(numfeature):
    """
    Generate all possible binary input combinations for a given number of features.
    Args:
        numfeature (int): The number of binary features.
    Returns:
        np.ndarray: A 2D array where each row represents a unique binary input combination.
    """
    inputs = []
    for i in range(pow(2,numfeature)):
        code = '{0:0' + str(numfeature) + 'b}'
        bistr = code.format(i)
        # print(bistr)
        # oneinput = np.array(list(bistr))
        inputs.append(bistr)
    
    binary_array = np.array([list(map(int, list(s))) for s in inputs])
    # print(binary_array.shape)

    return binary_array
    # return inputs 

def identity_masker(data):
    return data

def genRandomInputForSHAP(numfeature, numsample):
    """
    Generate random binary input combinations for a given number of features.
    Args:
        numfeature (int): The number of binary features.
        numsample (int): The number of random samples to generate.
    Returns:
        np.ndarray: A 2D array where each row represents a unique binary input combination.
    """
    inputs = []
    for i in range(numsample):
        randint = random.randint(0, pow(2, numfeature) - 1) 
        code = '{0:0' + str(numfeature) + 'b}'
        bistr = code.format(randint)
        # print(bistr)
        # oneinput = np.array(list(bistr))
        inputs.append(bistr)
    
    binary_array = np.array([list(map(int, list(s))) for s in inputs])
    # print(binary_array.shape)

    return binary_array
    # return inputs 


def evaluation2(calkorows, calkirows, calkovalues, calkivalues, gtkorows, gtkirows, gtkovalues, gtkivalues, inputnames, outname, debug=False):
    wrongrows = 0
    examinednodes = 0
    for node, rows in calkorows.items():
        print(f"----Report rows for node {node} ------")
        korow = sorted(list(calkorows.get(node, [])))
        kirow = sorted(list(calkirows.get(node, [])))
        gtkorow = sorted(list(gtkorows.get(node, [])))
        gtkirow = sorted(list(gtkirows.get(node, [])))
        kolack = sorted(list(set(gtkorow).difference(set(korow))))
        koextra = sorted(list(set(korow).difference(set(gtkorow)))) 
        print("KO Lack {:5}: {}".format(len(kolack), sorted(list(kolack))))
        print("KO Extr {:5}: {}".format(len(koextra), sorted(list(koextra))))
        print("KO Real {:5}: {}".format(len(gtkorow), sorted(list(gtkorow))))
        print("KO Calc {:5}: {}\n".format(len(korow), sorted(list(korow))))
        kilack = sorted(list(set(gtkirow).difference(set(kirow))))
        kiextra = sorted(list(set(kirow).difference(set(gtkirow)))) 
        print("KI Lack {:5}: {}".format(len(kilack), sorted(list(kilack))))
        print("KI Extr {:5}: {}".format(len(kiextra), sorted(list(kiextra))))
        print("KI Real {:5}: {}".format(len(gtkirow), sorted(list(gtkirow))))
        print("KI Calc {:5}: {}\n".format(len(kirow), sorted(list(kirow))))

        wrongrows += len(kolack) + len(koextra) + len(kilack) + len(kiextra)
        examinednodes += 1

    errorko = 0.0
    errorki = 0.0 
    num = 0
    for node, value in calkovalues.items():
        kovalue = round(calkovalues[node], 4)
        kivalue = round(calkivalues[node], 4)
        
        gtkovalue = round(gtkovalues[node], 4)
        gtkivalue = round(gtkivalues[node], 4)
        num += 1

        if abs(gtkovalue - kovalue) >= 0.005 or abs(gtkivalue - kivalue) >= 0.005:
            print("{:20} - Correct: KO: {:10} | KI : {:10} - Incorrect KO: {:10} | KI: {:10}".format(node, gtkovalue, gtkivalue, kovalue, kivalue))
        errorko += (gtkovalue - kovalue)*(gtkovalue - kovalue)
        errorki += (gtkivalue - kivalue)*(gtkivalue - kivalue)

    keys = sorted(gtkovalues.keys()) 

    y_true_ko = np.array([[abs(gtkovalues[k]) for k in keys]])
    y_pred_ko = np.array([[abs(calkovalues[k]) for k in keys]])

    y_true_ki = np.array([[abs(gtkivalues[k]) for k in keys]])
    y_pred_ki = np.array([[abs(calkivalues[k]) for k in keys]]) 

    dncg_ko = round(ndcg_score(y_true_ko, y_pred_ko),4)
    dncg_ki = round(ndcg_score(y_true_ki, y_pred_ki),4)
    ave_dncg = (dncg_ko + dncg_ki)/2

        

    orispecies = gtkovalues.keys()

    print("-----KO RANKING COMPARISON-----")
    print_ranked_keys(calkovalues, orispecies)
    print("REFERENCE:")
    print_ranked_keys(gtkovalues)
    print("=======")
    print("-----KI RANKING COMPARISON-----")
    print_ranked_keys(calkivalues, orispecies)
    print("REFERENCE:")
    print_ranked_keys(gtkivalues)
    print("============")

    print("Relative_RMSE_KO_is ", (math.sqrt(errorko/num)/(2.0*((len(inputnames))+1))))
    print("Relative_RMSE_KI_is ", (math.sqrt(errorki/num))/(2.0*((len(inputnames))+1)))
    print("NDCG_KO_is ", dncg_ko)
    print("NDCG_KI_is ", dncg_ki)

    korankingres = compare_rankings(gtkovalues, calkovalues)
    kirankingres = compare_rankings(gtkivalues, calkivalues)
    print_rank_comparison(korankingres, kirankingres)

    avewrongrows = wrongrows/examinednodes 
    tabsize = 2**len(inputnames)
    averightrows = ((tabsize) - avewrongrows)/((tabsize))
    print("AVERAGE_CORRECT_ROWS_ALL_NODES: ", round(averightrows,4))
    return averightrows, ave_dncg



def evaluation1(proprows, simrows, simtable, outname):
    numinput = math.log2(len(simtable))

    propko, propki = rowstovalues(proprows, simtable, outname)
    simko, simki = rowstovalues(simrows, simtable, outname) 

    # print(f"propko {propko}")
    # print(f"propki {propki}")
    # print(f"simko {simko}")
    # print(f"simki {simki}")


    unreachablekeys = simrows.keys() - proprows.keys() 
    print(f"There are {len(unreachablekeys)} unreachable nodes, they are:")
    print(unreachablekeys)

    errorko = 0.0
    errorki = 0.0 
    num = 0

    for node, value in simko.items():
        num += 1
        if node not in propko:
            propko[node] = 0.0
        if node not in propki:
            propki[node] = 0.0 

        if abs(value - propko[node]) >= 0.005 or abs(simki[node] - propki[node]) >= 0.005:
            print("{:20} - Correct: KO: {:10} | KI : {:10} - Incorrect KO: {:10} | KI: {:10}".format(node, value, simki[node], propko[node], propki[node]))
        errorko += (value - propko[node])*(value - propko[node])
        errorki += (simki[node] - propki[node])*(simki[node] - propki[node])

    wrongrows = 0
    examinednodes = 0

    keys = sorted(simko.keys()) 

    y_true_ko = np.array([[abs(simko[k]) for k in keys]])
    y_pred_ko = np.array([[abs(propko[k]) for k in keys]])

    y_true_ki = np.array([[abs(simki[k]) for k in keys]])
    y_pred_ki = np.array([[abs(propki[k]) for k in keys]]) 

    dncg_ko = round(ndcg_score(y_true_ko, y_pred_ko),4)
    dncg_ki = round(ndcg_score(y_true_ki, y_pred_ki),4)
    ave_dncg = (dncg_ko + dncg_ki)/2

    for node in sorted(list(simrows.keys())):
        if node not in proprows:
            proprows[node] = set()
        lack = simrows[node].difference(proprows[node]) 
        extra = proprows[node].difference(simrows[node])
        print(node)
        print("Lack  {:5}: {}".format(len(lack), sorted(list(lack))))
        print("Extra {:5}: {}".format(len(extra), sorted(list(extra))))
        print("Real  {:5}: {}".format(len(simrows[node]), sorted(list(simrows[node]))))
        print("Prop  {:5}: {}\n".format(len(proprows[node]), sorted(list(proprows[node]))))
        wrongrows += len(lack) + len(extra)
        examinednodes += 1

    orispecies = simko.keys()

    print("-----KO RANKING COMPARISON-----")
    print_ranked_keys(propko, orispecies)
    print("REFERENCE:")
    print_ranked_keys(simko)
    print("=======")
    print("-----KI RANKING COMPARISON-----")
    print_ranked_keys(propki, orispecies)
    print("REFERENCE:")
    print_ranked_keys(simki)
    print("============")

    print("Relative_RMSE_KO_is ", (math.sqrt(errorko/num)/(2.0*((numinput)+1))))
    print("Relative_RMSE_KI_is ", (math.sqrt(errorki/num))/(2.0*((numinput)+1)))
    print("NDCG_KO_is ", dncg_ko)
    print("NDCG_KI_is ", dncg_ki)

    korankingres = compare_rankings(simko, propko)
    kirankingres = compare_rankings(simki, propki)
    print_rank_comparison(korankingres, kirankingres)

    avewrongrows = wrongrows/examinednodes 
    averightrows = (len(simtable) - avewrongrows)/(len(simtable))
    print("AVERAGE_CORRECT_ROWS_ALL_NODES: ", round(averightrows,4))
    return averightrows, ave_dncg
    

def estimateErrorFromInputs(siminrows, proprows, nodes):
    average_wrong_rows = 0
    for inname in nodes:
        wrongrow = 0
        if inname in proprows and inname in siminrows:
            # print("Comparing node {}".format(inname))
            # print("Simulated rows: {}".format(sorted(list(siminrows[inname]))))
            # print("Proposed rows: {}".format(sorted(list(proprows[inname]))))
            
            lack = set(siminrows[inname]).difference(set(proprows[inname]))
            extra = set(proprows[inname]).difference(set(siminrows[inname]))
            wrongrow += len(lack) + len(extra)
            # print("{:20}: lack: {} - extra: {}".format(inname, lack, extra))
        else: 
            print("Node {} is not considerred".format(inname))
        if len(siminrows[inname]) > 0:
            average_wrong_rows += wrongrow/len(siminrows[inname])
        else:
            if wrongrow > 0:
                average_wrong_rows += 1
        
    average_wrong_rows = (average_wrong_rows/len(nodes))/(2**len(nodes))
    return average_wrong_rows 


def rowstovalues(rowdict, simtable, outname):
    """
    Calculate the KO and KI values for each node based on the provided row dictionary and simulation table.
    Args:
        rowdict (dict): A dictionary mapping node names to sets of row IDs where they are True.
        simtable (dict): A dictionary mapping row IDs to dictionaries of species states.
        outname (str): The name of the output species to consider.
    Returns: 
        dict: A dictionary mapping node names to their KO values.
        dict: A dictionary mapping node names to their KI values.
    """
    # print(rowdict)
    kovalues = dict()
    kivalues = dict()
    for node, rowids in rowdict.items():
        kovalues[node] = 0.0
        kivalues[node] = 0.0
        for inputid in rowids:
            row = simtable[inputid]
            if row[outname]:
                if row[node]:
                    kovalues[node] += row['PROP']
                else:
                    kivalues[node] += row['PROP']
            else:
                if row[node]:
                    kovalues[node] -= row['PROP']
                else:
                    kivalues[node] -= row ['PROP']
        kovalues[node] = round(kovalues[node], 4)
        kivalues[node] = round(kivalues[node], 4)
    return kovalues, kivalues

def topDown(net, formulas, burows, inrows, node_layers):   
    """
    Top-down approach to refine the rows associated with each node in a directed acyclic graph (DAG).  
    Args:
        net (nx.DiGraph): The directed graph representing the network.
        formulas (dict): A dictionary with species names as keys and their boolean formulas as values.
        burows (dict): A dictionary mapping node names to sets of row IDs that are upper bounds for those nodes.
        inrows (dict): A dictionary mapping input node names to sets of row IDs that are directly associated with those input nodes.
        node_layers (dict): A dictionary mapping node names to their respective layers in the DAG.
    Returns:
        dict: A dictionary mapping node names to their refined sets of row IDs. 
    """
    rowsofnodes = copy.deepcopy(burows) # node: rows that are counted by this node, final! 
    for node, layer in node_layers.items():
        # get rows from inrows:
        try:
            rowsfrominput = inrows[node]
        except:
            # print(f"Cannot find {node} in list of rows for input nodes")
            try:
                rowsfrominput = rowsofnodes[node]
            except:
                rowsofnodes[node] = set()
                rowsfrominput = rowsofnodes[node]
        if layer == 0: # input rows is ground truth
            rowsofnodes[node] = rowsfrominput 
        else:
            try:
                form = formulas[node]
            except:
                print("Pass the input node {}".format(node))
                continue       
            op = form.val
            if op == 'OR' or op == "AND":
                left, right = form.left.val, form.right.val
                try:
                    leftinrows = rowsofnodes[left]
                except:
                    print("Cannot find node {} in rows of input nodes".format(left))
                try:
                    rightinrows = rowsofnodes[right]
                except:
                    print("Cannot find node {} in rows of input nodes".format(right))
                sumrows = leftinrows.union(rightinrows)
                # sumrows is the upper bound of the rows that node can have 
                rowsofnodes[node] = sumrows.intersection(rowsofnodes[node])
            else:
                print("Unary operator, take the rows of incoming node")
                if op == 'NOT':
                    parent = form.right.val
                else:
                    parent = form.val
                outedges = list(net.out_edges(parent))
                if len(outedges) > 1:
                    parentinrows = rowsofnodes[parent]
                    rowsofnodes[node] = parentinrows.intersection(rowsofnodes[node]) 
    return rowsofnodes 

def mergeUnrollNodes(unrollednodes, rowsofnodes): 
    """
    Merge the rows of unrolled nodes back to their original nodes.
    Args:
        unrollednodes (dict): A dictionary mapping unrolled node names to their original node names.
        rowsofnodes (dict): A dictionary mapping node names to sets of row IDs associated with those nodes.
    Returns:
        dict: A dictionary mapping original node names to their merged sets of row IDs.
    """
    print("\n\n-------MERGING unrolled nodes-------")
    finalrowsofnodes = dict() 
    for node, rows in rowsofnodes.items():
        print(f"Node {node} has rows \n{sorted(list(rows))}")
        originalnode = node
        while originalnode in unrollednodes:
            originalnode = unrollednodes[originalnode]
        print(f"Node {node} is unrolled to {originalnode}, merge rows")
        if originalnode not in finalrowsofnodes:
            finalrowsofnodes[originalnode] = set()
        finalrowsofnodes[originalnode].update(rows)
    return finalrowsofnodes

def runmetric(bikoinshapss, bikiinshapss, koshaps, kishaps, outname, propko, propki, binet, inputnames, korows, kirows, koinrows, kiinrows, rowsofnodes, orispeciesnames):
    simko = bikoinshapss[outname] | koshaps[outname]
    simki = bikiinshapss[outname] | kishaps[outname]

    if orispeciesnames:
        print("Consider only original nodes")
        oricalko, oricalki, oripropko, oripropki = dict(), dict(), dict(), dict()
        for spe in orispeciesnames:
            if spe != outname:
                try:
                    oricalko[spe] = simko[spe]
                    oricalki[spe] = simki[spe]
                    oripropko[spe] = propko[spe]
                    oripropki[spe] = propki[spe]
                except:
                    print(f"Cannot find node {spe} in the result lists, set zero")
                    oricalko[spe] = 0.0
                    oricalki[spe] = 0.0
                    oripropko[spe] = 0.0
                    oripropki[spe] = 0.0

    unreachablekeys = simko.keys() - propko.keys() 
    print(f"There are {len(unreachablekeys)} unreachable nodes, they are:")
    print(unreachablekeys)

    errorko = 0.0
    errorki = 0.0 
    num = 0
    showNetwork(binet, bikoinshapss[outname], bikiinshapss[outname], koshaps[outname], kishaps[outname], "binary.html")
    for node, value in simko.items():
        if orispeciesnames:
            if node not in orispeciesnames:
                continue
        num += 1
        if node not in propko:
            propko[node] = 0.0
        if node not in propki:
            propki[node] = 0.0 
        if abs(value - propko[node]) >= 0.005 or abs(simki[node] - propki[node]) >= 0.005:
            print("{:20} - Correct: KO: {:10} | KI : {:10} - Incorrect KO: {:10} | KI: {:10}".format(node, value, simki[node], propko[node], propki[node]))
        errorko += (value - propko[node])*(value - propko[node])
        errorki += (simki[node] - propki[node])*(simki[node] - propki[node])

    del propko[outname]
    del propki[outname]

    wrongrows = 0
    examinednodes = 0

    if not orispeciesnames:
        # ensure both have the same keys
        keys = sorted(simko.keys())

        y_true_ko = np.array([[abs(simko[k]) for k in keys]])
        y_pred_ko = np.array([[abs(propko[k]) for k in keys]])

        y_true_ki = np.array([[abs(simki[k]) for k in keys]])
        y_pred_ki = np.array([[abs(propki[k]) for k in keys]])


        dncg_ko = ndcg_score(y_true_ko, y_pred_ko)
        dncg_ki = ndcg_score(y_true_ki, y_pred_ki)

        if korows and kirows:
            for input in sorted(list(inputnames)):
                if input in rowsofnodes: 
                    gt = koinrows[input].union(kiinrows[input])
                    lack = gt.difference(rowsofnodes[input])
                    extra = rowsofnodes[input].difference(gt)
                    print(input)
                    print("Lack  {:5}: {}".format(len(lack), sorted(list(lack))))
                    print("Extra {:5}: {}".format(len(extra), sorted(list(extra))))
                    print("Real  {:5}: {}".format(len(gt), sorted(list(gt))))
                    print("Prop  {:5}: {}\n".format(len(rowsofnodes[input]), sorted(list(rowsofnodes[input]))))
                    wrongrows += len(lack) + len(extra)
                    examinednodes += 1


            for inter, row in dict(sorted(korows.items(),key=lambda item: item[0])).items():
                if inter not in rowsofnodes:
                    rowsofnodes[inter] = set()
                
                gt = row.union(kirows[inter])
                lack = gt.difference(rowsofnodes[inter])
                extra = rowsofnodes[inter].difference(gt)
                print(inter)
                print("Lack  {:5}: {}".format(len(lack), sorted(list(lack))))
                print("Extra {:5}: {}".format(len(extra), sorted(list(extra))))
                print("Real  {:5}: {}".format(len(gt), sorted(list(gt))))
                print("Prop  {:5}: {}\n".format(len(rowsofnodes[inter]), sorted(list(rowsofnodes[inter]))))
                wrongrows += len(lack) + len(extra)
                examinednodes += 1
                

        print("-----KO RANKING COMPARISON-----")
        print_ranked_keys(propko)
        print("REFERENCE:")
        print_ranked_keys(simko)
        print("=======")
        print("-----KI RANKING COMPARISON-----")
        print_ranked_keys(propki)
        print("REFERENCE:")
        print_ranked_keys(simki)

        print("Number of nodes is {}".format(num))
        print("Number of input nodes is {}".format(len(inputnames)))
        print("Relative_RMSE_KO_is ", (math.sqrt(errorko/num)/(2.0*(len(inputnames)+1))))
        print("Relative_RMSE_KI_is ", (math.sqrt(errorki/num))/(2.0*(len(inputnames)+1)))
        print("NDCG_KO_is ", dncg_ko)
        print("NDCG_KI_is ", dncg_ki)

        korankingres = compare_rankings(simko, propko)
        kirankingres = compare_rankings(simki, propki)
        print_rank_comparison(korankingres, kirankingres)
    else:
        # ensure both have the same keys
        keys = sorted(oricalko.keys())

        y_true_ko = np.array([[abs(oricalko[k]) for k in keys]])
        y_pred_ko = np.array([[abs(oripropko[k]) for k in keys]])

        y_true_ki = np.array([[abs(oricalki[k]) for k in keys]])
        y_pred_ki = np.array([[abs(oripropki[k]) for k in keys]])


        dncg_ko = ndcg_score(y_true_ko, y_pred_ko)
        dncg_ki = ndcg_score(y_true_ki, y_pred_ki)

        if korows and kirows:
            for input in sorted(list(inputnames)):
                if input in rowsofnodes: 
                    gt = koinrows[input].union(kiinrows[input])
                    lack = gt.difference(rowsofnodes[input])
                    extra = rowsofnodes[input].difference(gt)
                    print(input)
                    print("Lack  {:5}: {}".format(len(lack), sorted(list(lack))))
                    print("Extra {:5}: {}".format(len(extra), sorted(list(extra))))
                    print("Real  {:5}: {}".format(len(gt), sorted(list(gt))))
                    print("Prop  {:5}: {}\n".format(len(rowsofnodes[input]), sorted(list(rowsofnodes[input]))))
                    wrongrows += len(lack) + len(extra)
                    examinednodes += 1

            for inter, row in dict(sorted(korows.items(),key=lambda item: item[0])).items():
                if inter not in rowsofnodes:
                    rowsofnodes[inter] = set()

                if inter not in orispeciesnames:
                    continue
                
                gt = row.union(kirows[inter])
                lack = gt.difference(rowsofnodes[inter])
                extra = rowsofnodes[inter].difference(gt)
                print(inter)
                print("Lack  {:5}: {}".format(len(lack), sorted(list(lack))))
                print("Extra {:5}: {}".format(len(extra), sorted(list(extra))))
                print("Real  {:5}: {}".format(len(gt), sorted(list(gt))))
                print("Prop  {:5}: {}\n".format(len(rowsofnodes[inter]), sorted(list(rowsofnodes[inter]))))
                wrongrows += len(lack) + len(extra)
                examinednodes += 1
                

        print("-----KO RANKING COMPARISON-----")
        print_ranked_keys(oripropko)
        print("REFERENCE:")
        print_ranked_keys(oricalko)
        print("=======")
        print("-----KI RANKING COMPARISON-----")
        print_ranked_keys(oripropki)
        print("REFERENCE:")
        print_ranked_keys(oricalki)

        print("Number of nodes is {}".format(num))
        print("Number of input nodes is {}".format(len(inputnames)))
        print("Relative_RMSE_KO_is ", (math.sqrt(errorko/num)/(2.0*(len(inputnames)+1))))
        print("Relative_RMSE_KI_is ", (math.sqrt(errorki/num))/(2.0*(len(inputnames)+1)))
        print("NDCG_KO_is ", dncg_ko)
        print("NDCG_KI_is ", dncg_ki)

        korankingres = compare_rankings(oricalko, oripropko)
        kirankingres = compare_rankings(oricalki, oripropki)
        print_rank_comparison(korankingres, kirankingres)

    avewrongrows = wrongrows/examinednodes 
    averightrows = (pow(2, len(inputnames)) - avewrongrows)/(pow(2, len(inputnames)))
    return averightrows 


