import argparse as ap
from typing import Dict, Any
import hashlib
import json
import math
import copy 
import random 
import numpy as np
from Shapley.exceptions import InforError

from math import isfinite

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



def parseArguments():
    """
    Parse command line arguments.

    Returns:
        argparse.ArgumentParser: The argument parser with defined options. 
    """
    
    parser = ap.ArgumentParser(description='Input the expression file and ')
    parser.add_argument('-e', '--expression', type=str, help="path to expression file", required=True)
    parser.add_argument('-o', '--output', type=str, help='list of interested components to examine')
    parser.add_argument('-k', '--knockout', help='Perform knockout or not', \
                        action='store_true')
    parser.add_argument('-i', '--knockin', help='Perform knockin or not', \
                        action='store_true')
    parser.add_argument('-b', '--binary', help='Work with binary network or not', \
                        action='store_true')
    parser.add_argument('-a', '--acyclic', help='Extract acyclic network with respect to output node', \
                        action='store_true')
    parser.add_argument('-p', '--propagate', help='Run propagation for binary network or not', \
                        action='store_true')
    parser.add_argument('-d', '--debug', help="Print log to debug or not", \
                        action='store_true')

    return parser

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

def toDecimal(state, sortedinput, sortedinter, debug=False):
    """
    Convert the state of a Boolean network to decimal representations for inputs and intermediate nodes.
    Args:
        state (dict): A dictionary with species names as keys and their boolean states as values.
        sortedinput (list): A sorted list of input species names.
        sortedinter (list): A sorted list of intermediate species names.
        debug (bool): If True, print debug information.
    Returns:   
        tuple: A tuple containing two integers - decimal representation of inputs and intermediate nodes.
    """
    try:
        inputbitstr = ''
        for input in sortedinput:
            if state[input]:
                inputbitstr += '1'
            else:
                inputbitstr += '0'
        try:
            interbitstr = ''
            for inter in sortedinter:
                if state[inter]:
                    interbitstr += '1'
                else:
                    interbitstr += '0'
            
            intinputs = int(inputbitstr, 2)
            intinters = int(interbitstr, 2)
            return intinputs, intinters 
        except InforError:
            print("Cannot find intermediate node {} in the state dictionary".format(inter))
    except InforError:
        print("Cannot find input node {} in the state dictionary".format(input))


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
    for item in state1.keys():
        # print(state1[item], state2[item])
        # state1[item] = state1[item] or state2[item] 
        state1[item] = state1[item] and state2[item]
    return state1


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
              
def genTableFromOutput(simoutputs, inputnames, sortedinput, sortedinter, outputnames, debug=False):
    """
    Generate a table from simulation outputs, indexing rows by species states.
    Args:
        simoutputs (list): A list of dictionaries, each representing an output state with species names as keys and their boolean states as values.
        inputnames (set): A set of input species names.
        sortedinput (list): A sorted list of input species names.
        sortedinter (list): A sorted list of intermediate species names.
        outputnames (set): A set of output species names.
        debug (bool): If True, print debug information.
    Returns:
        dict: A dictionary mapping row IDs to dictionaries of species states.
        dict: A dictionary mapping species names to sets of row IDs where they are True.
        dict: A dictionary mapping species names to sets of row IDs where they are False.
    """
    print("----Generate table from output----")
    index = dict() # for each node, save the IDs of row that the node is TRUE
    aindex = dict() # for each node, save the IDs of row that the node is TRUE
    dictresult = dict()
    for inter in sortedinter:
        index[inter] = set()
        aindex[inter] = set()
    for input in inputnames:
        index[input] = set()
        aindex[input] = set()
    for outputname in outputnames:
        index[outputname] = set() 
        aindex[outputname] = set() 
    
    # print("Integer verson of output:")
    for id, line in enumerate(simoutputs):
        # inp, inter = toDecimal(line, sortedinput, sortedinter)
        # print(inp, inter)
        # print(line)
        # each line is a dictionary with keys are name of species and value is true or false
        
        # for test 
        if id in [0, 1]:
            print(f"-------TESTING LINE: {id}-------")
            print(line)

        size = 0
        for input in inputnames:
            if line[input]:
                size += 1
        line['SIZE'] = size
        line['PROP'] = round(math.factorial(size)*math.factorial(len(inputnames) - size)/math.factorial(len(inputnames)),4)
        for inter in sortedinter:
            if line[inter]:
                index[inter].add(id)
            else:
                aindex[inter].add(id)
        for input in inputnames:
            if line[input]:
                index[input].add(id)
            else:
                aindex[input].add(id)
        for outputname in outputnames: 
            if line[outputname]:
                index[outputname].add(id)
            else:
                aindex[outputname].add(id)

        dictresult[id] = line
    return dictresult, index, aindex 

# a wrapper for SHAP, taking a dataset of samples and computes the output of the model for those samples
class BNmodel:
    """ A wrapper class for a Boolean network model to be used with SHAP.
    Attributes:
        map (list): A list mapping input feature indices to species names.
        species (dict): A dictionary with species names as keys and their boolean states as values.
        formulas (list): A list of dictionaries, each containing 'term' and 'formula' keys for the binary formulas.
        outputname (str): The name of the output species to be predicted.
    Methods:
        predict(inputs): Predict the output for a given set of input samples.
        returns: np.ndarray: An array of predicted outputs for the input samples.
    """
    def __init__(self, map, species, formulas, outputname):
        self.map = map # list of input in order 
        self.species = species 
        self.formulas = formulas 
        self.outputname = outputname
    
    def predict(self, inputs): # input is 2 dimension numpy array 
        # print(inputs.shape)
        outputs = np.zeros((inputs.shape[0],1))
        # outputs = []
        # print(outputs.shape)

        # print(len(outputs))
        # print(len(inputs[0]))
        for i in range(inputs.shape[0]):
            # print(inputs[i])
            inputstate = copy.deepcopy(self.species)
            for j in range(inputs.shape[1]):
                if inputs[i][j] == 1:
                    inputstate[self.map[j]] = True
                if inputs[i][j] == 0:
                    inputstate[self.map[j]] = False
            # print(inputstate)
            output = getOutput(self.formulas, inputstate, False, 1000, False)
            # print(output)
            if output[self.outputname]:
                outputs[i][0] = 1
            else:
                outputs[i][0] = 0
        return outputs

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
    kovalues = dict()
    kivalues = dict()
    for node, rowids in rowdict.items():
        kovalues[node] = 0.0
        kivalues[node] = 0.0
        for id in rowids:
            row = simtable[id]
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