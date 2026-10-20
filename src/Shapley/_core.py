import copy 
# import shap 
import numpy as np
import math 
import time
import random
import json

from pathlib import Path

from Shapley.utilities import toDecimal, random_percentage_selection
from Shapley.networkUtilities import getOutputCombined, processMasterDiamonds, simulateDiamond


def simulateNetworkCombined(inputstates, inputnames, internames, formulas, extranodes=None, biformulas=None, debug=False, kolist = None, kilist = None):
    """
    Simulate the network for input states and return the output states, work with both orignal and binarized network.
    Parameters:
        inputstates: List of input states (dictionaries)
        inputnames: list of inputs sorted alphabetical order 
        internames: list of internal nodes sorted alphabetical order
        formulas: Dictionary of boolean formulas for the network
        extranodes: LIST of extra nodes in the network (optional, in case of propagation)
        biformulas: List of dictionaries of binary boolean formulas in the exact order to simulate for the network (optional, in case of propagation)
        debug: Boolean flag for debugging 
        kolist: list of nodes to be knocked out 
        kilist: list of nodes to be knocked in 
    Returns:
        outputs: List of output states (dictionaries), if extranodes and biformulas are provided, the output states will include the values of the extra nodes as well
        blinkings: List of blinking states (dictionaries)
        decimalPairs: Dictionary mapping input decimal representation to output decimal representation
    """
    if extranodes or biformulas:
        print("Simulating original network then inferring values for extranodes and binary formulas")

    outputs = []
    blinkings = []
    decimalPairs = dict()
    for inputstate in inputstates:
        # this output is the output of the original network, has no extra or binary nodes 
        output, blinking, extravalues = getOutputCombined(formulas, inputstate, 1000, debug, extranodes, biformulas, kolist, kilist) 

        inp, inter = toDecimal(output, inputnames, internames)
        decimalPairs[inp] = inter

        # if extravalues:
        output.update(extravalues)

        outputs.append(output)
        blinkings.append(blinking)

    return outputs, blinkings, decimalPairs 


def propagateCyclicVersion1(net, simtable, index, aindex, outname, formulas, simformulas, originalspecies, extranodes, masterDiamonds, inputnames, orderedBiNodes, simpercent = 0, debug=False):
    """
    Propagate values through the network considering diamond structures to determine the rows of nodes.
    Parameters:
        net: The directed graph
        simtable: The simulation table with input-output mappings
        index: Dictionary mapping nodes to rows where they are True 
        aindex: Dictionary mapping nodes to rows where they are False
        outname: The output node name
        formulas: Dictionary of boolean formulas for the network
        extranodes: List of extra nodes added to the network
        diamonds: Dictionary of diamond structures in the network
    Returns:
        countedrowsofnode: Dictionary mapping each node to the set of rows it is associated with
    """
    print("\n\n--------PROPAGATION WITH CYCLIC NETWORK--------")
    if simpercent > 0.0:
        print("Simulate {} percent of node".format(simpercent))
        internames = originalspecies.difference(inputnames)
        if outname in internames:
            internames.remove(outname)
        selected = random_percentage_selection(internames, simpercent)
        print("Selected nodes to simulate:")
        print(selected)
    else:
        selected = set()

    MAX_PROCESS_PER_NODE = 100
    curs = [outname]
    countedrowsofnode = dict() # countedrowsofnode[node] = set of rows that make node count
    rowsofsinknodes = dict() # rows that are about to be filter for a diamond 
    solddiamonds = set() # save the diamonds that are already converged for a node (key) to avoid doing it again
    alreadysimulated = set() # save the nodes that are already simulated to avoid doing it again
    counterofnodes = dict() # count how many times a node is processed to avoid infinite loop 
    changed = True

    while curs and changed:
        changed = False # flag to indicate if any node's counted rows changed in the last iteration
        if debug:
            print("\n\n---Processing layers of {}---".format(curs))
        nextlayer = []
        for cur in curs:
            if '_to_' in cur: # if cur is an extranode that is wired back from an extra node
                root, des = cur.split("_to_")[0], cur.split("_to_")[1] 

                if cur not in countedrowsofnode: # count all rows, should be impossible 
                    # countrows = set(range(len(simtable))) 
                    countrows = set(simtable.keys())
                    countedrowsofnode[cur] = countrows 
                    changed = True 
                    print("WARNING: Extra node {} not counted yet, count all rows to wire back to root node {}".format(cur, root)) 
                else: # count only rows that cur is counted
                    countrows = countedrowsofnode[cur] 
                cur = root 

                # FROM HERE CUR IS WIRED BACK TO ROOT NODE 
                if cur not in countedrowsofnode: # 
                    countedrowsofnode[cur] = countrows 
                    changed = True # cur is touched for the first time, so even countrows is empty, still count as changed
                else: 
                    if countrows.difference(countedrowsofnode[cur]): 
                        changed = True # cur is changed if new rows are added 
                    countedrowsofnode[cur].update(countrows) 
            else:
                if cur not in countedrowsofnode: # count all rows, case of target node 
                    print("At TARGET node {}, count all rows".format(cur))
                    # countrows = set(range(len(simtable))) 
                    countrows = set(simtable.keys())
                    countedrowsofnode[cur] = countrows
                    changed = True # cur is touched for the first time
                else: # count only rows that cur is counted
                    countrows = countedrowsofnode[cur] 
            try:
                form = formulas[cur]
                # if cur == 'PP2A':
                #     print("FORMULA for PP2A is:")
                #     form.display()
            except:
                if debug:
                    print("\nReach node {} without in-comming edges".format(cur))
                continue
            

            # get incoming edge to cur node
            inedges = list(net.in_edges(cur)) # cur here is already wired back if it is an extranode, otherwise normal node
            assert len(inedges) <= 2, print("Support only binary network")
            
            if debug:
                print("\nCounted row for current node {} is \n{}".format(cur, sorted(list(countrows))))
            # if len(inedges) == 2: 
            if form.left and form.right: # to be more robust to support non-binary network, as long as the formula is in the form of (A AND B) or (A OR B)
                if debug:
                    print(f"{cur} ==== {form.left.val} {form.val} {form.right.val}")
                
                # get operator
                op = form.val 
                # check if one of the node is extranode 
                leftselfloop, rightselfloop = False, False
                leftroot, rightroot = None, None
                leftdes, rightdes = None, None
                if form.left.val in extranodes: 
                    if debug:
                        print(f"Left node {form.left.val} is an extra node")
                    # get the root node of the extra node 
                    leftroot = form.left.val.split("_to_")[0]
                    leftdes = form.left.val.split("_to_")[1]
                    if leftroot == leftdes:
                        leftselfloop = True

                if form.right.val in extranodes:
                    if debug:
                        print(f"Right node {form.right.val} is an extra node")
                    # get the root node of the extra node
                    rightroot = form.right.val.split("_to_")[0]
                    rightdes = form.right.val.split("_to_")[1]
                    if rightroot == rightdes:
                        rightselfloop = True

                # now start to process 
                if op == 'OR': 
                    # rows that left counted are rows that right = False 
                    if not rightselfloop:
                        leftrows = aindex[form.right.val]
                    else:
                        leftrows = countedrowsofnode[cur]

                    # rows that right counted are rows that left = False 
                    if not leftselfloop:
                        rightrows = aindex[form.left.val]
                    else:
                        rightrows = countedrowsofnode[cur]
                elif op == 'AND':
                    leftrows = index[form.right.val]
                    # rows that right counted are rows that left = True 
                    rightrows = index[form.left.val]
                else:
                    print(f"Do not support operator {op}")
                    break
                    
                # rowsofsinknodes[cur] = countedrowsofnode[cur] # save the rows of the sink node
                rowsofsinknodes[cur] = countrows
                
                # if form.left.val not in inputnames: 
                if 1: # propagate also to input nodes to estimate error 
                    if form.left.val not in countedrowsofnode:
                        countedrowsofnode[form.left.val] = set() 
                        changed = True # node is touched for the first time
                    
                    if form.left.val in selected and form.left.val not in alreadysimulated: # case to simulate
                        alreadysimulated.add(form.left.val)
                        if debug:
                            print(f"Left node {form.left.val} is a node to SIMULATE, simulate it with all rows")
                        ressimrows = simulateDiamond(simtable, form.left.val, outname, simformulas, set(range(len(simtable))), extranodes, inputnames)
                        countedrowsofnode[form.left.val] = ressimrows
                        if ressimrows.difference(countedrowsofnode[form.left.val]):
                            changed = True

                    if form.left.val not in alreadysimulated:   
                        if form.left.val not in masterDiamonds:
                            # intersect with rows that cur is count 
                            leftrows = leftrows.intersection(countrows)
                            if debug:
                                print(f"After operator {op}, {form.left.val} count only rows: \n{sorted(list(leftrows))}")
                            countedrowsofnode[form.left.val].update(leftrows)
                            if leftrows.difference(countedrowsofnode[form.left.val]):
                                changed = True
                        else: 
                            if debug:
                                print("PROCESS DIAMONDS for {}".format(form.left.val))
                            rowsfromdiamonds = processMasterDiamonds(net, simtable, form.left.val, masterDiamonds, rowsofsinknodes, index, aindex, formulas, solddiamonds, orderedBiNodes, debug)
                            # rowsfromdiamonds = processDiamond(net, simtable, form.left.val, masterDiamonds, diamonds, rowsofsinknodes, index, aindex, extranodes, formulas, solddiamonds)
                            countedrowsofnode[form.left.val].update(rowsfromdiamonds) 
                            if rowsfromdiamonds:
                                changed = True
               
                # if form.right.val not in inputnames: # value for input is already calculated, no need to bother 
                if 1: # propagate also to input nodes to estimate error 
                    if form.right.val not in countedrowsofnode:
                        countedrowsofnode[form.right.val] = set()
                        changed = True
                    
                    if form.right.val in selected and form.right.val not in alreadysimulated:
                        # check probablity to simulate this node 
                        alreadysimulated.add(form.right.val)
                        if debug:
                            print(f"Right node {form.right.val} is a node to SIMULATE, simulate it with all rows")
                        ressimrows = simulateDiamond(simtable, form.right.val, outname, simformulas, set(range(len(simtable))), extranodes, inputnames)
                        countedrowsofnode[form.right.val] = ressimrows
                        if ressimrows:
                            changed = True

                    if form.right.val not in alreadysimulated:
                        if form.right.val not in masterDiamonds:
                            # intersect with rows that cur is count 
                            rightrows = rightrows.intersection(countrows)
                            if debug:
                                print(f"After operator {op}, {form.right.val} count only rows: \n{sorted(list(rightrows))}")
                            countedrowsofnode[form.right.val].update(rightrows)
                            if rightrows.difference(countedrowsofnode[form.right.val]):
                                changed = True
                        else:
                            if debug:
                                print("PROCESS DIAMONDS for {}".format(form.right.val))
                            rowsfromdiamonds = processMasterDiamonds(net, simtable, form.right.val, masterDiamonds, rowsofsinknodes, index, aindex, formulas, solddiamonds, orderedBiNodes, debug)
                            # rowsfromdiamonds = processDiamond(net, simtable, form.right.val, masterDiamonds, diamonds, rowsofsinknodes, index, aindex, extranodes, formulas, solddiamonds)
                            countedrowsofnode[form.right.val].update(rowsfromdiamonds)
                            if rowsfromdiamonds:
                                changed = True

                if form.left.val not in nextlayer and not leftroot: 
                    nextlayer.append(form.left.val)
                if form.right.val not in nextlayer and not rightroot:
                    nextlayer.append(form.right.val)

                # if leftroot and (leftroot not in countedrowsofnode) and (leftroot not in nextlayer):
                if leftroot and (leftroot not in nextlayer):
                    if leftroot in countedrowsofnode:
                        # if leftroot is already counted, increment its counted rows with those from left extra node
                        counterofnodes[leftroot] = counterofnodes.get(leftroot, 0) + 1
                    else:
                        counterofnodes[leftroot] = 1
                    if counterofnodes[leftroot] > MAX_PROCESS_PER_NODE:
                        if debug:
                            print(f"Node {leftroot} has been processed {counterofnodes[leftroot]} times, skip further processing to avoid infinite loop")
                    else:
                        # print(f"Extranodes {form.left.val}, add {form.left.val} to nextlayer")
                        nextlayer.append(form.left.val)
                # if rightroot and (rightroot not in countedrowsofnode) and (rightroot not in nextlayer):
                if rightroot  and (rightroot not in nextlayer):
                    if rightroot in countedrowsofnode:
                        # if rightroot is already counted, increment its counted rows with those from right extra node
                        counterofnodes[rightroot] = counterofnodes.get(rightroot, 0) + 1    
                    else:
                        counterofnodes[rightroot] = 1
                    if counterofnodes[rightroot] > MAX_PROCESS_PER_NODE:
                        if debug:
                            print(f"Node {rightroot} has been processed {counterofnodes[rightroot]} times, skip further processing to avoid infinite loop")
                    else:
                        # print(f"Extranodes {form.right.val}, add  {form.right.val} to nextlayer")
                        nextlayer.append(form.right.val)
            else:
                if form.val == 'NOT':
                    singlemom = form.right.val
                    if debug:
                        print(f"{cur} === NOT {singlemom}")
                else:
                    singlemom = form.val
                    if debug:
                        print(f"{cur} === {singlemom}")

                singleroot, singledes = None, None 
                if singlemom in extranodes:
                    singleroot = singlemom.split("_to_")[0]
                    singledes = singlemom.split("_to_")[1]
                
                
                # if singlemom not in inputnames: 
                if 1: # propagate also to input nodes to estimate error 
                    if singlemom in selected and singlemom not in alreadysimulated:
                        # check probablity to simulate this node 
                        alreadysimulated.add(singlemom)
                        if debug:
                            print(f"Node {singlemom} is a node to simulate, SIMULATE it with all rows")
                        ressimrows = simulateDiamond(simtable, singlemom, outname, simformulas, set(range(len(simtable))), extranodes, inputnames)
                        countedrowsofnode[singlemom] = ressimrows
                        if ressimrows:
                            changed = True

                    if singlemom not in alreadysimulated:
                        if singlemom not in countedrowsofnode:
                            countedrowsofnode[singlemom] = set()
                            changed = True # node is touched for the first time
                        # check if current node cur is in a diamond of singlemom
                        if singlemom not in masterDiamonds:
                            countedrowsofnode[singlemom].update(countrows)
                            if countrows.difference(countedrowsofnode[singlemom]):
                                changed = True
                        else:
                            if debug:
                                print("PROCESS DIAMONDS for {}".format(singlemom))
                            rowsfromdiamonds = processMasterDiamonds(net, simtable, singlemom, masterDiamonds, rowsofsinknodes, index, aindex, formulas, solddiamonds, orderedBiNodes, debug)
                            # rowsfromdiamonds = processDiamond(net, simtable, singlemom, masterDiamonds, diamonds, rowsofsinknodes, index, aindex, extranodes, formulas, solddiamonds)
                            countedrowsofnode[singlemom].update(rowsfromdiamonds)
                            if rowsfromdiamonds:
                                changed = True

                if singlemom not in nextlayer and not singleroot:
                    nextlayer.append(singlemom)
                if singleroot and singleroot not in nextlayer:
                    if singleroot in countedrowsofnode:
                        # if singleroot is already counted, increment its counted rows with those from single extra node
                        counterofnodes[singleroot] = counterofnodes.get(singleroot, 0) + 1
                    else:
                        counterofnodes[singleroot] = 1
                    if counterofnodes[singleroot] > MAX_PROCESS_PER_NODE:
                        if debug:
                            print(f"Node {singleroot} has been processed {counterofnodes[singleroot]} times, skip further processing to avoid infinite loop")
                    else:
                        # print(f"Extranodes {singlemom}, add {singlemom} to nextlayer")
                        nextlayer.append(singlemom)
            if debug:
                print("\t\t\t---------------")
        curs = nextlayer
    if changed:
        print("Warning: Propagation ended but some nodes may not have converged due to cycles.")
        print("Consider increasing MAX_PROCESS_PER_NODE or checking for cycles.")
        print("Counter of processed nodes:")
        for node, count in counterofnodes.items():
            print(f"Node {node}: processed {count} times")
    else:
        print("Propagation converged for all nodes.")
        print("Counter of processed nodes:")
        for node, count in counterofnodes.items():
            print(f"Node {node}: processed {count} times")

    print("SIMULATE {} nodes over {} nodes".format(len(alreadysimulated), len(originalspecies) - len(inputnames)- 1))
    return countedrowsofnode #, potentiallylack, potentiallyextra

