import networkx as nx
import copy
from Shapley.exceptions import InforError
from Shapley.utilities import filterrows

from Shapley.visualization import showNetwork 
from Shapley.booleanFormulaHandler import sim1step, sim1bistep, parseFormula
from Shapley.utilities import dict_hash, merge2states, toDecimal
from Shapley.speciesHandler import initSpeciesStates, genInput

# inputstate is a list of all species with the predefine value for input species
# this function is the simulation process 
def getOutput(formulas, inputstate, isbi = False, maxStep = 10000, debug=False, extranodes=None) -> dict: 
    """
    Simulate the Boolean network until it converges or reaches the maximum number of steps. 
    Args:
        formulas (dict): A dictionary with species names as keys and their boolean formulas as values.
        inputstate (dict): A dictionary with species names as keys and their current boolean states as values.
        isbi (bool): If True, use binary simulation step; if False, use standard simulation step.
        maxStep (int): Maximum number of simulation steps.
        debug (bool): If True, print debug information.
        extranodes (list): A list of extra nodes to be considered in the simulation.
    Returns:
        dict: The final state of the network after simulation.
    """
    oldstate = dict()
    numstep = 0
    # print("New inputstate")
    for _ in range(maxStep):
        # print(dict(sorted(inputstate.items())))
        numstep += 1
        if isbi:
            if extranodes:
                for extranode in extranodes:
                    try:
                        rootname = extranode.split("_to_")[0] 
                        if rootname not in inputstate:
                            print("There is no node named {} as the root node of {}".format(rootname, extranode))
                        inputstate[extranode] = inputstate[rootname]
                    except InforError:
                        print("Cannot find the root of the extra node {}".format(extranode)) 
                        return
            inputstate = sim1bistep(formulas, inputstate, debug)
        else:
            # here, if the list of extranodes is not None, 
            # assign extranodes with the value of the root nodes
            if extranodes:
                for extranode in extranodes:
                    try:
                        rootname = extranode.split("_to_")[0] 
                        if rootname not in inputstate:
                            print("There is no node named {} as the root node of {}".format(rootname, extranode))
                        inputstate[extranode] = inputstate[rootname]
                    except InforError:
                        print("Cannot find the root of the extra node {}".format(extranode)) 
                        return

            inputstate = sim1step(formulas, inputstate, debug)
        hash = dict_hash(inputstate) 

        if hash not in oldstate:
            oldstate[hash] = numstep
        else:
            # merge all the state inside the loop 
            returnstate = copy.deepcopy(inputstate)
            for i in range(numstep-oldstate[hash] - 1):
                if isbi:
                    if extranodes:
                        for extranode in extranodes:
                            try:
                                rootname = extranode.split("_to_")[0] 
                                if rootname not in inputstate:
                                    print("There is no node named {} as the root node of {}".format(rootname, extranode))
                                inputstate[extranode] = inputstate[rootname]
                            except InforError:
                                print("Cannot find the root of the extra node {}".format(extranode)) 
                                return
                    inputstate = sim1bistep(formulas, inputstate, debug)
                else:
                    if extranodes:
                        for extranode in extranodes:
                            try:
                                rootname = extranode.split("_to_")[0] 
                                if rootname not in inputstate:
                                    print("There is no node named {} as the root node of {}".format(rootname, extranode))
                                inputstate[extranode] = inputstate[rootname]
                            except InforError:
                                print("Cannot find the root of the extra node {}".format(extranode)) 
                                return
                    inputstate = sim1step(formulas, inputstate, debug)
                returnstate = merge2states(returnstate, inputstate)
            # if toshow:
            #     print("Converge at step {}".format(numstep)) 
            #     print(returnstate)
            return returnstate 
    # print(type(inputstate))
    print("Cannot converge after {} steps".format(maxStep))
    return inputstate 

def getKnockoutOutput(formulas, inputstate, knockoutlist, isbi = False,  \
                      maxStep=1000, debug=False, extranodes = None, isKnockin=False) -> dict:
    """
    Simulate the Boolean network with certain nodes knocked out or knocked in 
    until it converges or reaches the maximum number of steps. 
    Args:
        formulas (dict): A dictionary with species names as keys and their boolean formulas as values.
        inputstate (dict): A dictionary with species names as keys and their current boolean states as values.
        knockoutlist (list): A list of species names to be knocked out (set to False) or knocked in (set to True).
        isbi (bool): If True, use binary simulation step; if False, use standard simulation step.
        maxStep (int): Maximum number of simulation steps.
        debug (bool): If True, print debug information.
        extranodes (list): A list of extra nodes to be considered in the simulation.
        isKnockin (bool): If True, species in knockoutlist are set to True; if False, they are set to False.
    Returns:
        dict: The final state of the network after simulation.
    """
    oldstate = dict()
    numstep = 0
    # print("New inputstate")
    if isKnockin:
        for node in knockoutlist:
            inputstate[node] = True
               
    for _ in range(maxStep):
        # print(dict(sorted(inputstate.items())))
        numstep += 1
        if isbi:
            if extranodes:
                for extranode in extranodes:
                    try:
                        rootname = extranode.split("_to_")[0] 
                        if rootname not in inputstate:
                            print("There is no node named {} as the root node of {}".format(rootname, extranode))
                        inputstate[extranode] = inputstate[rootname]
                    except InforError:
                        print("Cannot find the root of the extra node {}".format(extranode)) 
                        return
            inputstate = sim1bistep(formulas, inputstate, debug, knockoutlist, isKnockin)
        else:
            # here, if the list of extranodes is not None, 
            # assign extranodes with the value of the root nodes
            if extranodes:
                for extranode in extranodes:
                    try:
                        rootname = extranode.split("_to_")[0] 
                        if rootname not in inputstate:
                            print("There is no node named {} as the root node of {}".format(rootname, extranode))
                        inputstate[extranode] = inputstate[rootname]
                    except InforError:
                        print("Cannot find the root of the extra node {}".format(extranode)) 
                        return

            inputstate = sim1step(formulas, inputstate, debug)

            # for extranode in extranodes:
                # print("Extranode {} get value meanwhile rootnode have value {}")

        if knockoutlist:
            for node in knockoutlist:
                if isKnockin:
                    inputstate[node] = True
                    # print("Setting {} is True as it is knocked in".format(node))
                else:
                    inputstate[node] = False
                # print("Setting {} is False as it is knocked out".format(node))

        hash = dict_hash(inputstate) 

        # if 'CI_p' in knockoutlist or 'ci' in knockoutlist:
        #     all_vars = sorted(inputstate.keys())
        #     row = "".join(f"{var: <14}:{int(inputstate.get(var, False)): <2}" for var in all_vars)
        #     print(numstep)
        #     print(row)

        if hash not in oldstate:
            oldstate[hash] = numstep
        else:
            # if debug:
            #     print("Number of iteration {} and first loop point is {}".format(numstep, oldstate[hash]))

            # merge all the state inside the loop 
            returnstate = copy.deepcopy(inputstate)
            for i in range(numstep-oldstate[hash] - 1):
                if isbi:
                    if extranodes:
                        for extranode in extranodes:
                            try:
                                rootname = extranode.split("_to_")[0] 
                                if rootname not in inputstate:
                                    print("There is no node named {} as the root node of {}".format(rootname, extranode))
                                inputstate[extranode] = inputstate[rootname]
                            except InforError:
                                print("Cannot find the root of the extra node {}".format(extranode)) 
                                return
                    inputstate = sim1bistep(formulas, inputstate, debug, knockoutlist, isKnockin)
                else:
                    if extranodes:
                        for extranode in extranodes:
                            try:
                                rootname = extranode.split("_to_")[0] 
                                if rootname not in inputstate:
                                    print("There is no node named {} as the root node of {}".format(rootname, extranode))
                                inputstate[extranode] = inputstate[rootname]
                            except InforError:
                                print("Cannot find the root of the extra node {}".format(extranode)) 
                                return
                    inputstate = sim1step(formulas, inputstate, debug)
                
                if knockoutlist:
                    for node in knockoutlist:
                        if isKnockin:
                            inputstate[node] = True
                        else:
                            inputstate[node] = False
                        # print("Setting {} is False as it is knocked out".format(node))

                returnstate = merge2states(returnstate, inputstate)
            # if toshow:
            #     print("Converge at step {}".format(numstep)) 
            #     print(returnstate)
            return returnstate 
    # print(type(inputstate))
    print("Cannot converge after {} steps".format(maxStep))
    return inputstate 


def simBinaryNetwork(biformulas, inputnames, speciesnames, sortedinput, sortedinter, debug=False, maxstep=1000, extranodes=None):
    """
    Simulate a binary network given its boolean formulas and input states.
    Args:
        biformulas (list): A list of dictionaries, each containing 'term' and 'formula' keys for the binary formulas.
        inputnames (set): A set of input species names.
        speciesnames (set): A set of all species names.
        sortedinput (list): A sorted list of input species names.
        sortedinter (list): A sorted list of intermediate species names.
        debug (bool): If True, print debug information.
        maxstep (int): Maximum number of simulation steps.
        extranodes (list): A list of extra nodes to be considered in the simulation.
    Returns:
        tuple: A tuple containing a list of output states and a dictionary mapping input decimal values to intermediate decimal values.
    """
    print("----Simulate binary network----")
    species = initSpeciesStates(speciesnames, debug) # take list the names of all the species and set them to be None 
    inputstates = genInput(species, inputnames, debug) # each inputstate contains all species, input species are set to a certain value 

    # if debug:
        # print("Integer version of state of Binary network")
    decimalPairs = dict() 
    if extranodes:
        print("Extranodes are not empty, also need to assign value to extranodes")
        for inputstate in inputstates:
            for node in extranodes:
                try:
                    # print(node)
                    rootname = node.split("_to_")[0] 
                    if rootname not in inputstate:
                        print("There is no node named {} as the root node of {}".format(rootname, node))
                    inputstate[node] = inputstate[rootname]
                except:
                    print("Cannot find the root of the extra node {}".format(node)) 
                    return

    outputs = []
    for inputstate in inputstates:
        output = getOutput(biformulas, inputstate, True, 1000, debug, extranodes) 
        
        # this is to test the correctness of the getKnockoutOutput function 
        # output = getKnockoutOutput(biformulas, inputstate, None, True, 1000, False, extranodes)
        
        outputs.append(output)
        # if debug:
            # print(inp, inter)
        inp, inter = toDecimal(output, sortedinput, sortedinter)
        decimalPairs[inp] = inter



    return outputs, decimalPairs

    # run the simulation with the binary network 


def diamondDigger(binet, outname, biformulas):
    """
    Identify and process diamond structures in the binary network starting from the output node.
    Parameters:
        binet: Binary network
        outname: The output node name from which to start the diamond digging
        biformulas: Dictionary of binary boolean formulas for the network
    Returns:
        dict: A dictionary containing information about the diamond structures and their processing results
    """
    print("-------Digging diamonds-------")
    reversenet = binet.reverse(copy=True) # reverse the network to get the in-comming edges
    relevantnodes = nx.descendants(reversenet, outname) # get all the nodes that are relevant to the output node
    relevantnodes.add(outname) # add the output node to the relevant nodes
    processed = set() # save all the processed formulas 
    # relevantgraph = binet.subgraph(relevantnodes) # get the subgraph of the relevant nodes 
    relevantgraph = binet
    showNetwork(relevantgraph, None, None, None, None, "relevantgraph.html")
    curs = [outname] # start from the output node
    carryOns = dict() # save the nodes that are carried on to the next layer
    sinks = dict() # save the sink nodes that are converged
    while curs:
        print("\n\n---Processing layers of {}---".format(curs))
        nextlayer = []
        for cur in curs:
            try:
                form = biformulas[cur]
            except:
                print("\nReach node {} without in-comming edges".format(cur))
                continue
            # get incoming edge to cur node
            inedges = list(relevantgraph.in_edges(cur))
            assert len(inedges) <= 2, print("Support only binary network")

            if cur not in carryOns:
                    carryOns[cur] = dict()
            curset = set(carryOns[cur].keys())

            if len(inedges) == 2:
                print(f"{cur} ==== {form.left.val} {form.val} {form.right.val}")
                stringform = cur + " = " + form.left.val + " " + form.val +  " " + form.right.val
                if stringform in processed:
                    continue
                else:
                    processed.add(stringform)
                
                if form.left.val not in nextlayer:
                    nextlayer.append(form.left.val)
                if form.right.val not in nextlayer:
                    nextlayer.append(form.right.val)

                if form.left.val not in carryOns:
                    carryOns[form.left.val] = copy.deepcopy(carryOns.get(cur, dict()))
                    print(f"Node {form.left.val} carries on {carryOns[form.left.val]} of parent {cur}")
                else:
                    print(f"Node {form.left.val} already carries on {carryOns[form.left.val]}")
                    # if there is common node, add to sinks 
                    if form.left.val not in sinks:
                        sinks[form.left.val] = set()
                    
                    common = set(carryOns[form.left.val].keys()).intersection(curset)
                    foundedsinks = set()
                    for ele in common:
                        if carryOns[form.left.val][ele] == 'W' or carryOns[cur][ele] == 'W':
                            sinks[form.left.val].add(ele)
                            foundedsinks.add(ele)
                        else:
                            if carryOns[form.left.val][ele] == 'L' and carryOns[cur][ele] == 'R':
                                sinks[form.left.val].add(ele)
                                foundedsinks.add(ele)
                            elif carryOns[form.left.val][ele] == 'R' and carryOns[cur][ele] == 'L':
                                sinks[form.left.val].add(ele)
                                foundedsinks.add(ele)

                    print("Found sink nodes for {}: {}".format(form.left.val, foundedsinks))

                    for ele in carryOns[cur].items():
                        if ele[0] not in foundedsinks:
                            carryOns[form.left.val][ele[0]] = ele[1]
                        else:
                            carryOns[form.left.val][ele[0]] = 'W' # if it is sink, then it is converged

                if form.right.val not in carryOns:
                    carryOns[form.right.val] = copy.deepcopy(carryOns.get(cur, dict()))
                    print(f"Node {form.right.val} carries on {carryOns[form.right.val]} of parent {cur}")
                else:
                    print(f"Node {form.right.val} already carries on {carryOns[form.right.val]}")
                    # if there is common node, add to sinks
                    if form.right.val not in sinks:
                        sinks[form.right.val] = set()

                    common = set(carryOns[form.right.val].keys()).intersection(curset)
                    foundedsinks = set()
                    for ele in common:
                        if carryOns[form.right.val][ele] == 'W' or carryOns[cur][ele] == 'W':
                            sinks[form.right.val].add(ele)
                            foundedsinks.add(ele)
                        else:
                            if carryOns[form.right.val][ele] == 'L' and carryOns[cur][ele] == 'R':
                                sinks[form.right.val].add(ele)
                                foundedsinks.add(ele)
                            elif carryOns[form.right.val][ele] == 'R' and carryOns[cur][ele] == 'L':
                                sinks[form.right.val].add(ele)
                                foundedsinks.add(ele)

                    print("Found sink nodes for {}: {}".format(form.right.val, foundedsinks))
                    for ele in carryOns[cur].items():
                        if ele[0] not in foundedsinks:
                            carryOns[form.right.val][ele[0]] = ele[1]
                        else: 
                            carryOns[form.right.val][ele[0]] = 'W'
                   
                
                if cur in carryOns[form.left.val]:
                    if carryOns[form.left.val][cur] == 'R':
                        sinks[form.left.val].add(cur)
                        print("Find {} is sink node for {}".format(cur, form.left.val))
                else:
                    carryOns[form.left.val][cur] = 'L'

                if cur in carryOns[form.right.val]:
                    if carryOns[form.right.val][cur] == 'L':
                        sinks[form.right.val].add(cur)
                        print("Find {} is sink node for {}".format(cur, form.right.val))
                else:
                    carryOns[form.right.val][cur] = 'R'

            else:
                if form.val == 'NOT':
                    stringform = cur + " = NOT " + form.right.val
                    if stringform in processed:
                        continue
                    else:
                        processed.add(stringform)
                    singlemom = form.right.val
                    print(f"{cur} === NOT {singlemom}")
                else:
                    stringform = cur + " = " + form.val
                    if stringform in processed:
                        continue
                    else:
                        processed.add(stringform)
                    singlemom = form.val
                    print(f"{cur} === {singlemom}")
                if singlemom not in nextlayer:
                    nextlayer.append(singlemom)
                
                if singlemom not in carryOns:
                    carryOns[singlemom] = copy.deepcopy(carryOns.get(cur, dict()))
                    print(f"Node {singlemom} carries on {carryOns[singlemom]} of parent {cur}")
                else:
                    print(f"Node {singlemom} already carries on {carryOns[singlemom]}")
                    # if there is common node, add to sinks
                    if singlemom not in sinks:
                        sinks[singlemom] = set()
                    common = set(carryOns[singlemom].keys()).intersection(curset)
                    foundedsinks = set()
                    for ele in common:
                        if carryOns[singlemom][ele] == 'W' or carryOns[cur][ele] == 'W':
                            sinks[singlemom].add(ele)
                            foundedsinks.add(ele)
                        else:
                            if carryOns[singlemom][ele] == 'L' and carryOns[cur][ele] == 'R':
                                sinks[singlemom].add(ele)
                                foundedsinks.add(ele)
                            elif carryOns[singlemom][ele] == 'R' and carryOns[cur][ele] == 'L':
                                sinks[singlemom].add(ele)
                                foundedsinks.add(ele)
                    print("Found sink nodes for {}: {}".format(singlemom, foundedsinks))

                    for ele in carryOns[cur].items():
                        if ele[0] not in foundedsinks:
                            carryOns[singlemom][ele[0]] = ele[1]
                        else:
                            carryOns[singlemom][ele[0]] = 'W'
                    

        for node, carryon in carryOns.items():
            print(f"Node {node} carries on {carryon}")
           
        curs = nextlayer
    for node, diamondends in sinks.items():
        if len(diamondends) == 0:
            continue
        print(f"Node {node} has diamond ends {diamondends}")
    
    return sinks


def convergediamond(node, tobefilteredhalfforeigners, rowsofsinknodes, formulas, index, aindex, extranodes, goodhalfforeignersofnodes, countedrowsofnode, foreignerpartners):
    """
    Converge the rows for a node with branching (diamond structure) by filtering rows based on the states of half-foreigners.
    Args:
        node (str): The name of the node to be processed.
        tobefilteredhalfforeigners (dict): A dictionary mapping nodes to dictionaries of their half-foreigners and their sides.
        rowsofsinknodes (dict): A dictionary mapping sink node names to sets of row IDs associated with those sink nodes.
        formulas (dict): A dictionary with species names as keys and their boolean formulas as values.
        index (dict): A dictionary mapping node names to sets of row IDs where they are True.
        aindex (dict): A dictionary mapping node names to sets of row IDs where they are False.
        extranodes (list): A list of extra nodes to be considered in the filtering.
        goodhalfforeignersofnodes (dict): A dictionary mapping nodes to sets of their good half-foreigners. 
        countedrowsofnode (dict): A dictionary mapping node names to sets of row IDs that have been counted for those nodes.
        foreignerpartners (dict): A dictionary mapping nodes to dictionaries of their half-foreigners and the corresponding foreigner nodes.
    Returns:
        None: The function updates countedrowsofnode in place.
    """
    if node in tobefilteredhalfforeigners:
        print(f"Node {node} has branching and needs to be converge these sink node:")
        # print(tobefilteredofnodes[node])
        for node_, sides in tobefilteredhalfforeigners[node].items():
            print(node_)
    else:
        print(f"Node {node} has no branching, do nothing")
        return
    
    sinks = tobefilteredhalfforeigners[node] 
    for sink, dicts in sinks.items():
        if sink not in rowsofsinknodes:
            continue
        else:
            if len(rowsofsinknodes[sink]) == 0:
                print(f"Sink node {sink} has no rows to be filtered, continue")
                continue
            else:
                print('--------')
                print(f"Sink node {sink} has rows to be filtered:")
                print(list(rowsofsinknodes[sink]))
        
        allrows = copy.deepcopy(rowsofsinknodes[sink])
        leftor = dicts['leftor']
        rightor = dicts['rightor']
        leftand = dicts['leftand']
        rightand = dicts['rightand']
        print(f"Processing sink {sink} with:\n \tleftor {leftor} \n\trightor {rightor} \n\t leftand {leftand} \n\trightand {rightand}")
        sinkop = formulas[sink].val 
        if sinkop == 'OR':
            # filter out all the OR operators 
            # union of left and right rows when filter out and operator 
            for node_ in leftor.union(rightor):
                if node_ not in goodhalfforeignersofnodes[node]:
                    foreigner = foreignerpartners[node][node_]
                    allrows = filterrows('OR', foreigner, node_, index, aindex, extranodes, allrows)

            leftrows = copy.deepcopy(allrows)
            rightrows = copy.deepcopy(allrows) 

            for node_ in leftand.union(rightand): 
                # if a half foreigner is in both side that mean 
                foreigner = foreignerpartners[node][node_]
                allrows = filterrows('AND', foreigner, node_, index, aindex, extranodes, allrows)
            
            for node_ in leftand:
                if node_ not in goodhalfforeignersofnodes[node]:
                    foreigner = foreignerpartners[node][node_]
                    leftrows = filterrows('AND', foreigner, node_, index, aindex, extranodes, leftrows)
            for node_ in rightand:
                if node_ not in goodhalfforeignersofnodes[node]:
                    foreigner = foreignerpartners[node][node_]
                    rightrows = filterrows('AND', foreigner, node_, index, aindex, extranodes, rightrows)

            countedrowsofnode[node].update(leftrows.union(rightrows))
                    
        elif sinkop == 'AND':
            # filter out all the AND operators
            # for OR op, if belong to only one side then no filter 
            for node_ in leftand.union(rightand):
                # if node_ not in goodforeignersofnodes[node]:
                foreigner = foreignerpartners[node][node_]
                allrows = filterrows('AND', foreigner, node_, index, aindex, extranodes, allrows)

            if leftor and rightor:
                print("Both leftor and rightor are not empty, need to filter")
                for node_ in leftor.union(rightor):
                    if node_ not in goodhalfforeignersofnodes[node]:
                        foreigner = foreignerpartners[node][node_]
                        allrows = filterrows('OR', foreigner, node_, index, aindex, extranodes, allrows)

            countedrowsofnode[node].update(allrows) 
        else:
            print("Dot not support operator {}".format(sinkop))

def refineDiamonds(binet, biformulas, sinks):
    """
    Refine diamond structures in a binary network by identifying hybrid children and relevant nodes.
    Args:
        binet (nx.DiGraph): The directed graph representing the binary network.
        biformulas (dict): A dictionary with species names as keys and their boolean formulas as values.
        sinks (dict): A dictionary mapping begin nodes to lists of end nodes that form diamonds with the begin node.
    Returns:
        dict: A dictionary containing the refined diamond structures.
    """
    print("-------Refine diamonds-------")
    reversednet = binet.reverse(copy=True) # reverse the network to get the in-comming edges 
    diamonds = dict() # save the diamonds that are found
    purediamonds = dict() # save the pure diamonds that are found
    dnomaids = dict() # save the diamonds but flipped, end first begin later 
    for begin, ends in sinks.items():
        diamonds[begin] = dict() # save the diamonds for each begin node
        diamonds[begin]['goodchild'] = set() # save the good children of the diamond
        for end in ends:
            if end == 'goodchild':
                continue
            print(f"Processing diamond from {begin} to {end}")
            diamonds[begin][end] = dict() # save the diamond for each end node
            # get end formulas 
            if end not in biformulas:
                print(f"End node {end} is not in the binary formulas, pass")
                continue
            endform = biformulas[end]
            try:
                left = endform.left.val
                right = endform.right.val
                op = endform.val
            except:
                print(f"End node {end} is unary operator, pass")
                continue
            leftdes = nx.descendants(reversednet, left)
            leftdes.add(left) # add the left node to the descendants
            rightdes = nx.descendants(reversednet, right) 
            rightdes.add(right) # add the right node to the descendants
            begindes = nx.descendants(binet, begin) # get all the descendants of the end node 
            begindes.add(begin) # add the begin node to the descendants
            leftrelevant = leftdes.intersection(begindes) # get the relevant nodes of the left side
            rightrelevant = rightdes.intersection(begindes) # get the relevant nodes of the right side
            relevantnodes = leftrelevant.union(rightrelevant) # get the relevant nodes of the diamond
            relevantnodes.add(end) # add the end node to the relevant nodes
            relevantnodes.add(begin) # add the begin node to the relevant nodes
            diamonds[begin][end]['relevant'] = relevantnodes # save the relevant nodes of the diamond

            leftor = set() # save the left OR nodes
            rightor = set() # save the right OR nodes
            leftand = set() # save the left AND nodes
            rightand = set() # save the right AND nodes

            # look for hybidchildren (one parent is insider one parent is outsider) 
            for node in leftrelevant: 
                print(f"Processing node {node} in leftrelevant")
                if node not in biformulas:
                    continue
                form = biformulas[node]
                if form.val == 'OR' or form.val == 'AND':
                    if form.left.val in leftrelevant and form.right.val not in leftrelevant:
                        print(f"Node {node} is a hybrid child of {form.left.val} and {form.right.val}, left is insider, right is outsider")
                        # add the node to the ends
                        diamonds[begin][end][node] = (form.val, form.right.val) 
                        if form.val == 'OR':
                            leftor.add(node)
                        elif form.val == 'AND':
                            leftand.add(node)
                    elif form.right.val in leftrelevant and form.left.val not in leftrelevant:
                        print(f"Node {node} is a hybrid child of {form.left.val} and {form.right.val}, left is outsider, right is insider")
                        diamonds[begin][end][node] = (form.val, form.left.val) 
                        if form.val == 'OR':
                            leftor.add(node)
                        elif form.val == 'AND':
                            leftand.add(node)
                    # elif form.left.val in leftrelevant and form.right.val in leftrelevant:
                        # case both parents are insider 
                else:
                    print(f"Node {node} is unary operator, pass")
                    continue
            
            for node in rightrelevant:
                print(f"Processing node {node} in rightrelevant")
                if node not in biformulas:
                    continue
                form = biformulas[node]
                if form.val == 'OR' or form.val == 'AND':
                    if form.left.val in rightrelevant and form.right.val not in rightrelevant:
                        print(f"Node {node} is a hybrid child of {left} and {right}, left is insider, right is outsider")
                        # add the node to the ends
                        diamonds[begin][end][node] = (form.val, form.right.val) 
                        if form.val == 'OR':
                            rightor.add(node)
                        elif form.val == 'AND':
                            rightand.add(node)
                    elif form.right.val in rightrelevant and form.left.val not in rightrelevant:
                        print(f"Node {node} is a hybrid child of {left} and {right}, left is outsider, right is insider")
                        diamonds[begin][end][node] = (form.val, form.left.val) 
                        if form.val == 'OR':
                            rightor.add(node)
                        elif form.val == 'AND':
                            rightand.add(node)
                    
                else:
                    print(f"Node {node} is unary operator, pass")
                    continue
            if endform.val == 'OR': 
                diamonds[begin][end]['leftor'] = leftor
                diamonds[begin][end]['rightor'] = rightor
                diamonds[begin][end]['leftand'] = leftand
                diamonds[begin][end]['rightand'] = rightand
                diamonds[begin][end]['op'] = 'OR'
               


            elif endform.val == 'AND':  
                diamonds[begin][end]['leftand'] = leftand
                diamonds[begin][end]['rightand'] = rightand
                diamonds[begin][end]['leftor'] = leftor
                diamonds[begin][end]['rightor'] = rightor
                diamonds[begin][end]['op'] = 'AND'
               
                if leftor and rightor:
                    print(f"Node {end} has both leftor and rightor, no good child")
                else:
                    if leftor:
                        diamonds[begin]['goodchild'].update(leftor)
                    elif rightor:
                        diamonds[begin]['goodchild'].update(rightor)
            print("-----------")
            if end not in dnomaids:
                dnomaids[end] = dict() # save the diamonds but flipped, end first begin later
   
            dnomaids[end][begin] = dict() # save the diamonds but flipped, end first begin later
            if endform.val == 'AND': 
                dnomaids[end][begin]['leftand'] = leftand
                dnomaids[end][begin]['rightand'] = rightand
                dnomaids[end][begin]['leftor'] = leftor
                dnomaids[end][begin]['rightor'] = rightor
                dnomaids[end][begin]['op'] = 'AND'
                if leftor and rightor:
                    print(f"Node {end} has both leftor and rightor, no good child")
                else:
                    if leftor:
                        dnomaids[end][begin]['goodchild'] = leftor
                    elif rightor:
                        dnomaids[end][begin]['goodchild'] = rightor
            elif endform.val == 'OR':
                dnomaids[end][begin]['leftor'] = leftor
                dnomaids[end][begin]['rightor'] = rightor
                dnomaids[end][begin]['leftand'] = leftand
                dnomaids[end][begin]['rightand'] = rightand
                dnomaids[end][begin]['op'] = 'OR'
        # print("-------------")
         
    return diamonds, dnomaids


def processDnomaids(nodetoconverge, childtoprocess, foreignerparent, dnomaids, rowstoprocess, index, aindex, extranodes, diamonds, formulas):
    """
    Process the flipped diamonds to refine the rows for a specific node that needs to be converged.
    Args:
        nodetoconverge (str): The name of the node that needs to be converged.
        childtoprocess (str): The name of the child node to be processed.
        foreignerparent (str): The name of the foreigner parent node.
        dnomaids (dict): A dictionary containing the flipped diamond structures.
        rowstoprocess (set): A set of row IDs that need to be processed.
        index (dict): A dictionary mapping node names to sets of row IDs where they are True.
        aindex (dict): A dictionary mapping node names to sets of row IDs where they are False.
        extranodes (list): A list of extra nodes to be considered in the filtering.
        diamonds (dict): A dictionary containing the original diamond structures.
        formulas (dict): A dictionary with species names as keys and their boolean formulas as values.
    Returns:
        set: The refined set of row IDs after processing the flipped diamonds.
    """
    print ("Node {} is in the flipped diamonds".format(childtoprocess)) 

    for begin, dicts in dnomaids[childtoprocess].items(): 
        print("\tProcessing flipped diamond from {} to {}".format(childtoprocess, begin))
        if nodetoconverge in diamonds[begin]['goodchild']:
            print("\tNode {} is a good child of node {}, no need to process".format(childtoprocess, begin))
            continue
        print('\t',dicts) # this one contain leftor, rightor, leftand rightand 
        # get the dnomaid op 
        op = dicts['op']
        typeofinfluence = dict() # save the type of influence of the node to be converged
        for leftor in dicts['leftor']:
            print("\tForeigner of node {} is {}".format(leftor, diamonds[begin][childtoprocess][leftor][1]))
            if nodetoconverge == diamonds[begin][childtoprocess][leftor][1]:
                print("\tNeed to find rows that {} influence {}".format(nodetoconverge, leftor)) 
                if 'leftor' not in typeofinfluence:
                    typeofinfluence['leftor'] = set() 
                typeofinfluence['leftor'].add(leftor)
              
        for rightor in dicts['rightor']:
            print("\tForeigner of node {} is {}".format(rightor, diamonds[begin][childtoprocess][rightor][1]))
            if nodetoconverge == diamonds[begin][childtoprocess][rightor][1]:
                print("\tNeed to find rows that {} influence {}".format(nodetoconverge, rightor))
                if 'rightor' not in typeofinfluence:
                    typeofinfluence['rightor'] = set()
                typeofinfluence['rightor'].add(rightor)
               
        for leftand in dicts['leftand']:
            print("\tForeigner of node {} is {}".format(leftand, diamonds[begin][childtoprocess][leftand][1]))
            if nodetoconverge == diamonds[begin][childtoprocess][leftand][1]:
                print("\tNeed to find rows that {} influence {}".format(nodetoconverge, leftand))
                if 'leftand' not in typeofinfluence:
                    typeofinfluence['leftand'] = set()
                typeofinfluence['leftand'].add(leftand)

        for rightand in dicts['rightand']:
            print("\tForeigner of node {} is {}".format(rightand, diamonds[begin][childtoprocess][rightand][1]))
            if nodetoconverge == diamonds[begin][childtoprocess][rightand][1]:
                print("\tNeed to find rows that {} influence {}".format(nodetoconverge, rightand))
                if 'rightand' not in typeofinfluence:
                    typeofinfluence['rightand'] = set()
                typeofinfluence['rightand'].add(rightand)

        print("\t----Node {} influence node {} in type:----".format(nodetoconverge, childtoprocess))
        for intype, nodes in typeofinfluence.items():
            print("\t\t{}: {}".format(intype, nodes))
        
        if not typeofinfluence:
            print("Node {} also under the influence of node {}, do normal filtering".format(nodetoconverge, begin))
            if op == 'OR':
                print("Take only the rows that {} is False".format(foreignerparent))
                rowstoprocess = rowstoprocess.intersection(aindex[foreignerparent])
            elif op == 'AND':
                print("Take only the rows that {} is True".format(foreignerparent))
                rowstoprocess = rowstoprocess.intersection(index[foreignerparent])
            return rowstoprocess

        for intype, nodes in typeofinfluence.items():
            if intype == 'leftor': 
                for innode in nodes: 
                    # get the partner of the nodetoconverge 
                    partner = formulas[innode].right.val
                    print("Take only the rows that {} is False".format(partner))
                    rowstoprocess = rowstoprocess.intersection(aindex[partner])
            elif intype == 'rightor':
                for innode in nodes:
                    partner = formulas[innode].left.val
                    print("Take only the rows that {} is False".format(partner))
                    rowstoprocess = rowstoprocess.intersection(aindex[partner])

            elif intype == 'leftand':
                for innode in nodes:
                    partner = formulas[innode].right.val
                    print("Take only the rows that {} is True".format(partner))
                    rowstoprocess = rowstoprocess.intersection(index[partner])
                    if op == 'OR':
                        print("And rows that all right or are False")
                        for rightor in dicts['rightor']:
                            print("{} should be False".format(rightor))
                            rowstoprocess = rowstoprocess.intersection(aindex[diamonds[begin][childtoprocess][rightor][1]])
                        rightandrows = set()
                        print("And also rows that at least one rightand False") 
                        for rightand in dicts['rightand']:
                            # print("{} should be False".format(rightand))
                            temrows = rowstoprocess.intersection(aindex[formulas[rightand].left.val])
                            rightandrows.update(temrows)
                        if rightandrows:
                            rowstoprocess = rowstoprocess.intersection(rightandrows)

                    elif op == 'AND':
                        print("With AND operator, only need to take rows that {} is True".format(foreignerparent))
                        rowstoprocess = rowstoprocess.intersection(index[foreignerparent])

            elif intype == 'rightand':
                for innode in nodes:
                    partner = formulas[innode].left.val
                    print("Take only the rows that {} is True".format(partner))
                    rowstoprocess = rowstoprocess.intersection(index[partner])
                    print("After taking rows that {} is True, rows are: \n{}".format(partner, sorted(list(rowstoprocess))))

                    if op == 'OR':
                        print("And rows that all left OR are False")
                        for leftor in dicts['leftor']:
                            print("{} should be False".format(leftor))
                            rowstoprocess = rowstoprocess.intersection(aindex[diamonds[begin][childtoprocess][leftor][1]])
                        leftandrows = set()
                        print("And also only rows that at least one left AND False")
                        for leftand in dicts['leftand']:
                            print("{} should be False".format(leftand))
                            temrows = rowstoprocess.intersection(aindex[formulas[leftand].right.val])
                            leftandrows.update(temrows)
                        if leftandrows:
                            rowstoprocess = rowstoprocess.intersection(leftandrows)
                        
                    elif op == 'AND':
                        print("With AND operator, only need to take rows that {} is True".format(foreignerparent))
                        rowstoprocess = rowstoprocess.intersection(index[foreignerparent])
    print("Rows of node {} after processing flipped diamond {} are: \n{}".format(nodetoconverge, childtoprocess, sorted(list(rowstoprocess))))
    return rowstoprocess

def diamondDigger_old(net, outname, formulas):
    """
    Identify diamond structures in the network starting from the output node.
    A diamond structure is defined as a node that has multiple paths leading to it from a common ancestor node.
    Args:
        net (nx.DiGraph): The directed graph representing the network.
        outname (str): The name of the output node to start the search from.
        formulas (dict): A dictionary with species names as keys and their boolean formulas as values.
    Returns:
        dict: A dictionary mapping nodes to sets of ancestor nodes that form diamond structures.
    """
    print(f"\n\nFinding diamonds for output node {outname}")
    # first round, locate all the diamonds and their endpoints, starting from output node 
    curs = [outname]
    carryonofnodes = dict() 

    while curs:
        print(f'\n\n Processing layer of {curs}')
        nextlayer = []
        for cur in curs: 
            try:
                form = formulas[cur]
            except:
                print(f"Encounter node {cur} without incoming edges")
            inedges = list(net.in_edges(cur))
            assert len(inedges) <= 2, print("Support only binary network")

            if len(inedges) == 2: # care only about binary operator 
                left = form.left.val 
                right = form.right.val 
                # first, carryon  current node cur 
                if left not in carryonofnodes:
                    carryonofnodes[left] = dict() 
                
                if cur not in carryonofnodes[left]:
                    carryonofnodes[left][cur] = []
                carryonofnodes[left][cur].append('L')


                if right not in carryonofnodes:
                    carryonofnodes[right] = dict() 

                if cur not in carryonofnodes[right]:
                    carryonofnodes[right][cur] = []
                carryonofnodes[right][cur].append('R')

                if cur in carryonofnodes:
                    print(f"Now passing carry on of current node {cur} to its parents {left} and {right}")
                    # distribute the carryon of current node to its parents 
                    for node, carryon in carryonofnodes[cur].items():
                        # carryonofnodes[left][node] = copy.deepcopy(carryonofnodes[cur][node]) 
                        # print(node, carryon)
                        if node not in carryonofnodes[left]:
                            carryonofnodes[left][node] = []
                        else:
                            print("Node {} already in carry on of {}".format(node, left))
                        
                        if node not in carryonofnodes[right]:
                            carryonofnodes[right][node] = []
                        else:
                            print("Node {} already in carry on of {}".format(node, right)) 

                        for side in carryon:
                            carryonofnodes[left][node].append(side + 'L')
                            carryonofnodes[right][node].append(side + 'R')
                        
                        # carryonofnodes[right][node].extend(carryon)
                        # carryonofnodes[right][node].append('R')
                        # for side in carryonofnodes[left][node]: 
                        #     # print(side)
                        #     carryonofnodes[left][node].append('L')

                        # carryonofnodes[right][node] = copy.deepcopy(carryonofnodes[cur][node]) 
                        # for side in carryonofnodes[right][node]: 
                        #     # print(side)
                        #     carryonofnodes[right][node].append('R')

                if left not in nextlayer:
                    nextlayer.append(left)
                if right not in nextlayer:
                    nextlayer.append(right)
                
                print("Carryon of {:20}: {}".format(left,carryonofnodes[left]))
                print("Carryon of {:20}: {}".format(right,carryonofnodes[right]))
            elif len(inedges) == 1:
                singlemon = inedges[0][0] 

                #just pass all the carry on of current node cur to its single mom 
                if singlemon not in carryonofnodes:
                    carryonofnodes[singlemon] = dict()

                if cur in carryonofnodes:
                    print(f"Now passing carry on of current node {cur} to its singlemon {singlemon}")
                    for node, carryon in carryonofnodes[cur].items():
                        if node not in carryonofnodes[singlemon]:
                            carryonofnodes[singlemon][node] = []
                        else:
                            print(f"Node {node} already in carry on of node {singlemon}")
                        
                        for side in carryon:
                            carryonofnodes[singlemon][node].append(side)
                        

                if singlemon not in nextlayer:
                    nextlayer.append(singlemon)
                print("Carryon of {:20}: {}".format(singlemon,carryonofnodes[singlemon]))
            else:
                print(f"Node {cur} has no input, pass") 
        curs = nextlayer

    diamonds = dict()
    print ("-----------------Now check DIAMONDS from the resulted carryons--------------------------")
    for node, carryon in carryonofnodes.items():
        print(f"Carryon of {node}")
        if len(net.out_edges(node)) > 1:
            for one, onelist in carryon.items():
                print("{:20}: {}".format(one, onelist))
                L, R = False, False
                for astr in onelist:
                        tostop = False 
                        side = astr[0] 
                        if side == 'L':
                            L = True 
                        if side == 'R':
                            R = True
                if L and R:
                    print("\tNode {} should inherit the carryon rows of node {}".format(node, one))
                    if node not in diamonds:
                        diamonds[node] = set()
                    diamonds[node].add(one)

    return diamonds


def findConstraints(net, output, formulas, debug = False):
    """
    Find constraints for all nodes in the network that influence the given output node.
    Args:
        net (nx.DiGraph): The directed graph representing the network.
        output (str): The name of the output node to start the search from.
        formulas (dict): A dictionary with species names as keys and their boolean formulas as values.
        debug (bool): If True, print debug information.
    Returns: 
        dict: A dictionary mapping nodes to their constraints.
    """
    print("\nFinding constraints for output {}".format(output))
    # convert formulas from list to dictionary
    # formulas pass here is a list of string
    fordict = copy.deepcopy(formulas)
    # for term, formula in formulas.items():
    #     # fordict[formula['left']] = formula['right']
    #     fordict[term] = formula
    # if debug:
    #     print("Formulas are converted to dictionary form")
    #     # print(fordict)

    # traverse the graph from the output to the top to get the filter list of all the nodes 
    filters = dict()
    filters[output] = dict()

    curs = [output] # curs is a list of nodes at the same level (layer)
    while curs:
        if debug:
            print("Current layer includes {}".format(curs))
        nextlayer = []
        for i in range(len(curs)):
            cur = curs.pop(0) 
            # print("Working with {}".format(cur))
            # find constraints for the parent nodes of the current node
            inedges = list(net.in_edges(cur))
            num = len(inedges) # can be only 1 or 2 because of binary tree 
            assert num <= 2, print("Propagate function only supports binary tree at the moment")
            if num == 2: # filter is needed only in case of binary operators 
                edge1 = inedges.pop(0)
                edge2 = inedges.pop(0)
                com1 = edge1[0]
                com2 = edge2[0]

                # check if parent nodes influences multiple nodes from this current layer
                # need to release all the contraints inherit for the dependent node
                if com1 in filters:
                    # print(com1, ":", filters[com1])
                    if cur in filters[com1]:
                        filters[com1].pop(cur)

                if com2 in filters:
                    # print(com2, ":", filters[com2])
                    if cur in filters[com2]:
                        filters[com2].pop(cur)
                
                # if fordict[cur] is a string, convert to Node 
                if isinstance(fordict[cur], str):
                    temdict = dict()
                    temdict['right'] = fordict[cur]
                    temdict['left'] = cur
                    fordict[cur] = parseFormula(temdict, debug) 

                            
                op = fordict[cur].val
                if debug:
                    print(com1, op, com2)
                if op.upper() == "AND":
                    if com1 not in filters:
                        filters[com1] = copy.deepcopy(filters[cur])
                    filters[com1][com2] = True

                    if com2 not in filters:
                        filters[com2] = copy.deepcopy(filters[cur])
                    filters[com2][com1] = True

                elif op.upper() == "OR":
                    if com1 not in filters:
                        filters[com1] = copy.deepcopy(filters[cur])
                    filters[com1][com2] = False
                    
                    if com2 not in filters:
                        filters[com2] = copy.deepcopy(filters[cur])
                    filters[com2][com1] = False
        
                if com1 not in nextlayer:
                    nextlayer.append(com1)
                if com2 not in nextlayer:
                    nextlayer.append(com2)
            elif num == 1: 
                onlyedge = inedges.pop(0)
                onlycom = onlyedge[0]

                if onlycom in filters:
                    print(filters[onlycom])
                    if cur in filters[onlycom]:
                        filters[onlycom].pop(cur)

                if onlycom not in filters:
                    filters[onlycom] = copy.deepcopy(filters[cur])

                if onlycom not in nextlayer:
                    nextlayer.append(onlycom)
                
            else:
                print("Reach node {} without in-comming edges".format(cur))
                continue 
        curs = nextlayer

    if debug:
        for id, item in filters.items():
            print(id,":", item)

    return filters

