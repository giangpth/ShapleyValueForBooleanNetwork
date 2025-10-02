import networkx as nx
import copy
from Shapley.exceptions import InforError
from Shapley.utilities import filterrows

from Shapley.visualization import showNetwork 
from Shapley.booleanFormulaHandler import sim1step, sim1bistep, checkrow
from Shapley.utilities import dict_hash, merge2states, toDecimal, filterrows
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

def processPossitiveDiamond(beg, target, infordicts, allrows, goodchilds, index, aindex, extranodes):
    """
    Find rows that flow through the inner part of a diamond structure.
    Parameters:
        allrows: The set of all rows to be filtered
        leftors: Set of left OR nodes in the diamond
        rightors: Set of right OR nodes in the diamond
        leftands: Set of left AND nodes in the diamond
        rightands: Set of right AND nodes in the diamond
        externalnodes: Dict of external nodes influencing the diamond through hybrid childs (used as key of dictionary)
        index: Dictionary mapping nodes to rows where they are True
        aindex: Dictionary mapping nodes to rows where they are False
        extranodes: List of extra nodes added to the network
    Returns:
        rowstoreturn: Set of rows that flow through the inner part of the diamond
    """
    rowstofilter = copy.deepcopy(allrows)
    rowstoreturn = set()
    leftors = infordicts['leftor']
    rightors = infordicts['rightor']
    leftands = infordicts['leftand']
    rightands = infordicts['rightand']
    diaop = infordicts['op']
    print(f"Information of diamond from {beg} to {target}: \
                \n\t leftor: {leftors} \n\t rightor: {rightors} \
                \n\t leftand: {leftands} \n\t rightand: {rightands} \
                \n\t op: {diaop} \n\t goodchild: {goodchilds}")
    if diaop == 'OR':
        print("Working with OR diamond")
        # filter out all the OR first
        for node_ in leftors.union(rightors):
            if node_ not in goodchilds:
                try:
                    tofilter = infordicts[node_][1]
                except KeyError:
                    print(f"Cannot find foreigner parent for node {node_}, skip")
                    continue
                print(f"Node {node_} is not a good child, filter out")
                rowstofilter = filterrows('OR', tofilter, node_, index, aindex, extranodes, rowstofilter)
                # print("Rows after filtering out OR {} are: \n{}".format(node_, sorted(list(rowstofilter))))
            else:
                print(f"Node {node_} is a good child, do not filter out")
        if not rowstofilter:
            return set()
        # now filter out and in the left and in the right the merge them since it is needed only one side to pass the signal
        leftrows = copy.deepcopy(rowstofilter)
        rightrows = copy.deepcopy(rowstofilter)
        print("'\n'iltering AND operators on the left")
        for node_ in leftands:
            tofilter = infordicts[node_][1] 
            leftrows = filterrows('AND', tofilter, node_, index, aindex, extranodes, leftrows)
            

        print("\nFiltering AND operators on the right")
        for node_ in rightands:
            tofilter = infordicts[node_][1]
            rightrows = filterrows('AND', tofilter, node_, index, aindex, extranodes, rightrows)

        print("\nMerging left and right rows for AND operators in an OR diamond")
        rowstoreturn.update(leftrows.union(rightrows))
        return rowstoreturn
    
    elif diaop == 'AND':
        print("Working with AND diamond")
        print("Filter out all the AND operator first")
        for node_ in leftands.union(rightands): 
            try:
                tofilter = infordicts[node_][1]
            except KeyError:
                print(f"Cannot find foreigner parent for node {node_}, skip")
                continue
            rowstofilter = filterrows('AND', tofilter, node_, index, aindex, extranodes, rowstofilter)
            # print("Rows after filtering out AND {} are: \n{}".format(node_, sorted(list(rowstofilter))))
        
        # filter OR in only one side, then merge them since both side need to pass the signal
        if not rowstofilter:
            return set()
        leftrows = copy.deepcopy(rowstofilter)
        rightrows = copy.deepcopy(rowstofilter)
        print("\nFiltering OR operators on the left")
        for node_ in leftors:
            if node_ not in goodchilds:
                tofilter = infordicts[node_][1]
                leftrows = filterrows('OR', tofilter, node_, index, aindex, extranodes, leftrows)
                # print("Rows after filtering out OR {} are: \n{}".format(node_, sorted(list(leftrows))))
            else:
                print(f"Node {node_} is a good child, do not filter out")
        print("\nFiltering OR operators on the right")
        for node_ in rightors:
            if node_ not in goodchilds:
                tofilter = infordicts[node_][1]
                rightrows = filterrows('OR', tofilter, node_, index, aindex, extranodes, rightrows)
                # print("Rows after filtering out OR {} are: \n{}".format(node_, sorted(list(rightrows))))
            else:
                print(f"Node {node_} is a good child, do not filter out")
        print("\nMerging left and right rows for OR operators in an AND diamond")
        rowstoreturn.update(leftrows.union(rightrows))
        # if leftors and rightors:
        #     print("Both leftor and rightor are not empty, need to filter")
        #     for node_ in leftors.union(rightors):
        #         if node_ not in goodchilds:
        #             tofilter = infordicts[node_][1]
        #             rowstofilter = filterrows('OR', tofilter, node_, index, aindex, extranodes, rowstofilter)
        #             print("Rows after filtering out {} are: \n{}".format(node_, sorted(list(rowstofilter))))
                # rowstoreturn.update(rowstofilter)
        return rowstoreturn
    else:
        print("Do not support {} operator".format(diaop))
        return set()


def processDiamond(net, table, node, masterDiamonds, diamonds, rowsofsinknodes, index, aindex, extranodes, formulas, solddiamonds):
    """
    Process the biggest diamond structures in the network to infer rows for the node of interest, considering negation.
    Parameters:
        node: The current node being processed
        diamonds: Dictionary of diamond structures in the network
        rowsofsinknodes: Dictionary mapping sink nodes to their associated rows
        index: Dictionary mapping nodes to rows where they are True
        aindex: Dictionary mapping nodes to rows where they are False
        extranodes: List of extra nodes added to the network
        solddiamonds: Set of already processed diamond pairs
        formulas: Dictionary of boolean formulas for the network
    Returns:
        rowsofnodes: Set of rows associated with the current node after processing diamonds 
    """
    # print(f"PROCESS DIAMOND for node {node}") 
    if node in masterDiamonds:
        target = masterDiamonds[node] 
        if target not in rowsofsinknodes:
            print(f"End node {target} has no rows yet to be filtered, continue")
            return set()

        if (node, target) in solddiamonds:
            print(f"Diamond from {node} to {target} has been processed already, continue")
            return set()
        
        if not rowsofsinknodes[target]:
            print(f"End node {target} has no rows to be filtered, continue") 
            return set()
        
        print("Process DIAMOND from {} to {}".format(node, masterDiamonds[node])) 
        # first create set of rows to return first 
        rowstoreturn = set()
        goodchilds = diamonds[node]['goodchild']
        solddiamonds.add((node, target))
        
        dicts = diamonds[node][target] 
        allrows = copy.deepcopy(rowsofsinknodes[target])
        print(f"Rows of node {target} to be filtered are: \n{sorted(list(allrows))}")
        
        leftnots = dicts['leftnot']
        rightnots = dicts['rightnot']

        if not leftnots and not rightnots:
            print(f"Possitive diamond from {node} to {target}")
            setofrows = processPossitiveDiamond(node, target, dicts, allrows, goodchilds, index, aindex, extranodes) 
            rowstoreturn.update(setofrows)
        else: # there is negation 
            try:
                form = formulas[target]
            except KeyError:
                print(f"Cannot find formula for node {target}")
                return set()
            try:
                left = form.left.val
                right = form.right.val
            except:
                print(f"Cannot get left or right parent of node {target}")
            
            if leftnots and rightnots:
                print(f"Both-sided negative diamond from {node} to {target}")
                rowstoreturn = simulateDiamondOneStep(net, table, node, target, formulas, allrows)
            else:
                if leftnots and not rightnots: # only negation on the left 
                    negside = left 
                    nonnegside = right
                    sidename = 'left'
                    antisidename = 'right'
                elif rightnots and not leftnots: # only negation on the right 
                    negside = right 
                    nonnegside = left
                    sidename = 'right'
                    antisidename = 'left'
                print(f"One-sided negative diamond from {node} to {target} with neagation on the {sidename}")
                setofrows = processOnesidedDiamond(sidename, antisidename, node, negside, nonnegside, target, dicts, allrows, goodchilds, index, aindex, extranodes)
                rowstoreturn.update(setofrows)

    return rowstoreturn

def processOnesidedDiamond(side, antiside, beg, neg, nonneg, target, infordicts, allrows, goodchilds, index, aindex, extranodes):
    """
    Process the diamond with negation in only one side.
    Parameters:
        side: the side with negation, either left or right
        node: The current node being processed
        diamonds: Dictionary of diamond structures in the network
        rowsofsinknodes: Dictionary mapping sink nodes to their associated rows
        index: Dictionary mapping nodes to rows where they are True
        aindex: Dictionary mapping nodes to rows where they are False
        extranodes: List of extra nodes added to the network
        solddiamonds: Set of already processed diamond pairs
        formulas: Dictionary of boolean formulas for the network
    Returns:
        rowsofnodes: Set of rows associated with the current node after processing diamonds 
    """
    leftors = infordicts['leftor']
    rightors = infordicts['rightor']
    leftands = infordicts['leftand']
    rightands = infordicts['rightand']
    diaop = infordicts['op']
    print(f"Information of diamond from {beg} to {target}: \
                \n\t leftor: {leftors} \n\t rightor: {rightors} \
                \n\t leftand: {leftands} \n\t rightand: {rightands} \
                \n\t op: {diaop} \n\t goodchild: {goodchilds}")
    
    rowstoreturn = set()
    # first get the information for the infordicts 
    negandkey = side + 'and'
    negands = infordicts[negandkey]

    negorkey = side + 'or'
    negors = infordicts[negorkey]

    nonnegandkey = antiside + 'and'
    nonnegands = infordicts[nonnegandkey]

    nonnegorkey = antiside + 'or'
    nonnegors = infordicts[nonnegorkey] 


    # for both AND and OR opertor, take the rows that beg = neg and negative side has no external influence 
    # first take all rows that beg = neg 
    begnegTrue = allrows.intersection(index[beg]).intersection(index[neg])
    begnegFalse = allrows.intersection(aindex[beg]).intersection(aindex[neg])
    begnegagree = begnegTrue.union(begnegFalse) # can be because of even number of negations that cancel out each other or because of external influence 
    # now filter out the external influence of begnegagree by filter out all OR or AND in the negative side 
    print("Rows that BEGIN and NEG nodes agree are: :\n{}\n".format(sorted(list(begnegagree))))

    # now find rows that beg affect neg, no external OR influence 
    beginfnegnoor = copy.deepcopy(allrows) 
    for negorop in negors:
        if negorop not in goodchilds:
            external = infordicts[negorop][1] 
            beginfnegnoor = beginfnegnoor.intersection(aindex[external])
            print("{} is in OR with external {}, take only rows that {} is False, remaining: \n{}".format(negorop, external, external, sorted(list(beginfnegnoor))))
        else:
            print(f"Encounter goodchild {negorop}, no filter")
    # print("Rows: BEG-NEG-NOOR: \n{}\n".format(sorted(list(beginfnegnoor))))

    beginfnegnoand = copy.deepcopy(allrows)
    for negandop in negands:
        external = infordicts[negandop][1]
        beginfnegnoand = beginfnegnoand.intersection(index[external])
        print("{} is in AND with external {}, take only rows that {} is TRUE, remaining: \n{}".format(negandop, external, external, sorted(list(beginfnegnoand))))
    # print("Rows: BEG-NEG-NOAND: \n{}\n".format(sorted(list(beginfnegnoand))))

    beginfneg = beginfnegnoor.intersection(beginfnegnoand)
    begnoinfneg = allrows.difference(beginfneg)  

    print("BEG INF NEG: \n{}".format(sorted(list(beginfneg))))
    print("BEG NO INF NEG: \n{}".format(sorted(list(begnoinfneg))))
    
    # now find rows that beg affect nonneg, no external influence 
    beginfnonnegnoor = copy.deepcopy(allrows)
    for nonnegorop in nonnegors:
        if nonnegorop not in goodchilds:
            external = infordicts[nonnegorop][1] 
            beginfnonnegnoor = beginfnonnegnoor.intersection(aindex[external])
            print("{} is in OR with external {}, take only rows that {} is False, remaining: \n{}".format(nonnegorop, external, external, sorted(list(beginfnonnegnoor))))
        else:
            print(f"Encounter goodchild {nonnegorop}, no filter")
    # print("Rows: BEG-NONNEG-NOOR: \n{}\n".format(sorted(list(beginfnonnegnoor))))

    beginfnonnegnoand = copy.deepcopy(allrows)
    for nonnegandop in nonnegands:
        external = infordicts[nonnegandop][1] 
        beginfnonnegnoand = beginfnonnegnoand.intersection(index[external])
        print("{} is in AND with external {}, take only rows that {} is TRUE, remaining: \n{}".format(nonnegandop, external, external, sorted(list(beginfnonnegnoand))))
    # print("Rows: BEG-NONNEG-NOAND: \n{}\n".format(sorted(list(beginfnonnegnoand))))

    beginfnonneg = beginfnonnegnoor.intersection(beginfnonnegnoand)
    begnoinfnonneg = allrows.difference(beginfnonneg)

    print("BEG INF NONNEG: \n{}".format(sorted(list(beginfnonneg))))
    print("BEG NO INF NONNEG: \n{}".format(sorted(list(begnoinfnonneg))))

    


    if diaop == 'OR':
        print(f"OR one-sided negative diamnond from {beg} to {target}")

        pureinnerneg = beginfnegnoor.intersection(beginfnegnoand)
        pureinnernonneg = beginfnonnegnoor.intersection(beginfnonnegnoand) 

        innernoor = beginfnegnoor.intersection(beginfnonnegnoor)
        inner = (innernoor.intersection(beginfnegnoand)).union((innernoor.intersection(beginfnonnegnoand)))
        outter = allrows.difference(inner)
        print("INNER : \n{}".format(sorted(list(inner))))
        print("OUTTER: \n{}".format(sorted(list(outter))))

    
        innerbegnegagree = (pureinnerneg.intersection(beginfnonnegnoor)).intersection(begnegagree)
        
        print("First take rows flow INNER and BEG NEG AGREE: \n{}".format(sorted(list(innerbegnegagree))))
        rowstoreturn.update(innerbegnegagree) 

        # negside = outter.intersection(beginfneg.intersection(begnoinfnonneg).intersection(aindex[nonneg]))
        negside = beginfneg.intersection(begnoinfnonneg).intersection(aindex[nonneg])
        print("Take rows BEG inf only NEG and nonneg  is False: \n{}".format(sorted(list(negside)))) 
        rowstoreturn.update(negside)

        nonnegside = outter.intersection(beginfnonneg.intersection(begnoinfneg).intersection(aindex[neg]))
        print("Take rows BEG inf only NONNEG and neg is False: \n{}".format(sorted(list(nonnegside)))) 
        rowstoreturn.update(nonnegside)       

    elif diaop == 'AND':
        print(f"AND one-sided negative diamnond from {beg} to {target}")
        print("Take only rows flow INNER and BEG NEG AGREE: \n{}".format(sorted(list(innerbegnegagree))))
        rowstoreturn.update(innerbegnegagree)

    else:
        print(f"Do not support operator {diaop}")
    return rowstoreturn


def processMasterDiamonds(net, table, node, masterDiamonds, rowsofsinknodes, index, aindex, extranodes, formulas, solddiamond, orderedBiNodes, debug=False): 
    """
    Process the biggest diamond structures in the network to infer rows for the node of interest, considering negation.
    Parameters:
        node: The current node being processed
        diamonds: Dictionary of diamond structures in the network
        rowsofsinknodes: Dictionary mapping sink nodes to their associated rows
        index: Dictionary mapping nodes to rows where they are True
        aindex: Dictionary mapping nodes to rows where they are False
        extranodes: List of extra nodes added to the network
        solddiamonds: Set of already processed diamond pairs
        formulas: Dictionary of boolean formulas for the network
    Returns:
        rowsofnodes: Set of rows associated with the current node after processing diamonds 
    """
    if debug:
        print(f"PROCESS DIAMOND for node {node}") 
    rowsofnode = set() 
    if node in masterDiamonds: 
        if debug:
            print("Process diamond from {} to {}".format(node, masterDiamonds[node])) 
        target = masterDiamonds[node] 
        # goodchilds = diamonds[node]['goodchild']

        if target not in rowsofsinknodes:
            if debug:
                print(f"End node {target} has no rows yet to be filtered, continue")
            return set()

        if (node, target) in solddiamond:
            if debug:
                print(f"Diamond from {node} to {target} has been processed already, continue")
            return set()
        
        solddiamond.add((node, target))

        if not rowsofsinknodes[target]:
            if debug:
                print(f"End node {target} has no rows to be filtered, continue") 
            return set() 
        
        rowstoreturn = simulateDiamondOneStep(net, table, node, target, formulas, rowsofsinknodes[target], orderedBiNodes, debug)
        if debug:
            print("Rows after process diamond {} to {} are: \n{}".format(node, target, sorted(list(rowstoreturn))))
        return rowstoreturn 
        
        dicts = diamonds[node][target] 
        # print(dicts)
        allrows = copy.deepcopy(rowsofsinknodes[target])
        print(f"Rows of node {target} to be filtered are: \n{sorted(list(allrows))}")
        leftors = dicts['leftor']
        rightors = dicts['rightor']
        leftands = dicts['leftand']
        rightands = dicts['rightand']
        leftnots = dicts.get('leftnot', set())
        rightnots = dicts.get('rightnot', set())
        diop = dicts['op']
        print(f"\Information of diamond from {node} to {target} with:\n\tleftor: {leftors} \n\trightor: {rightors} \
                \n\t leftand: {leftands} \n\trightand: {rightands} \
                \n\t leftnot: {leftnots} \n\trightnot: {rightnots} \
                \n\t op: {diop} \n\tgoodchild: {goodchilds}")
        innerrows = processPossitiveDiamond(allrows, leftors, rightors, leftands, rightands, goodchilds, diop, dicts, index, aindex, extranodes)
        outterrows = allrows.difference(innerrows)
        print("Rows flow INNER the diamond are: \n{}".format(sorted(list(innerrows))))
        print("Rows flow OUTTER the diamond are: \n{}".format(sorted(list(outterrows))))
        # if diamond has no not influence 
        if not leftnots and not rightnots:
            print("Positive diamond, take inner rows")
            rowsofnode.update(innerrows)
        else: # in case there is negation
            # first get formula of end of the diamond 
            
            try:
                form = formulas[target]
            except KeyError:
                print(f"Cannot find formula for node {target}")
                return set()
            try:
                left = form.left.val
                right = form.right.val
            except:
                print(f"Cannot get left or right parent of node {target}")

            # get rows that begin node of the diamond is true 
            begintrue = allrows.intersection(index[node])

            # get rows that begin node of the diamond is false 
            beginfalse = allrows.intersection(aindex[node])

            print(f"PROCESS DIAMOND WITH NEGATION {node} to {target}")
            bothside = False 
            if leftnots and not rightnots: # only not in the left 
                negside = left 
                nonnegside = right
                sidename = 'LEFT'
            elif rightnots and not leftnots:
                negside = right 
                nonnegside = left
                sidename = 'RIGHT'
            else: # negation on both side 
                bothside = True 
                # just name side to process late 
                negside = left
                nonnegside = right
            
            # now get rows
            # get the rows that negside is True, False, respectively
            negtrue = allrows.intersection(index[negside])
            negfalse = allrows.intersection(aindex[negside]) 

            # get the rows that nonnegside is True, False, respectively 
            nonnegtrue = allrows.intersection(index[nonnegside])
            nonnegfalse = allrows.intersection(aindex[nonnegside])

            # find rows that begin == negside == true
            begnegtrue = negtrue.intersection(begintrue)
                # find rows that begin == negside == false
            begnegfalse = negfalse.intersection(beginfalse)

            # find rows that begin = false, negside = true
            begtruenegfalse = begintrue.intersection(negfalse)
            # find rows that begin = true, negside = false 
            begfalsenegtrue = beginfalse.intersection(negtrue)

            # find rows that begin == nonnegside = true
            begnonnegtrue = begintrue.intersection(nonnegtrue)
                # find rows that begin == nonnegside = false
            begnonnegfalse = beginfalse.intersection(nonnegfalse)

            # find rows that begin = true, nonnegside = false
            begtruenonnegfalse = begintrue.intersection(nonnegfalse)
            # find rows that begin = false, nonnegside = true
            begfalsenonnegtrue = beginfalse.intersection(nonnegtrue)
            
            print("Rows with {}=TRUE, {}=TRUE: \n{}".format(node, negside, sorted(list(begnegtrue))))
            print("Rows with {}=FALSE, {}=FALSE: \n{}".format(node, negside, sorted(list(begnegfalse))))
            print("Rows with {}=TRUE, {}=FALSE: \n{}".format(node, negside, sorted(list(begtruenegfalse))))
            print("Rows with {}=FALSE, {}=TRUE: \n{}".format(node, negside, sorted(list(begfalsenegtrue))))

            print("Rows with {}=TRUE, {}=TRUE: \n{}".format(node, nonnegside, sorted(list(begnonnegtrue))))
            print("Rows with {}=FALSE, {}=FALSE: \n{}".format(node, nonnegside, sorted(list(begnonnegfalse))))
            print("Rows with {}=TRUE, {}=FALSE: \n{}".format(node, nonnegside, sorted(list(begtruenonnegfalse))))
            print("Rows with {}=FALSE, {}=TRUE: \n{}".format(node, nonnegside, sorted(list(begfalsenonnegtrue))))
            
            if not bothside:
                print(f"{diop} diamond has negation on only {sidename} side of node {left}")
                # take inner rows that begin = negside 
                begnegagree = begnegtrue.union(begnegfalse)
                innerbegnegagree = innerrows.intersection(begnegagree)
                outterbegnegagree = outterrows.intersection(begnegagree)

                begnegdisagree = begtruenegfalse.union(begfalsenegtrue) 
                innerbegnegdisagree = innerrows.intersection(begnegdisagree) 
                outterbegnegdisagree = outterrows.intersection(begnegdisagree) 

                pureinnerbegnegagree = copy.deepcopy(innerbegnegagree)
                pureinnerbegnegdisagree = copy.deepcopy(innerbegnegdisagree)
                
                begnonnegagree = begnonnegtrue.union(begnonnegfalse) 
                innerbegnonnegagree = innerrows.intersection(begnonnegagree)
                outterbegnonnegagree = outterrows.intersection(begnonnegagree)

                begnonnegdisagree = begtruenonnegfalse.union(begfalsenonnegtrue)
                innerbegnonnegdisagree = innerrows.intersection(begnonnegdisagree)
                outterbegnonnegdisagree = outterrows.intersection(begnonnegdisagree)

                pureinnerbegnonnegagree = copy.deepcopy(innerbegnonnegagree)
                pureinnerbegnonnegdisagree = copy.deepcopy(innerbegnonnegdisagree)

                nonnegnosig = copy.deepcopy(allrows) # to save rows that signals CANNOT pass on the side without negation 
                negsig = copy.deepcopy(allrows) # to save rows that signals are passed on the side with negation 
                nonnegsig = copy.deepcopy(allrows) # to save rows that singals are passed on the side without negation 

                # First, get the key name and corresponding dictionary of the diamonds 
                andkeyname = sidename.lower() + 'and' 
                orkeyname = sidename.lower() + 'or'
                negands = dicts[andkeyname] 
                negors = dicts[orkeyname]

                # Then side without negation 
                if sidename.lower() == 'left':
                    nonnegandkeyname = 'rightand'
                    nonnegorkeyname = 'rightor'
                else:
                    nonnegandkeyname = 'leftor'
                    nonnegorkeyname = 'leftor'
                nonnegands = dicts[nonnegandkeyname]
                nonnegors = dicts[nonnegorkeyname]

                if diop == 'OR': # need to examize one-sided AND operators 
                    print(f"OR diamond from {node} to {target}")
                    print("Filter out AND operators on the side with negation")
                    for andop in negands: 
                        
                        # for INNER rows 
                        # filter out all AND operator on the side with negation to make sure that remaining rows are influenced by begin node only 
                        pureinnerbegnegagree = pureinnerbegnegagree.intersection(index[andop])
                        pureinnerbegnegdisagree = pureinnerbegnegdisagree.intersection(index[andop])
                        
                        # for ALL rows 
                        # find rows that signal is completely pass on the side of negation, passing all AND op first 
                        negsig = negsig.intersection(index[andop]) 

                    for orop in negors:
                        # for ALL rows 
                        # find rows that signal is completely pass on the side of negation, passing all OR op
                        negsig = negsig.intersection(aindex[orop]) 
                    
                    # now get rows that signal of begin node cannot pass, which is the complementary of negsig 
                    negnosig = allrows.difference(negsig) 

                    # print("Pure inner begin and neg disagree rows: \n{}".format(sorted(list(pureinnerbegnegdisagree))))

                    for nonandop in nonnegands:
                        # for INNER rows 
                        # filter out all AND op on the side without negation to make sure that remaining rows are influenced by begin only
                        pureinnerbegnonnegagree = pureinnerbegnonnegagree.intersection(index[nonandop])
                        pureinnerbegnonnegdisagree = pureinnerbegnonnegdisagree.intersection(index[nonandop])

                        # now for ALL rows 
                        # find rows that signal is completely pass on the side without negation, first passing all AND op first
                        nonnegsig = nonnegsig.intersection(index[nonandop])

                    # print("Pure inner begin and NON neg agree rows: \n{}".format(sorted(list(pureinnerbegnonnegagree)))) 
                    # print("Pure inner begin and NON neg disagree rows: \n{}".format(sorted(list(pureinnerbegnonnegdisagree))))

                    for nonorop in nonnegors:
                        # for ALL rows 
                        # find rows that signal is completely pass on the side WITHOUT negation, passing all OR op
                        nonnegsig = nonnegsig.intersection(aindex[nonorop])

                    # now get rows that signal of begin node cannot pass, which is the complementary of negsig 
                    negnosig = allrows.difference(negsig) 

                    onlysigneg = outterrows.intersection(negsig.intersection(nonnegnosig))
                    # but rows with only signal in negation side also need to couple with rows that non negation side is False 
                    onlysigneg = onlysigneg.intersection(aindex[nonnegside]) 

                    onlysignonneg = outterrows.intersection(nonnegsig.intersection(negnosig))
                    # same here, rows with only signal from non negation side also need to couple with rows that the negation side is False 
                    onlysignonneg = onlysignonneg.intersection(aindex[negside])

                    print("TAKING following set of rows:")
                    print("Pure inner begin and neg agree rows: \n{}".format(sorted(list(pureinnerbegnegagree))))
                    print("Rows signal only negative side are :\n{}".format(sorted(list(onlysigneg))))
                    print("Rows signal only non negative side are: \n{}".format(sorted(list(onlysignonneg))))
                    
                    rowsofnode.update(pureinnerbegnegagree)
                    rowsofnode.update(onlysigneg)
                    rowsofnode.update(onlysignonneg) 
                
                elif diop == 'AND': 
                    print(f"AND diamond from {node} to {target}") 
                    # for AND operator, allrows are filtered with all AND operators on the way already
                    # if there are or in both size then they are also filtered out, only need to check the case of OR in one side 
                    for negorop in negors: 
                        # filter out all OR op on the side with negation to make sure that remaining rows are influenced by begin node only 
                        pureinnerbegnegagree = pureinnerbegnegagree.intersection(aindex[negorop]) 
                        pureinnerbegnegdisagree = pureinnerbegnegdisagree.intersection(aindex[negorop]) 

                    for nonnegorop in nonnegors:
                        pureinnerbegnonnegagree = pureinnerbegnonnegagree.intersection(aindex[nonnegorop])
                        pureinnerbegnonnegdisagree = pureinnerbegnonnegdisagree.intersection(aindex[nonnegorop]) 

                    # for INNER rows, take only pureinner rows that beg and neg are agree (can caused by double negations or disable negation)
                    print("TAKING following set of rows:")
                    print("Pure inner begin and neg agree rows: \n{}".format(sorted(list(pureinnerbegnegagree))))
                    rowsofnode.update(pureinnerbegnegagree) 


            else:
                # do checkrow here 
                print(f"{diop} diamond from {node} to {target} having negation on BOTH side, SIMULATE a step to decide") 
                for row in allrows:
                    rawrow = copy.deepcopy(table[row])
                    changed = checkrow(net, rawrow, node, target, formulas)
                    if changed:
                        rowsofnode.add(row) 
            
            # print("DIAMOND with NEGATION, simulate a step to decide")
            # for row in allrows:
            #     rawrow = copy.deepcopy(table[row])
            #     changed = checkrow(net, rawrow, node, target, formulas)
            #     if changed:
            #         rowsofnode.add(row)

    print("Rows after process diamond {} to {} are: \n{}".format(node, target, sorted(list(rowsofnode))))
    return rowsofnode

def nodes_on_cycles_digraph(G: nx.DiGraph):
    """
    Identify nodes that are part of cycles in a directed graph.
    Parameters:
        G: A directed graph (networkx DiGraph)
    Returns:
        A set of nodes that are part of cycles in the graph
    """
    on = set()
    for comp in nx.strongly_connected_components(G):
        if len(comp) > 1:
            on |= comp
        else:
            n = next(iter(comp))
            if G.has_edge(n, n):  # self-loop => 1-cycle
                on.add(n)
    return on

def simulateDiamondOneStep(net, table, node, target, formulas, rowsofnode, orderedBiNodes, debug=False):
    """
    Simulate a step to decide the rows for the node of interest through the diamond structure.
    Parameters:
        net: The directed graph representing the network
        table: The truth table of the network
        node: The current node being processed
        target: The target node at the end of the diamond
        formulas: Dictionary of boolean formulas for the network
        rowsofnode: Set of rows associated with the current node before simulation
    Returns:
        rowsofnode: Set of rows associated with the current node after simulation 
    """
    silent = not debug
    if not silent:
        print("SIMULATE rows {} to decide rows for node {} in diamond to {}".format(sorted(list(rowsofnode)), node, target)) 
    newrowsofnode = set()
    node_descendants = nx.descendants(net, node)
    target_succesors = nx.ancestors(net, target)
    affected_nodes = set(node_descendants).intersection(set(target_succesors))
    affected_nodes.add(target)
    if not silent:
        print('Affected nodes between {} and {} are: {}'.format(node, target, affected_nodes))
    affected_nodes_ordered = sorted(affected_nodes, key=orderedBiNodes.__getitem__)

    affectedformulas = []  
    counter = 0
    for n in affected_nodes_ordered:
        if '_XTR_' not in n:
            counter += 1
        try:
            form = formulas[n]
            dictform = dict()
            dictform['term'] = n
            dictform['formula'] = form
            affectedformulas.append(dictform)
        except KeyError:
            print(f"Cannot find formula for node {n}, skip it")
            continue

    for row in rowsofnode:
        # print("Simulate row {}".format(row))
        rawrow = copy.deepcopy(table[row])
        if rawrow[target]:
            oldtarget = True
        else:
            oldtarget = False

        rawrow[node]  = not rawrow[node] # flip the begin node of the diamond

        oldstates = dict()
        numstep = 0
        for _ in range(100):
            numstep += 1
            rawrow = sim1bistep(affectedformulas, rawrow, True, None)
            # print(rawrow)
            hash = dict_hash(rawrow)
            
            if hash not in oldstates:
                oldstates[hash] = numstep
            else: # reach a loop 
                # print("Reach a loop after {} steps".format(numstep))
                returnstate = copy.deepcopy(rawrow)
                for i in range(numstep - oldstates[hash]-1):
                    rawrow = sim1bistep(affectedformulas, rawrow, True, None)
                    returnstate = merge2states(returnstate, rawrow)
                    break
                rawrow = copy.deepcopy(returnstate)
                break
        # print("Converge after {} steps".format(numstep))

        if oldtarget ^ rawrow[target]:
            newrowsofnode.add(row)
    if not silent:       
        print("Rows after simulate a step for diamond {} to {} are: \n{}".format(node, target, sorted(list(newrowsofnode))))
    return newrowsofnode 

def simulateOneNode(table, node, target, formulas, rowsofnode, extranodes, inputnames):
    """
    Simulate to decide the rows for the node of interest.
    Parameters:
        table: The truth table of the network
        node: The current node being processed
        formulas: Dictionary of boolean formulas for the network
        rowsofnode: Set of rows associated with the current node before simulation
        extranodes: List of extra nodes added to the network
    Returns:
        rowsofnode: Set of rows associated with the current node after simulation 
    """
    rowstoreturn = set()
    for row in rowsofnode:
        # print("Simulate row {}".format(row))
        rawrow = copy.deepcopy(table[row])
        oldvalue = rawrow[target]

        for k, value in rawrow.items():
            if k not in inputnames:
                rawrow[k] = False

        knockin = getKnockoutOutput(formulas, rawrow, [node], True, 1000, False, extranodes, True)
        knockout = getKnockoutOutput(formulas, rawrow, [node], True, 1000, False, extranodes, False)

        if (oldvalue ^ knockin[target]) or (oldvalue ^ knockout[target]): 
            rowstoreturn.add(row)

    return rowstoreturn

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
        if debug:
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


def diamondDigger(binet, outname, biformulas, debug=False):
    """
    Identify and process diamond structures in the binary network starting from the output node.
    Parameters:
        binet: Binary network
        outname: The output node name from which to start the diamond digging
        biformulas: Dictionary of binary boolean formulas for the network
    Returns:
        dict: A dictionary containing information about the diamond structures and their processing results
    """
    silent = not debug
    if not silent:
        print("-------Digging diamonds-------")
    # reversenet = binet.reverse(copy=True) # reverse the network to get the in-comming edges
    # relevantnodes = nx.descendants(reversenet, outname) # get all the nodes that are relevant to the output node
    # relevantnodes.add(outname) # add the output node to the relevant nodes
    # processed = set() # save all the processed formulas 
    # relevantgraph = binet.subgraph(relevantnodes) # get the subgraph of the relevant nodes 
    relevantgraph = binet
    # showNetwork(relevantgraph, None, None, None, None, "relevantgraph.html")
    curs = [outname] # start from the output node
    carryOns = dict() # save the nodes that are carried on to the next layer
    sinks = dict() # save the sink nodes that are converged
    while curs:
        if not silent:
            print("\n\n---Processing layers of {}---".format(curs))
        nextlayer = []
        for cur in curs:
            try:
                form = biformulas[cur]
            except:
                if not silent:
                    print("\nReach node {} without in-comming edges".format(cur))
                continue
            # get incoming edge to cur node
            inedges = list(relevantgraph.in_edges(cur))
            assert len(inedges) <= 2, print("Support only binary network")

            if cur not in carryOns:
                    carryOns[cur] = dict()
            curset = set(carryOns[cur].keys())

            if len(inedges) == 2:
                if not silent:
                    print(f"{cur} ==== {form.left.val} {form.val} {form.right.val}")
                # stringform = cur + " = " + form.left.val + " " + form.val +  " " + form.right.val
                # if stringform in processed: # this cause wrong in some case 
                #     continue
                # else:
                #     processed.add(stringform)
                
                if form.left.val not in nextlayer:
                    nextlayer.append(form.left.val)
                if form.right.val not in nextlayer:
                    nextlayer.append(form.right.val)

                if form.left.val not in carryOns:
                    carryOns[form.left.val] = copy.deepcopy(carryOns.get(cur, dict()))
                    if not silent:
                        print(f"Node {form.left.val} carries on {carryOns[form.left.val]} of parent {cur}")
                else:
                    if not silent:
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
                    if not silent:
                        print("Found sink nodes for {}: {}".format(form.left.val, foundedsinks))

                    for ele in carryOns[cur].items():
                        if ele[0] not in foundedsinks:
                            carryOns[form.left.val][ele[0]] = ele[1]
                        else:
                            carryOns[form.left.val][ele[0]] = 'W' # if it is sink, then it is converged

                if form.right.val not in carryOns:
                    carryOns[form.right.val] = copy.deepcopy(carryOns.get(cur, dict()))
                    if not silent:
                        print(f"Node {form.right.val} carries on {carryOns[form.right.val]} of parent {cur}")
                else:
                    if not silent:
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
                    if not silent:
                        print("Found sink nodes for {}: {}".format(form.right.val, foundedsinks))
                    for ele in carryOns[cur].items():
                        if ele[0] not in foundedsinks:
                            carryOns[form.right.val][ele[0]] = ele[1]
                        else: 
                            carryOns[form.right.val][ele[0]] = 'W'
                   
                
                if cur in carryOns[form.left.val]:
                    if carryOns[form.left.val][cur] == 'R':
                        sinks[form.left.val].add(cur)
                        if not silent:
                            print("Find {} is sink node for {}".format(cur, form.left.val))
                else:
                    carryOns[form.left.val][cur] = 'L'

                if cur in carryOns[form.right.val]:
                    if carryOns[form.right.val][cur] == 'L':
                        sinks[form.right.val].add(cur)
                        if not silent:
                            print("Find {} is sink node for {}".format(cur, form.right.val))
                else:
                    carryOns[form.right.val][cur] = 'R'

            else:
                if form.val == 'NOT':
                    # stringform = cur + " = NOT " + form.right.val
                    # if stringform in processed:
                    #     continue
                    # else:
                    #     processed.add(stringform)
                    singlemom = form.right.val
                    if not silent:
                        print(f"{cur} === NOT {singlemom}")
                else:
                    # stringform = cur + " = " + form.val
                    # if stringform in processed:
                    #     continue
                    # else:
                    #     processed.add(stringform)
                    singlemom = form.val
                    if not silent:
                        print(f"{cur} === {singlemom}")
                if singlemom not in nextlayer:
                    nextlayer.append(singlemom)
                
                if singlemom not in carryOns:
                    carryOns[singlemom] = copy.deepcopy(carryOns.get(cur, dict()))
                    if not silent:
                        print(f"Node {singlemom} carries on {carryOns[singlemom]} of parent {cur}")
                else:
                    if not silent:
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
                    if not silent:
                        print("Found sink nodes for {}: {}".format(singlemom, foundedsinks))

                    for ele in carryOns[cur].items():
                        if ele[0] not in foundedsinks:
                            carryOns[singlemom][ele[0]] = ele[1]
                        else:
                            carryOns[singlemom][ele[0]] = 'W'
                    

        # for node, carryon in carryOns.items():
        #     print(f"Node {node} carries on {carryon}")
           
        curs = nextlayer
    if debug:
        for node, diamondends in sinks.items():
            if len(diamondends) == 0:
                continue
            print(f"Node {node} has diamond ends {diamondends}")
    
    return sinks
           
def findMasterDiamond(binet, sinks, debug=False):
    masterDiamonds = dict() 
    for begin, ends in sinks.items():
        biggest = begin
        for end in ends:
            if biggest in nx.ancestors(binet, end):
                biggest = end 
        if debug:
            print("Master diamond of {} is {}".format(begin, biggest))
        masterDiamonds[begin] = biggest 
    return masterDiamonds 


def refineDiamonds(binet, biformulas, sinks, silent=True):
    """
    Processing diamond structures in a binary network by identifying hybrid children and relevant nodes.
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
    # dnomaids = dict() # save the diamonds but flipped, end first begin later 
    for begin, ends in sinks.items():
        diamonds[begin] = dict() # save the diamonds for each begin node
        diamonds[begin]['goodchild'] = set() # save the good children of the diamond, good childrens are the OR only one side of an AND diamond 
        for end in ends:
            if end == 'goodchild':
                continue
            if not silent:
                print(f"Processing diamond from {begin} to {end}")
            diamonds[begin][end] = dict() # save the diamond for each end node
            # get end formulas 
            if end not in biformulas:
                if not silent:
                    print(f"End node {end} is not in the binary formulas, pass")
                continue
            endform = biformulas[end]
            try:
                left = endform.left.val
                right = endform.right.val
                op = endform.val
            except:
                if not silent:
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
            leftnot = set() # save the left NOT nodes
            rightnot = set() # save the right NOT nodes

            # look for hybidchildren (one parent is insider one parent is outsider) 
            for node in leftrelevant: 
                if not silent:
                    print(f"Processing node {node} in leftrelevant")
                if node not in biformulas:
                    continue
                form = biformulas[node]
                if form.val == 'OR' or form.val == 'AND':
                    if form.left.val in leftrelevant and form.right.val not in leftrelevant:
                        if not silent:
                            print(f"Node {node} is a hybrid child of {form.left.val} and {form.right.val}, left is insider, right is outsider")
                        # add the node to the ends
                        diamonds[begin][end][node] = (form.val, form.right.val) 
                        if form.val == 'OR':
                            leftor.add(node)
                        elif form.val == 'AND':
                            leftand.add(node)
                    elif form.right.val in leftrelevant and form.left.val not in leftrelevant:
                        if not silent:
                            print(f"Node {node} is a hybrid child of {form.left.val} and {form.right.val}, left is outsider, right is insider")
                        diamonds[begin][end][node] = (form.val, form.left.val) 
                        if form.val == 'OR':
                            leftor.add(node)
                        elif form.val == 'AND':
                            leftand.add(node)
                    # elif form.left.val in leftrelevant and form.right.val in leftrelevant:
                        # case both parents are insider 
                else:
                    if not silent:
                        print(f"Node {node} is unary operator")
                    if form.val == 'NOT':
                        if node == begin:
                            if not silent:
                                print(f"Node {node} is the begin node, skip")
                            continue
                        if not silent:
                            print(f"Node {node} is a NOT node with parent {form.right.val}")
                        leftnot.add(node)
                    else:
                        if not silent:
                            print(f"Node {node} is unary operator {form.val}, pass")
                        continue
            
            for node in rightrelevant:
                if not silent:
                    print(f"Processing node {node} in rightrelevant")
                if node not in biformulas:
                    continue
                form = biformulas[node]
                if form.val == 'OR' or form.val == 'AND':
                    if form.left.val in rightrelevant and form.right.val not in rightrelevant:
                        if not silent:
                            print(f"Node {node} is a hybrid child of {left} and {right}, left is insider, right is outsider")
                        # add the node to the ends
                        diamonds[begin][end][node] = (form.val, form.right.val) 
                        if form.val == 'OR':
                            rightor.add(node)
                        elif form.val == 'AND':
                            rightand.add(node)
                    elif form.right.val in rightrelevant and form.left.val not in rightrelevant:
                        if not silent:
                            print(f"Node {node} is a hybrid child of {left} and {right}, left is outsider, right is insider")
                        diamonds[begin][end][node] = (form.val, form.left.val) 
                        if form.val == 'OR':
                            rightor.add(node)
                        elif form.val == 'AND':
                            rightand.add(node)
                    
                else:
                    if not silent:
                        print(f"Node {node} is unary operator")
                    if form.val == 'NOT':
                        if node == begin:
                            if not silent:
                                print(f"Node {node} is the begin node, skip")
                            continue
                        if not silent:
                            print(f"Node {node} is a NOT node with parent {form.right.val}")
                        rightnot.add(node)
                    else:
                        if not silent:
                            print(f"Node {node} is unary operator {form.val}, pass")
                        continue
                        
            
            # save the leftor, rightor, leftand, rightand, leftnot, rightnot to the diamonds
            diamonds[begin][end]['leftor'] = leftor
            diamonds[begin][end]['rightor'] = rightor
            diamonds[begin][end]['leftand'] = leftand
            diamonds[begin][end]['rightand'] = rightand
            diamonds[begin][end]['leftnot'] = leftnot
            diamonds[begin][end]['rightnot'] = rightnot
            if endform.val == 'OR': 
                diamonds[begin][end]['op'] = 'OR'
            elif endform.val == 'AND':  
                diamonds[begin][end]['op'] = 'AND'
               
                if leftor and rightor:
                    if not silent:
                        print(f"Node {end} has both leftor and rightor, no good child")
                else:
                    if not leftnot and not rightnot: # goodchild is only for possitive diamond 
                        if leftor:
                            diamonds[begin]['goodchild'].update(leftor)
                        elif rightor:
                            diamonds[begin]['goodchild'].update(rightor)
            # print("-----------")
            # if end not in dnomaids:
            #     dnomaids[end] = dict() # save the diamonds but flipped, end first begin later
   
            # dnomaids[end][begin] = dict() # save the diamonds but flipped, end first begin later
            # dnomaids[end][begin]['leftand'] = leftand
            # dnomaids[end][begin]['rightand'] = rightand
            # dnomaids[end][begin]['leftor'] = leftor
            # dnomaids[end][begin]['rightor'] = rightor
            # dnomaids[end][begin]['leftnot'] = leftnot
            # dnomaids[end][begin]['rightnot'] = rightnot
            # if endform.val == 'AND': 
            #     dnomaids[end][begin]['op'] = 'AND'
            #     if leftor and rightor:
            #         print(f"Node {end} has both leftor and rightor, no good child")
            #     else:
            #         if leftor:
            #             dnomaids[end][begin]['goodchild'] = leftor
            #         elif rightor:
            #             dnomaids[end][begin]['goodchild'] = rightor
            # elif endform.val == 'OR':
            #     dnomaids[end][begin]['op'] = 'OR'
        # print("-------------")
         
    return diamonds #, {'flipped': None} # dnomaids

