import copy 
import networkx as nx
import shap 
import numpy as np

from Shapley.visualization import showNetwork
from Shapley.utilities import parseArguments, readfile,  getOrderedList, rowstovalues, toDecimal, genTableFromOutput, rank_dict_values, filterrows
from Shapley.shapleyCalculation import extractPhe, calKSV4Input, calKSV 
from Shapley.networkUtilities import processDnomaids, refineDiamonds, diamondDigger, getOutput, getKnockoutOutput, simBinaryNetwork
from Shapley.speciesHandler import getSpeciesName, getInputNames, genInput, initSpeciesStates
from Shapley.booleanFormulaHandler import parseFormula, toBinaryFormulas, getFormula
from Shapley.networkHandler import convertBiBooleanFormulas2Network, convertBooleanFormulas2Network, manipulateNetwork

def BooleanShapleyAnalysis():
    """
    Main function to perform Boolean Shapley Analysis.
    """
    parser = parseArguments()
    args= vars(parser.parse_args())
    if len(args) < 2:
        parser.print_help()
    else:
        debug = args['debug']
        isko = args['knockout']
        iski = args['knockin']
        isbi = args['binary']
        isacyclic = args['acyclic']
        exp = args['expression']
        strouts = args['output'].split()
        isprop = args['propagate']
        outputnames = set()
        for strout in strouts:
            outputnames.add(strout)

        print("The expression file is {}".format(exp))
        print("The interested output is: ")
        print(outputnames)
        if isko:
            print("Perform normal procedure for input and knockout procedure for intermediate nodes")
        else:
            print("Perform only normal procedure")
        # now read the expression file and take the string formulas from this 
        lines = readfile(exp, debug)

        speciesnames = getSpeciesName(lines, debug) # get species names 
        print("List of all species:")
        print(speciesnames)

        if not outputnames.issubset(speciesnames):
            print("The list of interested output is not a subset of species")
            return -1
        
        inputnames = getInputNames(lines, speciesnames, debug)
        print("-----Input names-----")
        print(inputnames)

        internames = speciesnames.difference(inputnames).difference(outputnames)
        print("----Intermediate nodes-----")
        print(internames) 


        print("----Getting boolean formulas-----")
        strformulas = getFormula(lines, debug) # list of formulas in strings 
        # formulas = [] # list of formulas in nodes 
        formulas = dict() 
        for strformula in strformulas:
            root = parseFormula(strformula, debug)
            if debug:
                print("Parsing formula for {}".format(strformula['left']))
                root.display()
                print("\n")
            # root.display() 
            # thisfor = {formula['left']: root}
            # formulas.append(thisfor)
            formulas[strformula['left']] = root 
        
        # get the network of the original model, set acyclic = True to delete all the feedback arcs 
        # orinet, anet, aformulas, extranodes = convertBooleanFormulas2Network(formulas, inputnames, \
        #     speciesnames, "network", acyclic=isacyclic, debug=debug)
        orinet, orifas = convertBooleanFormulas2Network(formulas, inputnames, \
            speciesnames, debug=debug)
        
        
        # get the list of the extra node added to the network in case of cycle removement 
        

        # the function returns the simulation output 
        oridecimalpairs, oriinputshapss = workWithOriginalNetwork(orinet, \
            inputnames, speciesnames,outputnames, internames, formulas, isko, iski, debug)

        
        # print ("-----Simulation output is:----------")
        # print(orioutputs)
        # for row in orioutputs:
        #     print(dict(sorted(row.items()))) 

        # for outputname in outputnames:
        #     workWithSHAP(list(inputnames), speciesnames, outputname, formulas, debug)

        if isacyclic:
            print("----------WITH ACYCLIC FLAG-------------")
            # get the acyclic network 
            anet, aformulas, extranodes, nodes_positions = manipulateNetwork(orinet, inputnames, formulas, isacyclic, False, debug)
            if debug:
                for term, formula in aformulas.items():
                    print("{} = ".format(term))
                    formula.display()

            # now do the limiting procedure for each output node 
            adecimalpairs, ashapss = workwithAcyclicNetwork(anet, \
                inputnames, internames, outputnames, speciesnames, aformulas, extranodes, iski, debug)
            
            acount = 0

            for oriinp, oriinter in oridecimalpairs.items():
                if oriinter != adecimalpairs[oriinp]:
                    print(oriinp, oriinter, adecimalpairs[oriinp])
                    acount += 1
            print("\n----TOTAL NUMBER OF DIFFERENCES BETWEEN ORIGNAL NETWORK AND ACYCLIC NETWORK IS {}----\n".format(acount))

            if isbi: # convert the acyclic network to binary network 
                print("-------------Now convert the acyclic network to the binary network---------------")
                sortedinput, sortedinter = getOrderedList(inputnames, internames, True)
                bidecimalpairs = workwithBinaryNetwork(aformulas, inputnames, outputnames, \
                    speciesnames, "abinetwork", sortedinput, sortedinter, isko, iski, debug, extranodes, isprop)
                bicount = 0
                for ainp, ainter, in adecimalpairs.items():
                    if ainter != bidecimalpairs[ainp]:
                        print(ainp, ainter, bidecimalpairs[ainp])
                        bicount += 1
                print("\n----TOTAL NUMBER OF DIFFERENCES BETWEEN ACYCLIC NETWORK AND BINARIZED NETWORK IS {}----\n".format(bicount))


        '''
        # Here is to calculate degree  
        print("------Degree-------")
        btness = nx.degree(orinet)
        print(btness)
        
        # Here is to calculate betweenness 
        print("------Betweenness-------")
        btness = nx.betweenness_centrality(orinet)
        print(btness)

        # Here is to calculate pagerank 
        for outputname in outputnames:
            pagerankseed = dict()
            for innode in inputnames:
                pagerankseed[innode] = oriinputshapss[outputname][innode]
            print("------Pagerank-------")
            # pagerank = nx.pagerank(orinet, personalization=pagerankseed)
            pagerank = nx.pagerank(orinet)
            print(pagerank)
        '''



def workWithSHAP(inputnames, speciesnames, outputname, formulas, debug):
    """
    Work with SHAP to calculate Shapley values for a given output node. 
    Parameters:
        inputnames: List of input node names
        speciesnames: List of all species names
        outputname: The output node name for which Shapley values are calculated
        formulas: Dictionary of boolean formulas for the network
        debug: Boolean flag for debugging
    Returns:
        None
    """
    print("----------Work with SHAP-----------")
    species = initSpeciesStates(speciesnames, debug)

    # inputs = genRandomInputForSHAP(len(inputnames), 32)
    inputs = genAllInputForSHAP(len(inputnames))
    print("--------Input for SHAP------")
    print(inputs)

    model = BNmodel(inputnames, species, formulas, outputname)
    print(model.predict(inputs))

    explainer = shap.Explainer(model.predict, inputs)

    shap_values = explainer(inputs)

    print(inputnames)
    print(shap_values.base_values)
    print(np.round(shap_values.values,3))
    # print(np.round(shap_values.values,4))
    print("----------End working with SHAP-----------")

    # for outputname in outputnames:



def workwithAcyclicNetwork(anet, inputnames, internames, outputnames, speciesnames, aformulas, extranodes, iski, debug):
    """
    Work with the acyclic network to calculate Shapley values.
    Parameters:
        anet: Acyclic network
        inputnames: List of input node names
        internames: List of intermediate node names
        outputnames: List of output node names
        speciesnames: List of all species names
        aformulas: Dictionary of boolean formulas for the acyclic network
        extranodes: List of extra nodes added to the network
        iski: Boolean flag for knock-in analysis
        debug: Boolean flag for debugging
    Returns:
        decimalPairs: Dictionary of input-output pairs in decimal representation
        koinshapss: Dictionary of Shapley values for input nodes
    """
    species = initSpeciesStates(speciesnames, debug)
    inputstates = genInput(species, inputnames, debug)

    sortedinput, sortedinter = getOrderedList(inputnames, internames, True)
    decimalPairs = dict()
    koinshapss = None 

    # initial the extranode with the value of the node it reflect 
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

    cloneinputstates = copy.deepcopy(inputstates)

    outputs = []
    for inputstate in cloneinputstates:
        output = getOutput(aformulas, inputstate, False, 1000, False, extranodes) 
        outputs.append(output)
        inp, inter = toDecimal(output, sortedinput, sortedinter)
        decimalPairs[inp] = inter


    genphe = extractPhe(inputnames, outputnames, outputs)
    koinshapss = calKSV4Input(genphe, inputnames, outputnames)

    kiinshapss = None
    if iski: 
        kiinshapss = calKSV4Input(genphe, inputnames, outputnames, True)

    
    print("----Input shapley value of acyclic network-----")
    for item in koinshapss.items():
        print(item)   
        print('\n')

    for outname in outputnames:
        if kiinshapss:
            showNetwork(anet, koinshapss[outname], kiinshapss[outname], None, None, "acyclic.html")
        else:
            showNetwork(anet, koinshapss[outname], kiinshapss, None, None, "acyclic.html")
    return decimalPairs, koinshapss


def workWithOriginalNetwork(net, inputnames, speciesnames, outputnames, internames, \
                            formulas, isko, iski, debug): 
    """
    Work with the original network to calculate Shapley values.
    Parameters:
        net: Original network
        inputnames: List of input node names
        speciesnames: List of all species names
        outputnames: List of output node names
        internames: List of intermediate node names
        formulas: Dictionary of boolean formulas for the original network
        isko: Boolean flag for knock-out analysis
        iski: Boolean flag for knock-in analysis
        debug: Boolean flag for debugging
    Returns:
        decimalPairs: Dictionary of input-output pairs in decimal representation
        koinputshapss: Dictionary of Shapley values for input nodes
    """
    species = initSpeciesStates(speciesnames, debug)
    inputstates = genInput(species, inputnames, debug)
    cloneinputstates = copy.deepcopy(inputstates)

    sortedinput, sortedinter = getOrderedList(inputnames, internames, True)

    decimalPairs = dict()

    outputs = []
    for inputstate in cloneinputstates:
        output = getOutput(formulas, inputstate, False, 1000, debug)        
        outputs.append(output)
        inp, inter = toDecimal(output, sortedinput, sortedinter) 
        decimalPairs[inp] = inter 

    genphe = extractPhe(inputnames, outputnames, outputs)
    koinputshapss = calKSV4Input(genphe, inputnames, outputnames)
    kiinputshapss = calKSV4Input(genphe, inputnames, outputnames, True)

    
    print("----Input knockout shapley value-----")
    for item in koinputshapss.items():
        print((item))
        print('\n')
    
    print("----Input knockIN shapley value-----")
    for item in kiinputshapss.items():
        print((item))
        print('\n') 
   
    koshaps, kishaps = None, None
    if isko:
        print("-------Now perform the knockout procedure to original network-------")
        # can use the inputstates 
        vs = {}
        for internode in internames:
            if internode not in outputnames:
                print("Knockout {}".format(internode))
                clone2inputstates = copy.deepcopy(inputstates)
                koouputs = []
                for inputstate in clone2inputstates:
                    output = getKnockoutOutput(formulas, inputstate, [internode], False, 1000, False)
                    koouputs.append(output)
                kogenphe = extractPhe(inputnames, outputnames, koouputs)
                vs[internode] = kogenphe 
        # for outname in outputnames:
        koshaps, korows = calKSV(genphe, vs, outputnames, len(inputnames))
        for outname in outputnames:
            print("----KNOCKOUT VALUE for output {}----".format(outname))
            print(dict(sorted(koshaps[outname].items())))
            print("\n")
    
    if iski:
        print("-------Now perform the KNOCKIN procedure to intermediate nodes-------")
        vski = {}
        for internode in internames:
            if internode not in outputnames:
                print("Knockin {}".format(internode))
                clone3inputstates = copy.deepcopy(inputstates)
                kiouputs = []
                for inputstate in clone3inputstates:
                    output = getKnockoutOutput(formulas, inputstate, [internode],\
                                                False, 1000, False, None, True)
                    kiouputs.append(output)
                kigenphe = extractPhe(inputnames, outputnames, kiouputs)
                vski[internode] = kigenphe 
        
        kishaps, kirows = calKSV(genphe, vski, outputnames, len(inputnames))
        for outname in outputnames:
            print("----KNOCKIN VALUE for output {}----".format(outname))
            print(dict(sorted(kishaps.items())))
            print("\n")

    for outname in outputnames:
        if koshaps and kishaps:
            showNetwork(net, koinputshapss[outname], kiinputshapss[outname], \
                    koshaps[outname], kishaps[outname], "original.html")
        else:
            showNetwork(net, koinputshapss[outname], kiinputshapss[outname], \
                    None, None, "original.html")
    return decimalPairs, koinputshapss
           


def propagateUnrollNetwork(unrollednet, simtable, outname, unrolledformulas, index, aindex, extranodes, unrollednodes):
    """
    Propagate values through the unrolled network to determine the rows of nodes.
    Parameters:
        unrollednet: The unrolled directed graph
        simtable: The simulation table with input-output mappings
        outname: The output node name
        unrolledformulas: Dictionary of boolean formulas for the unrolled network
        index: Dictionary mapping nodes to rows where they are True
        aindex: Dictionary mapping nodes to rows where they are False
        extranodes: List of extra nodes added to the network
        unrollednodes: Dictionary mapping unrolled nodes to their original nodes
    Returns:
        rowsofnodes: Dictionary mapping each node to the set of rows it is associated with
    """
    print("\n\n-------PROPAGATING unrolled network-------")
    curs = [outname] # start from the output node
    rowsofnodes = dict() # save the rows of nodes that are counted by this node
    while curs:
        print("\n\n---Processing layers of {}---".format(curs))
        nextlayer = []
        for cur in curs:
            try:
                form = unrolledformulas[cur][0]
                form.display()
            except:
                print("\nReach node {} without in-comming edges".format(cur))
                continue

            # get incoming edge to cur node
            inedges = list(unrollednet.in_edges(cur))
            assert len(inedges) <= 2, print("Support only binary network")
            # while cur in unrollednodes:
            #     cur = unrollednodes[cur]
            #  if cur  is unrolled, then get the original node
            if cur not in rowsofnodes:
                rowsofnodes[cur] = set(range(len(simtable)))
            curset = rowsofnodes[cur]
            print(f"Current rows of node {cur} are: \n{sorted(list(curset))}")

            if len(inedges) == 2:
                print(f"{cur} ==== {form.left.val} {form.val} {form.right.val}")
                if form.left.val not in nextlayer:
                    nextlayer.append(form.left.val)
                if form.right.val not in nextlayer:
                    nextlayer.append(form.right.val)

               
                left = form.left.val
                while left in unrollednodes:
                    left = unrollednodes[left]
                print(f"Left node {form.left.val} is unrolled to {left}")
               
                right = form.right.val
                while right in unrollednodes:
                    right = unrollednodes[right]
                print(f"Right node {form.right.val} is unrolled to {right}") 

                print(f"Rows of node {form.left.val} before inferred from {cur} are: \n{sorted(list(rowsofnodes.get(form.left.val, set()))) }")
                print(f"Rows of node {form.right.val} before inferred from {cur} are: \n{sorted(list(rowsofnodes.get(form.right.val, set()))) }")
                
                if form.val == 'OR':
                    if form.right.val not in rowsofnodes:
                        rowsofnodes[form.right.val] = set()
                    # right takes all the rows that left is False 
                    rowsofnodes[form.right.val].update(curset.intersection(aindex[left]))

                    if form.left.val not in rowsofnodes:
                        rowsofnodes[form.left.val] = set()
                    # left takes all the rows that right is False
                    rowsofnodes[form.left.val].update(curset.intersection(aindex[right]))

                if form.val == 'AND':
                    if form.right.val not in rowsofnodes:
                        rowsofnodes[form.right.val] = set()
                    # right takes all the rows that left is True 
                    rowsofnodes[form.right.val].update(curset.intersection(index[left]))

                    if form.left.val not in rowsofnodes:
                        rowsofnodes[form.left.val] = set()
                    # left takes all the rows that right is True
                    rowsofnodes[form.left.val].update(curset.intersection(index[right]))
                
                print(f"Rows of node {form.left.val} inferred from {cur} are: \n{sorted(list(rowsofnodes[form.left.val]))}")
                print(f"Rows of node {form.right.val} inferred from {cur} are: \n{sorted(list(rowsofnodes[form.right.val]))}")
            else:
                if form.val == 'NOT':
                    singlemom = form.right.val
                    print(f"{cur} === NOT {singlemom}")
                else:
                    singlemom = form.val
                    print(f"{cur} === {singlemom}")

                print(f"Rows of node {singlemom} before inferred from {cur} are: \n{sorted(list(rowsofnodes.get(singlemom, set()))) }")

                if singlemom not in rowsofnodes:
                    rowsofnodes[singlemom] = set()
                rowsofnodes[singlemom].update(curset)
                
                print(f"Rows of node {singlemom} inferred from {cur} are: \n{sorted(list(rowsofnodes[singlemom]))}")

                if singlemom not in nextlayer:
                    nextlayer.append(singlemom)
        curs = nextlayer
    return rowsofnodes


def processDiamonds(node, diamonds, rowsofsinknodes, index, aindex, extranodes, solddiamonds, dnomaids, formulas):
    """
    Process diamond structures in the network to infer rows for nodes.
    Parameters:
        node: The current node being processed
        diamonds: Dictionary of diamond structures in the network
        rowsofsinknodes: Dictionary mapping sink nodes to their associated rows
        index: Dictionary mapping nodes to rows where they are True
        aindex: Dictionary mapping nodes to rows where they are False
        extranodes: List of extra nodes added to the network
        solddiamonds: Set of already processed diamond pairs
        dnomaids: Set of nodes that are part of Dnomaid structures
        formulas: Dictionary of boolean formulas for the network
    Returns:
        rowsofnodes: Set of rows associated with the current node after processing diamonds
    """
    print("-------Processing diamonds for node {}-------".format(node))
    rowsofnodes = set() # save the rows for node infered from diamonds 
    if node in diamonds:
        print(f"Node {node} has diamonds")
        goodchild = diamonds[node]['goodchild'] # get the good children of the diamond
        for end, dicts in diamonds[node].items():
            if end == 'goodchild':
                # print(f"Node {node} has good children {goodchild}, continue")
                continue
            print(f"Processing the diamon {node} to {end}")
            if end not in rowsofsinknodes:
                print(f"Node {end} has no rows yet to be filtered, continue")
                continue
            if (node, end) in solddiamonds:
                print(f"Node {node} to {end} has been processed, continue")
                continue

            solddiamonds.add((node, end))
            if len(rowsofsinknodes[end]) == 0:
                print(f"Node {end} has no rows to be filtered, continue")
                continue
            allrows = copy.deepcopy(rowsofsinknodes[end])
            print(f"Rows of node {end} to be filtered are: \n{sorted(list(allrows))}")
            leftor = dicts['leftor']
            rightor = dicts['rightor']
            leftand = dicts['leftand']
            rightand = dicts['rightand']
            diamondop = dicts['op']
            # goodchild = dicts['goodchild']
            print(f"Processing diamond {node} to {end} with:\n \tleftor {leftor} \n\trightor {rightor} \n\t leftand {leftand} \n\trightand {rightand} \n\top {diamondop} \n\tgoodchild {goodchild}")
            if diamondop == 'OR':
                for node_ in leftor.union(rightor):
                    if node_ not in goodchild:
                        tofilter = diamonds[node][end][node_][1]
                        # if node_ in dnomaids:
                        if False:
                            allrows = processDnomaids(node, node_, tofilter, dnomaids, allrows, index, aindex, extranodes, diamonds, formulas)
                        else:
                            print(f"Node {node_} is not a good child, filter out")
                            allrows = filterrows('OR', tofilter, node_, index, aindex, extranodes, allrows)
                            print("Rows after filtering out {} are: \n{}".format(node_, sorted(list(allrows))))
                    else:
                        print(f"Node {node_} is a good child, do not filter out")
                if not allrows:
                    continue

                leftrows = copy.deepcopy(allrows)
                rightrows = copy.deepcopy(allrows)
                print("Filtering AND operators on the left")
                for node_ in leftand:
                    tofilter = diamonds[node][end][node_][1]
                    # if node_ in dnomaids:
                    #     leftrows = processDnomaids(node, node_, tofilter, dnomaids, leftrows, index, aindex, extranodes, diamonds, formulas)
                    # else:
                    leftrows = filterrows('AND', tofilter, node_, index, aindex, extranodes, leftrows)
                    

                print("Filtering AND operators on the right")
                for node_ in rightand:
                    tofilter = diamonds[node][end][node_][1]
                    # if node_ in dnomaids:
                    #     rightrows = processDnomaids(node, node_, tofilter, dnomaids, rightrows, index, aindex, extranodes, diamonds, formulas)
                    # else:
                    rightrows = filterrows('AND', tofilter, node_, index, aindex, extranodes, rightrows)
                print("Merging left and right rows for AND operators in an OR diamond")
                print("Rows of node {} inferred from diamond {} to {} are: \n{}".format(node, node, end, sorted(list(leftrows.union(rightrows)))))
                rowsofnodes.update(leftrows.union(rightrows))

            elif diamondop == 'AND':

                for node_ in leftand.union(rightand):
                    tofilter = diamonds[node][end][node_][1]
                    # if node_ in dnomaids:
                    #     allrows = processDnomaids(node, node_, tofilter, dnomaids, allrows, index, aindex, extranodes, diamonds, formulas)
                    # else:
                    allrows = filterrows('AND', tofilter, node_, index, aindex, extranodes, allrows)
                if leftor and rightor:
                    print("Both leftor and rightor are not empty, need to filter")
                    for node_ in leftor.union(rightor):
                        if node_ not in goodchild:
                            tofilter = diamonds[node][end][node_][1]
                            # if node_ in dnomaids:
                            if False:
                                allrows = processDnomaids(node, node_, tofilter, dnomaids, allrows, index, aindex, extranodes, diamonds, formulas)
                            else:
                                allrows = filterrows('OR', tofilter, node_, index, aindex, extranodes, allrows)
                                print("Rows after filtering out {} are: \n{}".format(node_, sorted(list(allrows))))
                print("Rows of node {} inferred from diamond {} to {} are: \n{}".format(node, node, end, sorted(list(allrows))))
                rowsofnodes.update(allrows)
    return rowsofnodes

def process1side(net, op, node, side, rowstofilter, index, aindex, extranodes, goodchilds, formulas):
    """
    Process one side of a diamond structure to filter rows based on the path from node to side.
    Parameters:
        net: The directed graph
        op: The operator at the diamond node ('OR' or 'AND')
        node: The current node being processed
        side: The side node of the diamond being processed
        rowstofilter: The set of rows to be filtered
        index: Dictionary mapping nodes to rows where they are True
        aindex: Dictionary mapping nodes to rows where they are False
        extranodes: List of extra nodes added to the network 
        goodchilds: List of good children nodes that should not be filtered
        formulas: Dictionary of boolean formulas for the network
    Returns:
        rowstoreturn: Set of rows after filtering based on the path from node to side   
    """
    if side == node:
        if op == 'OR':
            return set() 
        if op == 'AND':
            return rowstofilter 
    paths = list(nx.all_shortest_paths(net, source=node, target=side))
    rowstoreturn = set ()
    for path in paths: 
        rowsofpath = copy.deepcopy(rowstofilter) 
        print("Process path {}".format(path)) 
        for i in range(len(path) - 2):
            node_ = path[i+1] 
            parent = path[i] 
            if node_ in goodchilds:
                print(f"Node {node_} is a good child, skip")
                continue
            # get partner node 
            op = formulas[node_].val
            if op not in ['OR', 'AND']:
                print(f"Node {node_} is unary operator, skip")
                continue

            partner = formulas[node_].right.val if formulas[node_].left.val == parent else formulas[node_].left.val
            if op == 'OR':
                print(f"Filtering rows of node {node_} with OR operator")
                print("Take only the rows that {} is False".format(partner))
                rowsofpath = filterrows('OR', partner, node_, index, aindex, extranodes, rowsofpath)
            elif op == 'AND':
                print(f"Filtering rows of node {node_} with AND operator")
                print("Take only the rows that {} is True".format(partner))
                rowsofpath = filterrows('AND', partner, node_, index, aindex, extranodes, rowsofpath)
        rowstoreturn.update(rowsofpath)
    return rowstoreturn

def processDiamonds2(net, node, diamonds, rowsofsinknodes, index, aindex, extranodes, solddiamonds, formulas):
    """
    Process diamond structures in the network to infer rows for nodes using path-based filtering.
    Parameters:
        net: The directed graph
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
    rowofnodes = set() # save the rows for node infered from diamonds 
    if node in diamonds: 
        print(f"Node {node} has diamonds")
        goodchild = diamonds[node]['goodchild']
        for end, dicts in diamonds[node].items():
            if end == 'goodchild':
                print(f"Node {node} has good children {goodchild}, continue")
                continue 
            print('------------')
            print(f"Processing the diamond {node} to {end}") 
            if end not in rowsofsinknodes:
                continue
            if (node, end) in solddiamonds:
                print(f"Node {node} to {end} has been processed, continue")
                continue
            solddiamonds.add((node, end))
            if len(rowsofsinknodes[end]) == 0:
                print(f"Node {end} has no rows to be filtered, continue")
                continue
            allrows = copy.deepcopy(rowsofsinknodes[end])
            print(f"Rows of node {end} to be filtered are: \n{sorted(list(allrows))}")
            
            left = formulas[end].left.val
            right = formulas[end].right.val
            op = formulas[end].val
            leftrows = process1side(net, op, node, left, allrows, index, aindex, extranodes, goodchild, formulas)
            rightrows = process1side(net, op, node, right, allrows, index, aindex, extranodes, goodchild, formulas)
            if op == 'OR':
                print("Rows of node {} inferred from diamond {} to {} are: \n{}".format(node, node, end, sorted(list(leftrows.union(rightrows)))))
                rowofnodes.update(leftrows.union(rightrows))
            elif op == 'AND':
                print("Rows of node {} inferred from diamond {} to {} are: \n{}".format(node, node, end, sorted(list(leftrows.intersection(rightrows)))))
                rowofnodes.update(leftrows.intersection(rightrows))
            
    return rowofnodes
        

def propagateWithDiamonds(net, simtable, index, aindex, outname, formulas, extranodes, diamonds, dnomaids):
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
        dnomaids: Set of nodes that are part of Dnomaid structures
    Returns:
        countedrowsofnode: Dictionary mapping each node to the set of rows it is associated with
    """
    print("--------PROPAGATION--------")
    curs = [outname]
    countedrowsofnode = dict() # countedrowsofnode[node] = set of rows that make node count
    rowsofsinknodes = dict() # rows that are about to be filter for a diamond 
    solddiamonds = set() # save the diamonds that are already converged for a node (key) to avoid doing it again
    while curs:
        print("\n\n---Processing layers of {}---".format(curs))
        nextlayer = []
        for cur in curs:
            try:
                form = formulas[cur]
            except:
                print("\nReach node {} without in-comming edges".format(cur))
                continue
            # get incoming edge to cur node
            inedges = list(net.in_edges(cur))
            assert len(inedges) <= 2, print("Support only binary network")
            if cur not in countedrowsofnode: # count all rows
                countrows = set(range(len(simtable)))
                countedrowsofnode[cur] = countrows
            else: # count only rows that cur is counted
                countrows = countedrowsofnode[cur]
            print("\nCounted row for current node {} is \n{}".format(cur, sorted(list(countrows))))
            if len(inedges) == 2: 
                print(f"{cur} ==== {form.left.val} {form.val} {form.right.val}")
                # get operator
                op = form.val 
                # check if one of the node is extranode 
                leftselfloop, rightselfloop = False, False
                rootname = None
                if form.left.val in extranodes:
                    print(f"Left node {form.left.val} is an extra node")
                    # get the root node of the extra node 
                    rootname = form.left.val.split("_to_")[0]
                    desname = form.left.val.split("_to_")[1]
                    if rootname == desname:
                        leftselfloop = True

                if form.right.val in extranodes:
                    print(f"Right node {form.right.val} is an extra node")
                    # get the root node of the extra node
                    rootname = form.right.val.split("_to_")[0]
                    desname = form.right.val.split("_to_")[1]
                    if rootname == desname:
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
                    
                rowsofsinknodes[cur] = countedrowsofnode[cur] # save the rows of the sink node
                
                #rows of current node cur are propagated to parents only in the case cur is not belong to a diamond of the parent 
                # check left parent first 
                if form.left.val not in countedrowsofnode:
                    countedrowsofnode[form.left.val] = set()
                if form.left.val not in diamonds:
                    # intersect with rows that cur is count 
                    leftrows = leftrows.intersection(countrows)
                    print(f"After operator {op}, {form.left.val} count only rows {sorted(list(leftrows))}")
                    countedrowsofnode[form.left.val].update(leftrows)
                else: 
                    print("PROCESS DIAMONDS for {}".format(form.left.val))
                    # rowsfromdiamonds = processDiamonds2(net, form.left.val, diamonds, rowsofsinknodes, index, aindex, extranodes, solddiamonds, formulas)
                    rowsfromdiamonds = processDiamonds(form.left.val, diamonds, rowsofsinknodes, index, aindex, extranodes, solddiamonds, dnomaids, formulas)
                    countedrowsofnode[form.left.val].update(rowsfromdiamonds) 
                
                # check right parent
                if form.right.val not in countedrowsofnode:
                    countedrowsofnode[form.right.val] = set()
                if form.right.val not in diamonds:
                    # intersect with rows that cur is count 
                    rightrows = rightrows.intersection(countrows)
                    print(f"After operator {op}, {form.right.val} count only rows {sorted(list(rightrows))}")
                    countedrowsofnode[form.right.val].update(rightrows)
                else:
                    print("PROCESS DIAMONDS for {}".format(form.right.val))
                    # rowsfromdiamonds = processDiamonds2(net, form.right.val, diamonds, rowsofsinknodes, index, aindex, extranodes, solddiamonds, formulas)
                    rowsfromdiamonds = processDiamonds(form.right.val, diamonds, rowsofsinknodes, index, aindex, extranodes, solddiamonds, dnomaids, formulas)
                    countedrowsofnode[form.right.val].update(rowsfromdiamonds)

                if form.left.val not in nextlayer:
                    nextlayer.append(form.left.val)
                if form.right.val not in nextlayer:
                    nextlayer.append(form.right.val)
            else:
                if form.val == 'NOT':
                    singlemom = form.right.val
                    print(f"{cur} === NOT {singlemom}")
                else:
                    singlemom = form.val
                    print(f"{cur} === {singlemom}")
                
                if singlemom not in countedrowsofnode:
                    countedrowsofnode[singlemom] = set()
                
                # check if current node cur is in a diamond of singlemom
                if singlemom not in diamonds:
                    countedrowsofnode[singlemom].update(countrows)
                else:
                    print("PROCESS DIAMONDS for {}".format(singlemom))
                    # rowsfromdiamonds = processDiamonds2(net, singlemom, diamonds, rowsofsinknodes, index, aindex, extranodes, solddiamonds, formulas)
                    rowsfromdiamonds = processDiamonds(singlemom, diamonds, rowsofsinknodes, index, aindex, extranodes, solddiamonds, dnomaids, formulas)
                    countedrowsofnode[singlemom].update(rowsfromdiamonds)

                if singlemom not in nextlayer:
                    nextlayer.append(singlemom)
        curs = nextlayer
    return countedrowsofnode


# do everything with binary network (convert, get speciesnames, simulate...)          
def workwithBinaryNetwork(formulas, inputnames, outputnames, orispeciesnames, networkname, sortedinput, sortedinter, isko = False, iski = False, debug=False, extranodes=None, isprop=False):
    """
    Work with the binary network to calculate Shapley values and perform knockout/knockin analyses.
    Parameters:
        formulas: List of boolean formulas for the original network
        inputnames: List of input node names
        outputnames: List of output node names      
        orispeciesnames: Set of original species names
        networkname: Name of the network
        sortedinput: List of sorted input node names
        sortedinter: List of sorted intermediate node names
        isko: Boolean flag for knockout analysis
        iski: Boolean flag for knockin analysis
        debug: Boolean flag for debug mode
        extranodes: List of extra nodes added to the network 
        isprop: Boolean flag to indicate if propagation with diamonds should be used
    Returns:
        None
    """
    print("\n\n-------WORKING WITH BINARY NETWORK-------")
    bistrformulas = toBinaryFormulas(formulas, True)
    binet, nodes_positions = convertBiBooleanFormulas2Network(bistrformulas, inputnames, orispeciesnames, "bi" + networkname, False, debug, extranodes) 
    
    for outname in outputnames:
        showNetwork(binet, None, None, None, None, 'binary.html')


    biformulas = [] # this list of dictionary is for the simulation, maintaining correctly the order for the simulation 
    biformulasdict = dict() # this is for the expanding function 
    bispeciesnames = set() 
    biformulasdictandorder = dict() # this is for the unroll function, to keep the order of the binary formulas
    # for term, bistrformula in bistrformulas.items():
    for i, formula in enumerate(bistrformulas):
        term = formula['term']
        bistrformula = formula['formula']
        thisfor = {'left': term, "right": bistrformula}

        thisbiformula = dict()
        thisbiformula['term'] = term
        thisbiformula['formula'] = parseFormula(thisfor, debug) 
        biformulasdict[term] = thisbiformula['formula'] 

        forunroll = copy.deepcopy(thisbiformula['formula']) # make a copy of the formula for unrolling later
        biformulasdictandorder[term] = (forunroll, i) # save the order of the binary formulas, also for the unrolling later

        biformulas.append(thisbiformula)

        bispeciesnames.add(term)
        coms = bistrformula.split()
        for com in coms:
            if com != '(' and com != ')' and com != 'AND' and com != 'OR' and com != 'NOT' and com != '=':
                if '-' in com: # this is for excluding '-' in the formulas which is confusing for the compiler 
                    com = com.replace('-','_')
                bispeciesnames.add(com)
        
    bispecies = initSpeciesStates(bispeciesnames, debug)

    biinternames = bispeciesnames.difference(inputnames).difference(outputnames)
    if debug:
        print(f"----There are {len(biinternames)} intermediate nodes in binarized network:-----")
        print(biinternames)

    # simulate binary network, ONLY TO TEST THE CONSISTENCY WITH THE ORIGINAL NETWORK 
    realbioutputs, bidecimalpairs = simBinaryNetwork(biformulas, inputnames, bispeciesnames, sortedinput, sortedinter, False, 1000, extranodes)

    table, index, aindex = genTableFromOutput(realbioutputs, inputnames, sortedinput, biinternames, outputnames, True) 

    intactbigenphe = extractPhe(inputnames, outputnames, realbioutputs)
    if debug:
        print("--------Binary network maps between genotype and phenotype-------")
        print("Uncomment to see the binary genotype-phenotype mapping")
        # print(intactbigenphe)
    

    bikoinshapss, koinrows = calKSV4Input(intactbigenphe, inputnames, outputnames, False, table)
    print("-----Knockout Shapley value of input nodes in the binary network-----")
    for out, item in bikoinshapss.items():
        print(out, ":", item)

    bikiinshapss, kiinrows = calKSV4Input(intactbigenphe, inputnames, outputnames, True, table)
    print("-----KnockIN Shapley value of input nodes in the binary network-----")
    for out, item in bikiinshapss.items():
        print(out, ":", item) 

    # show the vanilla binary network here  first 
    for outname in outputnames:
        showNetwork(binet, bikoinshapss[outname], bikiinshapss[outname], None, None, 'binarywithvalue.html')

    
    # ONLY TO TEST THE CONSISTENCY WITH THE ORIGINAL NETWORK 
    korows, kirows = None, None
    if isko:
        # now do the knockout procedure with the binary network, 
        print("-----Now do the knockout procedure with the binary network-----")

        # first generate all possible input states 
        inputstates = genInput(bispecies, inputnames, False)
        vsko = {} # to store genphe of knockout network 
        for internode in biinternames:
            if internode not in outputnames:
                print("Knockout {}".format(internode))
                inputstatescopy = copy.deepcopy(inputstates)
                kooutputs = []
                for inputstate in inputstatescopy:
                    output = getKnockoutOutput(biformulas, inputstate, [internode], True, 1000, False, extranodes)
                    kooutputs.append(output)

                genphe = extractPhe(inputnames, outputnames, kooutputs)
                vsko[internode] = genphe
        koshaps, korows = calKSV(intactbigenphe, vsko, outputnames, len(inputnames), table, inputnames)
        for outname in outputnames:
            print("---- Binary KNOCKOUT VALUE for output {}----".format(outname))
            print(dict(sorted(koshaps.items())))
            print("\n")
    if iski:
        print("-----Now do the knockIN procedure with the binary network-----")
        # now do the knockIN procedure with the binary network 
        vski = {}
        inputstates = genInput(bispecies, inputnames, False)
        for internode in biinternames:
            if internode not in outputnames:
                print("KnockIN {}".format(internode))
                inputstatescopy = copy.deepcopy(inputstates)
                kioutputs = []
                for inputstate in inputstatescopy:
                    output = getKnockoutOutput(biformulas, inputstate, [internode], True, 1000, False, extranodes, isKnockin=True)
                    kioutputs.append(output)
                genphe = extractPhe(inputnames, outputnames, kioutputs)
                vski[internode] = genphe

                # all_vars = sorted({var for inner in kioutputs for var in inner})
                # if internode == 'CI_p' or internode == 'ci':
                #     print("---- Binary KNOCKIN VALUE for node {}----".format(internode))
                #     for trow in kioutputs:
                #         prow = "".join(f"{var: <14}:{int(trow.get(var, False)): <2}" for var in all_vars)
                #         print(prow)

        kishaps, kirows = calKSV(intactbigenphe, vski, outputnames, len(inputnames), table, inputnames)
        for outname in outputnames:
            print("---- Binary KNOCKIN VALUE for output {}----".format(outname))
            print(dict(sorted(kishaps.items())))
            print("\n")
    
    propko, propki = None, None
    if isprop:        
        for outname in outputnames:
            print("-----Propageting to output {}-----".format(outname)) 
            diamonds = diamondDigger(binet, outname, biformulasdict) 
            finediamonds, dnomaids = refineDiamonds(binet, biformulasdict, diamonds)
            if debug:
                print("-----Diamonds found in the binary network-----")
                for begin, ends in diamonds.items():
                    print(f"Node {begin} has diamonds: {ends}")
            print("-----Flipped diamonds found in the binary network-----")
            for begin, ends in dnomaids.items():
                print(f"Node {begin} has flipped diamonds: {ends}")
            rowsofnodes = propagateWithDiamonds(binet, table, index, aindex, outname, biformulasdict, extranodes, finediamonds, dnomaids)
            # rowsofnodes = dict()
            # tobefilteredforeignersofnodes, goodforeignersofnodes, foreignersofnodes = preprocessBranching(binet, outname, biformulasdict)
            
            # # rowsofnodes = propagateBottomUp(binet, table, index, aindex, outname, biformulasdict, extranodes)
            # rowsofnodes = propageteBottomUpWithBranchedInformation(binet, table, index, aindex, outname, biformulasdict, extranodes, tobefilteredforeignersofnodes, goodforeignersofnodes, foreignersofnodes)

            # # correct the rows for input nodes 
            # for input in inputnames:
            #     rowsofnodes[input] = koinrows[input].union(kiinrows[input])


        ### Unroll the binary network for propagation
        # unrollnet = copy.deepcopy(binet) # make a copy of the binary network to unroll it later
        # for outname in outputnames:
        #     unrollednodes, unrolledformulas = unrollBinaryNetwork(unrollnet, outname, biformulasdictandorder, nodes_positions, debug)
        #     temrowsofnodes = propagateUnrollNetwork(unrollnet, table, outname, unrolledformulas, index, aindex, extranodes, unrollednodes)
        #     rowsofnodes = mergeUnrollNodes(unrollednodes, temrowsofnodes)
        
        # for node, rows in rowsofnodes.items():
            # print("Node {} counts rows: \n{}".format(node, sorted(list(rows))))
        ### End of unrolling and propagating

        propko, propki = rowstovalues(rowsofnodes, table, outname)
        for node, value in propko.items():
            propko[node] = round(propko[node], 4)
            propki[node] = round(propki[node], 4)
            print("{:20} \t\t\t KO: {:10} | KI: {:10}".format(node, propko[node], propki[node]))

    if isko and iski:
        errorko = 0.0
        errorki = 0.0 
        num = 0
        for outname in outputnames:
            showNetwork(binet, bikoinshapss[outname], bikiinshapss[outname], koshaps[outname], kishaps[outname], "binary.html")
            for node, tvalue in bikoinshapss[outname].items():
                try:
                    if abs(tvalue - propko[node]) >= 0.005 or abs(bikiinshapss[outname][node] - propki[node]) >= 0.005:
                        print("{:20} - Correct: KO: {:10} | KI: {:10} - Incorrect: KO: {:10} | KI: {:10}".format(node, tvalue, bikiinshapss[outname][node], propko[node], propki[node]))
                    errorko += (tvalue - propko[node])*(tvalue - propko[node])
                    errorki += (bikiinshapss[outname][node] - propki[node])*(bikiinshapss[outname][node] - propki[node])
                    num += 1
                except:
                    # print("Cannot find {} in the set".format(node))
                    continue
            for node, tvalue in koshaps[outname].items():
                try:
                    if abs(koshaps[outname][node] - propko[node]) >= 0.005 or abs(kishaps[outname][node] - propki[node]) >= 0.005:
                        print("{:20} - Correct: KO: {:10} | KI: {:10} - Incorrect: KO: {:10} | KI: {:10}".format(node, tvalue, kishaps[outname][node], propko[node], propki[node]))
                    errorko += (tvalue - propko[node])*(tvalue - propko[node])
                    errorki += (kishaps[outname][node] - propki[node])*(kishaps[outname][node] - propki[node])
                    num += 1
                except:
                    # print("Cannot find {} in the set".format(node))
                    continue
        
        # print("Number of nodes is {}".format(num))
        # print("Error KO: ", math.sqrt(errorko/num))
        # print("Error KI: ", math.sqrt(errorki/num))

        # also get the order of nodes to compare 

        calshapko = dict()
        calshapki = dict()
        for outname in outputnames:
            for node, tvalue in koshaps[outname].items():
                if node in propko:
                    calshapko[node] = tvalue
            for node, tvalue in kishaps[outname].items():
                if node in propko:
                    calshapki[node] = tvalue
            for node, tvalue in bikoinshapss[outname].items():
                if node in propko:
                    calshapko[node] = tvalue
            for node, tvalue in bikiinshapss[outname].items():
                if node in propko:
                    calshapki[node] = tvalue

        rankedcalko = rank_dict_values(calshapko)
        rankedcalki = rank_dict_values(calshapki)

        del propko[outname]
        del propki[outname]
        rankedko = rank_dict_values(propko)
        rankedki = rank_dict_values(propki)


        print("Ranked KO calculated: ", rankedcalko)
        print("Ranked KO prop: ", rankedko)
        print("Ranked KI calculated: ", rankedcalki)
        print("Ranked KI prop: ", rankedki)

        rankerrorko = 0
        rankerrorki = 0
        for node, rank in rankedko.items():
            if rank != rankedcalko[node]:
                rankerrorko += 1

        for node, rank in rankedki.items():
            if rank != rankedcalki[node]:
                rankerrorki += 1
        print("Rank error KO: {} over {}".format (rankerrorko, len(rankedko)))
        print("Rank error KI: {} over {}".format (rankerrorki, len(rankedki)))


        # print(korows)
        # print(kirows)
        # print(koinrows)
        # print(kiinrows)
        # print(rowsofnodes)
        if korows and kirows:
            for input in inputnames:
                if input in rowsofnodes: 
                    gt = koinrows[input].union(kiinrows[input])
                    lack = gt.difference(rowsofnodes[input])
                    extra = rowsofnodes[input].difference(gt)
                    # if len(lack) > 0 or len(extra) > 0:
                    if True:
                        print(input)
                        print("Lack  {:5}: {}".format(len(lack), sorted(list(lack))))
                        print("Extra {:5}: {}".format(len(extra), sorted(list(extra))))
                        print("Real  {:5}: {}".format(len(gt), sorted(list(gt))))
                        print("Prop  {:5}: {}\n".format(len(rowsofnodes[input]), sorted(list(rowsofnodes[input]))))

            for inter, row in korows.items():
                if inter in rowsofnodes:
                    gt = row.union(kirows[inter])
                    lack = gt.difference(rowsofnodes[inter])
                    extra = rowsofnodes[inter].difference(gt)
                    # if len(lack) > 0 or len(extra) > 0:
                    if True:
                        print(inter)
                        print("Lack  {:5}: {}".format(len(lack), sorted(list(lack))))
                        print("Extra {:5}: {}".format(len(extra), sorted(list(extra))))
                        print("Real  {:5}: {}".format(len(gt), sorted(list(gt))))
                        print("Prop  {:5}: {}\n".format(len(rowsofnodes[inter]), sorted(list(rowsofnodes[inter]))))

        # # # Get all variable names (assuming all inner dicts use the same keys)
        # all_vars = sorted({var for inner in table.values() for var in inner})

        # # Print header
        # # header = "".join(f"{var:<20}" for var in all_vars)
        # # print(f"{'ID':<10}{header}")

        # # Print rows
        # for id_, vars_dict in table.items():
        #     print(f"{id_:<10}")
        #     row = "".join(f"{var: <5}:{int(vars_dict.get(var, False)): <2}\t" for var in all_vars)
        #     print(f"{row}")
        #     print()
    return bidecimalpairs


             
            

