import math
import copy 
from booleanFormulaHandler import parseFormula, getResult, toBinaryFormulas, expandFunction, propagateFromTarget, expandFormula, unrollBinaryNetwork
from networkHandler import convertBiBooleanFormulas2Network, convertBooleanFormulas2Network, manipulateNetwork, showNetwork
import networkx as nx
from pyvis.network import Network 
import timeit
import shap 
import numpy as np
from utilities import parseArguments, readfile, getSpeciesName, getInputNames, getFormula, genInput, genSpecies, getOrderedList, getOutput, rowstovalues 
from utilities import toDecimal, extractPhe, calKSV4Input, getKnockoutOutput, simBinaryNetwork, genTableFromOutput, calKSV, diamondDigger, checkrow, rank_dict_values

def main():
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
            speciesnames, "network", debug=debug)
        
        
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
                inputnames, internames, outputnames, speciesnames, aformulas, extranodes, isko, isbi, iski, debug)
            
            acount = 0

            for oriinp, oriinter in oridecimalpairs.items():
                if oriinter != adecimalpairs[oriinp]:
                    print(oriinp, oriinter, adecimalpairs[oriinp])
                    acount += 1
            print("Total number of differences is {}".format(acount))

            if isbi: # convert the acyclic network to binary network 
                print("-------------Now conver the acyclic network to the binary network---------------")
                sortedinput, sortedinter = getOrderedList(inputnames, internames, True)
                bidecimalpairs = workwithBinaryNetwork(aformulas, inputnames, outputnames, \
                    speciesnames, "abinetwork", sortedinput, sortedinter, isko, iski, debug, extranodes, isprop)
                bicount = 0
                for ainp, ainter, in adecimalpairs.items():
                    if ainter != bidecimalpairs[ainp]:
                        print(ainp, ainter, bidecimalpairs[ainp])
                        bicount += 1
                print("Total number of differences between anetwork and binetwork is {}".format(bicount))


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
    print("----------Work with SHAP-----------")
    species = genSpecies(speciesnames, debug)

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



def workwithAcyclicNetwork(anet, inputnames, internames, outputnames, speciesnames, aformulas, extranodes, isko, isbi, iski, debug):
    species = genSpecies(speciesnames, debug)
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
            showNetwork(anet, outname, koinshapss[outname], kiinshapss[outname], None, None, "acyclic.html")
        else:
            showNetwork(anet, outname, koinshapss[outname], kiinshapss, None, None, "acyclic.html")
    return decimalPairs, koinshapss


def workWithOriginalNetwork(net, inputnames, speciesnames, outputnames, internames, \
                            formulas, isko, iski, debug): 
    species = genSpecies(speciesnames, debug)
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
            showNetwork(net, outname, koinputshapss[outname], kiinputshapss[outname], \
                    koshaps[outname], kishaps[outname], "original.html")
        else:
            showNetwork(net, outname, koinputshapss[outname], kiinputshapss[outname], \
                    None, None, "original.html")
    return decimalPairs, koinputshapss
           

def filterrows(op, nodetofilter, onenode, index, aindex, extranodes, allrows):
    print("FILTER")
    if op == 'OR':
        
        if nodetofilter not in extranodes:
            # print("NOT EXTRANODE")
            childrows = aindex[nodetofilter]
            # countedchild[onenode] = childrows
            print("Counting the rows that {} is FALSE for {} in operator {}".format(nodetofilter, onenode, op))
            allrows = allrows.intersection(childrows)
            print("Rows after filter are {}".format(sorted(list(allrows))))
        else:
            # print("EXTRANODE")
            coms = nodetofilter.split("_to_")
            if len(coms) == 2 and coms[0] == coms[1]:
                print("{} is selfloop, do not count it".format(nodetofilter))
            else:
                childrows = aindex[nodetofilter]
                # countedchild[onenode] = childrows
                print("Counting the rows that {} is FALSE for {} in operator {}".format(nodetofilter, onenode, op))
                allrows = allrows.intersection(childrows)
                print("Rows after filter are {}".format(sorted(list(allrows))))
    elif op == 'AND':
        if nodetofilter not in extranodes:
            # print("NOT EXTRANODE")
            childrows = index[nodetofilter]
            # countedchild[onenode] = childrows
            print("Counting the rows that {} is TRUE for {} in operator {}".format(nodetofilter, onenode, op))
            allrows = allrows.intersection(childrows)
            print("Rows after filter are {}".format(sorted(list(allrows))))
        else:
            # print("EXTRANODE")
            coms = nodetofilter.split("_to_")
            if len(coms) == 2 and coms[0] == coms[1]:
                print("{} is selfloop, do not count it".format(nodetofilter))
            else:
                childrows = index[nodetofilter]
                # countedchild[onenode] = childrows
                print("Counting the rows that {} is TRUE for {} in operator {}".format(nodetofilter, onenode, op))
                allrows = allrows.intersection(childrows)
                print("Rows after filter are {}".format(sorted(list(allrows))))
        # pass
    else:
        print("Do not support {} operator".format(op))
    
    return allrows
    

def topDown(net, formulas, burows, inrows, node_layers):   
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
                
def convergediamond(node, tobefilteredhalfforeigners, rowsofsinknodes, formulas, index, aindex, extranodes, goodhalfforeignersofnodes, countedrowsofnode, foreignerpartners):
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
        
def diamondDigger(binet, outname, biformulas):
    print("-------Digging diamonds-------")
    reversenet = binet.reverse(copy=True) # reverse the network to get the in-comming edges
    relevantnodes = nx.descendants(reversenet, outname) # get all the nodes that are relevant to the output node
    relevantnodes.add(outname) # add the output node to the relevant nodes
    processed = set() # save all the processed formulas 
    # relevantgraph = binet.subgraph(relevantnodes) # get the subgraph of the relevant nodes 
    relevantgraph = binet
    showNetwork(relevantgraph, outname, None, None, None, None, "relevantgraph.html")
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

def refineDiamonds(binet, biformulas, sinks):
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


def propagateUnrollNetwork(unrollednet, simtable, outname, unrolledformulas, index, aindex, extranodes, unrollednodes):
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

def mergeUnrollNodes(unrollednodes, rowsofnodes): 
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


# do everything with binary network (convert, get speciesnames, simulate...)          
def workwithBinaryNetwork(formulas, inputnames, outputnames, orispeciesnames, networkname, sortedinput, sortedinter, isko = False, iski = False, debug=False, extranodes=None, isprop=False):
    bistrformulas = toBinaryFormulas(formulas, True)
    binet, nodes_positions = convertBiBooleanFormulas2Network(bistrformulas, inputnames, orispeciesnames, "bi" + networkname, False, debug, extranodes) 
    
    for outname in outputnames:
        showNetwork(binet, outname, None, None, None, None, 'binary.html')


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
        biformulasdictandorder[term] = (forunroll, i) # save the order of the binary formulas

        biformulas.append(thisbiformula)

        bispeciesnames.add(term)
        coms = bistrformula.split()
        for com in coms:
            if com != '(' and com != ')' and com != 'AND' and com != 'OR' and com != 'NOT' and com != '=':
                if '-' in com: # this is for excluding '-' in the formulas which is confusing for the compiler 
                    com = com.replace('-','_')
                bispeciesnames.add(com)
        
    bispecies = genSpecies(bispeciesnames, debug)

    biinternames = bispeciesnames.difference(inputnames).difference(outputnames)
    if debug:
        print("----Intermediate nodes in binary network-----")
        print(biinternames)
        print("----Intermediate nodes size in binary network is {}-----".format(len(biinternames)))


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
        showNetwork(binet, outname, bikoinshapss[outname], bikiinshapss[outname], None, None, 'binarywithvalue.html')

    
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
            showNetwork(binet, outname, bikoinshapss[outname], bikiinshapss[outname], koshaps[outname], kishaps[outname], "binary.html")
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


             
            
if __name__ == "__main__":
    start = timeit.default_timer()
    main()
    print("-----------RUNNING TIME---------------")
    print(timeit.default_timer() - start)
