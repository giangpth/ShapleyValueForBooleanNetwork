import math
import copy 
from booleanFormulaHandler import parseFormula, getResult, toBinaryFormulas, expandFunction, propagateFromTarget, expandFormula
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
            anet, aformulas, extranodes, nodes_layer = manipulateNetwork(orinet, inputnames, formulas, isacyclic, False, debug)
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
    print("---FILTERING---")
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
    
                             
def  processBranching_v2(net, toconvergenode, carryon, countedrowsofnodes, rowsofpairs, resofpairs, index, aindex, formulas, extranodes, convergedpairs, simtable):
    if toconvergenode in carryon:
        des = nx.descendants(net, toconvergenode)
        desincluded = copy.deepcopy(des)
        desincluded.add(toconvergenode)
        print("BRANCHING: Node {} has branching having desendants {}".format(toconvergenode, des))
        countedpairs = set()

        print("Checking convergence for the pairs {}".format(carryon[toconvergenode]))
        for pair in carryon[toconvergenode]:
            if pair[0] in desincluded and pair[1] in desincluded:
                countedpairs.add(pair) 

        if len(countedpairs) > 0:
            print("Need convergence of the pairs {}".format(countedpairs))
        else:
            print("No pair to converge, continue")
            return

        for pair in countedpairs:
            if toconvergenode not in convergedpairs:
                convergedpairs[toconvergenode] = set()
            if pair in convergedpairs[toconvergenode]:
                print("Already counted the pair {}, continue".format(pair))
                continue
            if not rowsofpairs[pair]:
                print(f"Pair {pair} carries no rows, continue")
                continue
            convergedpairs[toconvergenode].add(pair)
            rowsoffactors = dict()
            print("Rows carried by pair {} are {}".format(pair, sorted(list(rowsofpairs[pair]))))
            for factor in pair:
                print(f"----Process factor {factor}----")
                if factor == toconvergenode:
                    rowsoffactors[factor] = rowsofpairs[pair]
                    continue
                rowsoffactors[factor] = set()
                # get all the path from toconvergenode to the resulted node
                rowsofpath = dict() 
                paths = nx.all_simple_paths(net, source=toconvergenode, target=factor)
                for id_, path in enumerate(paths):
                    initrows = rowsofpairs[pair] # get rows counted by the pair, this time is only to pass to enrichfilter to test
                    print("Path {}: {}".format(id_, path))
                    checked = set()
                    tocheck = list(path)
                    # tofilter = set() 
                    while tocheck:
                        node_ = tocheck.pop(0)
                        if node_ in checked:
                            continue
                        checked.add(node_)
                        if node_ == toconvergenode:
                            continue
                        try:
                            form = formulas[node_]
                        except:
                            print("Pass the input node {}".format(node_))
                            continue
                        op = form.val 
                        if op == 'OR' or op == 'AND':
                            if form.right.val in desincluded and form.left.val in desincluded:
                                if form.right.val not in checked:
                                    tocheck.append(form.right.val)
                                if form.left.val not in checked:
                                    tocheck.append(form.left.val)
                                if op == "OR":
                                    print('Both parents {} {} of {} in OR are in descendants, add to checklist'.format(form.right.val, form.left.val, node_))
                                continue
                            elif form.right.val in desincluded and form.left.val not in desincluded:
                                print('Parent {} of node {} is not in descendants'.format(form.left.val, node_))
                                if form.right.val not in checked:
                                    tocheck.append(form.right.val)
                        
                                initrows = filterrows(op, form.left.val, node_, index, aindex, extranodes, node_, initrows)
                                if len(initrows) == 0:
                                    print("EMPTY ROWS, stop!")
                                    break
                            elif form.right.val not in desincluded and form.left.val in desincluded:
                                print('Parent {} of node {} is not in descendants'.format(form.right.val, node_))
                                if form.left.val not in checked:
                                    tocheck.append(form.left.val)
               
                                initrows = filterrows(op, form.right.val, node_, index, aindex, extranodes, node_, initrows)
                                if len(initrows) == 0: 
                                    print("EMPTY ROWS, stop!")
                                    break
                            else:
                                print("BOTH parents {} and {} of node {} are NOT in DESCENDANTS".format(form.left.val, form.right.val, node_))
                        else:
                            print(f"Unary operator of {node_}, do not concern")
                    rowsofpath[id_] = initrows
                    # for node, op in tofilter:
                    #     print("Filter {} with operator {}".format(node, op))
                    #     initrows = enrichfilter(op, node, factor, index, aindex, extranodes, rowsoffactors[node], initrows)
                    #     print("Rows after filter are {}".format(sorted(list(initrows))))
                    # rowsofpath[id_] = initrows

                for id_, rows in rowsofpath.items():
                    print("Rows counted by path {} are {}".format(id_, sorted(list(rows))))
                    rowsoffactors[factor] = rowsoffactors[factor].union(rows)
            
            rowstoconverge = rowsoffactors[pair[0]].intersection(rowsoffactors[pair[1]]).intersection(countedrowsofnodes[resofpairs[pair]])
            # print("Rows counted by the pair {} are {}".format(pair, sorted(list(rowstoconverge))))

            rowstocheck = rowsofpairs[pair].difference(rowstoconverge)
            print("Rows left to check of pair {} are: \n {}".format(pair, rowstocheck))
            for id in rowstocheck:
                row = simtable[id] 
                if checkrow(net, row, toconvergenode, resofpairs[pair], formulas):
                    rowstoconverge.add(id)

            rowstoconverge = rowstoconverge.intersection(countedrowsofnodes[resofpairs[pair]])

            countedrowsofnodes[toconvergenode].update(rowstoconverge)
            print("--------")
    else:
        print("Carry nothing to converge")

def  processBranching(net, toconvergenode, carryon, countedrowsofnodes, rowsofpairs, resofpairs, index, aindex, formulas, extranodes, convergedpairs, simtable):
    if toconvergenode in carryon:
        des = nx.descendants(net, toconvergenode)
        desincluded = copy.deepcopy(des)
        desincluded.add(toconvergenode)
        print("BRANCHING: Node {} has branching having desendants {}".format(toconvergenode, des))
        countedpairs = set()
        

        print("Checking convergence for the pairs {}".format(carryon[toconvergenode]))
        for pair in carryon[toconvergenode]:
            if pair[0] in desincluded and pair[1] in desincluded:
                countedpairs.add(pair) 

        if len(countedpairs) > 0:
            print("Need convergence of the pairs {}".format(countedpairs))
        else:
            print("No pair to converge, continue")
            return

        for pair in countedpairs:
            if toconvergenode not in convergedpairs:
                convergedpairs[toconvergenode] = set()
            if pair in convergedpairs[toconvergenode]:
                print("Already counted the pair {}, continue".format(pair))
                continue
            if not rowsofpairs[pair]:
                print(f"Pair {pair} carries no rows, continue")
                continue
            convergedpairs[toconvergenode].add(pair)
            
            opofpair = formulas[resofpairs[pair]].val

            print("\nRows carried by pair {} are {}".format(pair, sorted(list(rowsofpairs[pair]))))
            if rowsofpairs[pair] == set():
                print(f"Pair {pair} carries no rows, continue")
                continue

            nodeandop = dict() # save the node and operator to be filtered for each factor 
            for factor in pair:
                print(f"----Process factor {factor}----")
                # for each factor, scan all the simple paths from factor to target node and save the node and operator to be filtered 
                nodeandop[factor] = dict()
                if factor == toconvergenode:
                    continue 

                paths = nx.all_simple_paths(net, source=toconvergenode, target=factor)

                for id_, path in enumerate(paths):
                    print("Path {}: {}".format(id_, path))
                    checked = set()
                    tocheck = list(path)
                    # tofilter = set() 
                    while tocheck:
                        node_ = tocheck.pop(0)
                        if node_ in checked:
                            continue
                        checked.add(node_)
                        if node_ == toconvergenode:
                            continue
                        try:
                            form = formulas[node_]
                        except:
                            print("Pass the input node {}".format(node_))
                            continue
                        op = form.val 

                        if op == 'OR' or op == 'AND':
                            if form.right.val in desincluded and form.left.val in desincluded:
                                if form.right.val not in checked:
                                    tocheck.append(form.right.val)
                                if form.left.val not in checked:
                                    tocheck.append(form.left.val)

                                if op == "OR":
                                    print('Both parents {} {} of {} in OR are in descendants, add to checklist'.format(form.right.val, form.left.val, node_))
                            elif form.right.val in desincluded and form.left.val not in desincluded:
                                print('Parent {} of node {} is not in descendants'.format(form.left.val, node_))
                                if form.right.val not in checked:
                                    tocheck.append(form.right.val)
                        
                                nodeandop[factor][form.left.val] = op
                            elif form.right.val not in desincluded and form.left.val in desincluded:
                                print('Parent {} of node {} is not in descendants'.format(form.right.val, node_))
                                if form.left.val not in checked:
                                    tocheck.append(form.left.val)
               
                                nodeandop[factor][form.right.val] = op
                            else:
                                print("BOTH parents {} and {} of node {} are NOT in DESCENDANTS".format(form.left.val, form.right.val, node_))
                        else:
                            print(f"Unary operator of {node_}, do not concern")
            
            if opofpair == 'OR':
                rowstoprocess = copy.deepcopy(rowsofpairs[pair])
                andset1 = set() 
                andset2 = set()
                for nodetofilter, op in nodeandop[pair[0]].items():
                    if op == 'AND':
                        andset1.add(nodetofilter)
                    elif op == 'OR':
                        rowstoprocess = filterrows(op, nodetofilter, pair[0], index, aindex, extranodes, rowstoprocess)
                for nodetofilter, op in nodeandop[pair[1]].items(): 
                    if op == 'AND':
                        andset2.add(nodetofilter)
                    elif op == 'OR':
                        rowstoprocess = filterrows(op, nodetofilter, pair[1], index, aindex, extranodes, rowstoprocess)

                if andset1 and andset2:
                    print("AND set 1 is {}".format(andset1))
                    print("AND set 2 is {}".format(andset2))
                    rowofandset1 = copy.deepcopy(rowstoprocess)
                    rowofandset2 = copy.deepcopy(rowstoprocess)
                    for node in andset1: 
                        rowofandset1 = filterrows('AND', node, toconvergenode, index, aindex, extranodes, rowofandset1)
                    for node in andset2:
                        rowofandset2 = filterrows('AND', node, toconvergenode, index, aindex, extranodes, rowofandset2)
                    rowstoprocess = rowofandset1.union(rowofandset2)

                countedrowsofnodes[toconvergenode].update(rowstoprocess)

            if opofpair == 'AND':
                rowstoprocess = copy.deepcopy(rowsofpairs[pair])
                orset1 = set()
                orset2 = set()
                for nodetofilter, op in nodeandop[pair[0]].items():
                    if op == 'OR':
                        orset1.add(nodetofilter)
                    elif op == 'AND':
                        rowstoprocess = filterrows(op, nodetofilter, pair[0], index, aindex, extranodes, rowstoprocess)
                for nodetofilter, op in nodeandop[pair[1]].items(): 
                    if op == 'OR':
                        orset2.add(nodetofilter)
                    elif op == 'AND':
                        rowstoprocess = filterrows(op, nodetofilter, pair[1], index, aindex, extranodes,  rowstoprocess)

                if orset1 and orset2:
                    print("OR set 1 is {}".format(orset1))
                    print("OR set 2 is {}".format(orset2))
                    orset = orset1.union(orset2)
                    for node in orset:
                        rowstoprocess = filterrows('OR', node, toconvergenode, index, aindex, extranodes, rowstoprocess)
                countedrowsofnodes[toconvergenode].update(rowstoprocess)
                    
    else:
        print("Carry nothing to converge")


def propagateBottomUp(net, simtable, index, aindex, outname, formulas, extranodes):
    curs = [outname] 
    countedrowsofnode = dict() # visited[node] = set(), rows that make node count 
    rowsofpairs = dict() # save the loss information of the OR operator 
    carryonofnodes = dict() # carryonofnodes[node] = set(), pairs that are carried on to the next layer
    resofpairs = dict() # the node that is the result of pairs 
    converged = dict() # save the pairs that are already converged for a node (key) to avoid doing it again 
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
                        
                    # get the rows that the root node is counted 

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
                        # leftrows = set()

                    rowsofthispair = index[form.left.val].intersection(index[form.right.val]).intersection(countrows)
                    
                elif op == 'AND':
                    leftrows = index[form.right.val]
                    # rows that right counted are rows that left = True 
                    rightrows = index[form.left.val]

                    rowsofthispair = aindex[form.left.val].intersection(aindex[form.right.val]).intersection(countrows)
     
                else:
                    print(f"Do not support operator {op}")
                    break
                    
                pair = (form.left.val, form.right.val)
                rowsofpairs[pair] = rowsofthispair
                resofpairs[pair] = cur
                print(f"Pair {pair} has rows {sorted(list(rowsofpairs[pair]))}")

                if form.left.val not in carryonofnodes:
                    carryonofnodes[form.left.val] = {pair}
                else:
                    carryonofnodes[form.left.val].add(pair)
                if cur in carryonofnodes:
                    carryonofnodes[form.left.val] = carryonofnodes[form.left.val].union(carryonofnodes[cur])
                if form.right.val not in carryonofnodes:
                    carryonofnodes[form.right.val] = {pair}
                else:
                    carryonofnodes[form.right.val].add(pair)
                if cur in carryonofnodes:
                    carryonofnodes[form.right.val] = carryonofnodes[form.right.val].union(carryonofnodes[cur])

                
                #left first
                # intersect with rows that cur is count 
                leftrows = leftrows.intersection(countrows)
                print(f"After operator {op}, {form.left.val} count only rows {sorted(list(leftrows))}")
                if form.left.val not in countedrowsofnode:
                    countedrowsofnode[form.left.val] = set()
                countedrowsofnode[form.left.val].update(leftrows)

                # intersect with rows that cur is count 
                rightrows = rightrows.intersection(countrows)
                print(f"After operator {op}, {form.right.val} count only rows {sorted(list(rightrows))}")
                if form.right.val not in countedrowsofnode:
                    countedrowsofnode[form.right.val] = set()
                countedrowsofnode[form.right.val].update(rightrows)
                
                # do the convergence here if needed 
                if len(list(net.out_edges(form.left.val))) > 1:
                    # print(list(net.out_edges(form.left.val)))
                    processBranching(net, form.left.val, carryonofnodes, countedrowsofnode, rowsofpairs, resofpairs, index, aindex, formulas, extranodes, converged, simtable)

                if len(list(net.out_edges(form.right.val))) > 1:
                    # print(list(net.out_edges(form.right.val)))
                    processBranching(net, form.right.val, carryonofnodes, countedrowsofnode, rowsofpairs, resofpairs, index, aindex, formulas, extranodes, converged, simtable)
                
                # add left and right to the next layers 
                if form.right.val not in nextlayer:
                    nextlayer.append(form.right.val)
                if form.left.val not in nextlayer:  
                    nextlayer.append(form.left.val)

            elif len(inedges) == 1:
                if form.val == "NOT":
                    singlemom = form.right.val 
                    print(f"{cur} === NOT {singlemom} inherits rows {sorted(list(countrows))}")
                else:
                    singlemom = form.val
                    print(f"{cur} === {singlemom} inherits rows {sorted(list(countrows))}")
                
                if cur in carryonofnodes:
                    carryonofnodes[singlemom] = carryonofnodes[cur]

                if singlemom not in countedrowsofnode:
                    countedrowsofnode[singlemom] = set()
                countedrowsofnode[singlemom].update(countrows)

                # check branching
                if len(list(net.out_edges(singlemom))) > 1:
                    processBranching(net, singlemom, carryonofnodes, countedrowsofnode, rowsofpairs, resofpairs, index, aindex, formulas, extranodes, converged, simtable)
                
                # add this node to the next layer 
                if singlemom not in nextlayer:
                    nextlayer.append(singlemom)
            else:
                print("Reach node {} without in-comming edges".format(cur))
                continue 
        curs = nextlayer
    return countedrowsofnode 

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
                
def preprocessBranching (binet, outname, biformulas):
    print("-------Preprocessing branching for binary network-------")
    reversenet = binet.reverse(copy=True) # reverse the network to get the in-comming edges
    relevantnodes = nx.descendants(reversenet, outname) # get all the nodes that are relevant to the output node
    relevantnodes.add(outname) # add the output node to the relevant nodes

    foreignersofnodes = dict()
    halfforeignersofnodes = dict() # nodes that have one parent is foreign, the other is not 
    tobeconvergedofnodes = dict() # save the nodes that both parents are in the family 
    foreingerop = dict()

    for node in relevantnodes:
        numoutedges = len(list(binet.out_edges(node)))
        if numoutedges > 1: # possibly branched node
            foreignersofnodes[node] = set() # save the foreign nodes that have child inside the family 
            halfforeignersofnodes[node] = set() # save the foreign nodes that have one parent inside the family, the other is not
            tobeconvergedofnodes[node] = set() # save the nodes that need to be converged
            foreingerop[node] = dict()

            des = nx.descendants(binet, node) # get all the descendants of the node
            desinclude = copy.deepcopy(des)
            desinclude.add(node) 
            middlenodes = relevantnodes.intersection(desinclude) # get the middle nodes lying between source and target 

            for mnode in middlenodes:
                if mnode not in biformulas:
                    print("Node {} is not in the binary formulas, pass".format(mnode))
                    continue
                if mnode == node: 
                    print("Node {} is the source node, do nothing".format(mnode))
                    continue

                form = biformulas[mnode]

                if form.val == 'OR' or form.val == 'AND':
                    # check if both parents are in descendants 
                    if form.left.val in desinclude and form.right.val in desinclude:
                        print("Node {} has branching with parents {} and {} are in the family".format(mnode, form.left.val, form.right.val))
                        tobeconvergedofnodes[node].add(mnode) # add the node to the set of nodes that need to be converged
                    elif form.left.val in desinclude and form.right.val not in desinclude:
                        print("Node {} has branching with parent {} is in the family, but parent {} is not".format(mnode, form.left.val, form.right.val))
                        # add the foreign node to the set of foreigners
                        foreignersofnodes[node].add(form.right.val)
                        foreingerop[node][form.right.val] = form.val # save the operator of the foreign node

                        halfforeignersofnodes[node].add(mnode) # add the foreign node to the set of half foreigners 
                        
                    elif form.left.val not in desinclude and form.right.val in desinclude:
                        print("Node {} has branching with parent {} is in the family, but parent {} is not".format(mnode, form.right.val, form.left.val))
                        # add the foreign node to the set of foreigners
                        foreignersofnodes[node].add(form.left.val)
                        foreingerop[node][form.left.val] = form.val # save the operator of the foreign node

                        halfforeignersofnodes[node].add(mnode)
                    else:
                        print("Node {} has no branching with parents {} and {} are both not in the family, IMPOSSIBLE".format(mnode, form.left.val, form.right.val))
                else:
                    print("Node {} is unary operator, do nothing".format(mnode)) 


    print("-------Demographics information-------")
    for node, foreigners in foreignersofnodes.items():
        if len(foreigners) > 0:
            print("Node {} has foreigners {}".format(node, foreigners))
        print('\n')

        if len(halfforeignersofnodes[node]) > 0:
            print("Node {} has half foreigners {}".format(node, halfforeignersofnodes[node]))
        print('\n') 

        if len(tobeconvergedofnodes[node]) > 0:
            print("Node {} has to be converged {}".format(node, tobeconvergedofnodes[node]))
        print('\n\n')

    # from here find foreigner nodes that do not need to be filtered out  
    goodforeignersofnodes = dict() # save the foreign nodes that do not need to be filtered out 
    tobefilteredforeignersofnodes = dict() # save the foreign nodes that need to be filtered out

    for node, tobeconvergednodes in tobeconvergedofnodes.items():
        print('Processing node {}'.format(node))
        goodforeignersofnodes[node] = set() # initialize the set of good foreigners for each node
        tobefilteredforeignersofnodes[node] = dict() # initialize the set of to be filtered foreigners for each node
        if len(tobeconvergednodes) == 0:
            continue
        
        for tobeconvergednode in tobeconvergednodes: 
            tobefilteredforeignersofnodes[node][tobeconvergednode] = dict() # initialize the set of to be filtered foreigners for each to be converged node
            convergeform = biformulas[tobeconvergednode]
            left = convergeform.left.val
            right = convergeform.right.val
            op = convergeform.val
            print('Converging node {} ==== {} {} {}'.format(tobeconvergednode, left, op, right))

            leftdes = nx.descendants(reversenet, left)
            rightdes = nx.descendants(reversenet, right)

            leftforeigners = foreignersofnodes[node].intersection(leftdes) 
            print("Left foreigners of node {} to node {} are {}".format(node, tobeconvergednode, leftforeigners))
          
            rightforeigners = foreignersofnodes[node].intersection(rightdes) 
            print("Right foreigners of node {} to node {} are {}".format(node, tobeconvergednode, rightforeigners))

            leftor = set() 
            rightor = set() 
            leftand = set()
            rightand = set()
            for foreigner in leftforeigners: 
                if foreigner in foreingerop[node]:
                    if foreingerop[node][foreigner] == 'OR':
                        leftor.add(foreigner)
                    elif foreingerop[node][foreigner] == 'AND':
                        leftand.add(foreigner)
            for foreigner in rightforeigners:
                if foreigner in foreingerop[node]:
                    if foreingerop[node][foreigner] == 'OR':
                        rightor.add(foreigner)
                    elif foreingerop[node][foreigner] == 'AND':
                        rightand.add(foreigner)

            tobefilteredforeignersofnodes[node][tobeconvergednode]['leftor'] = leftor
            tobefilteredforeignersofnodes[node][tobeconvergednode]['rightor'] = rightor
            tobefilteredforeignersofnodes[node][tobeconvergednode]['leftand'] = leftand
            tobefilteredforeignersofnodes[node][tobeconvergednode]['rightand'] = rightand

            if op == 'AND': 
                if leftor and rightor:
                    pass 
                else:
                    if leftor:
                        goodforeignersofnodes[node].update(leftor)
                        print("Add leftor {} to good foreigners of node {}".format(leftor, node))
                    elif rightor:
                        goodforeignersofnodes[node].update(rightor)
                        print("Add rightor {} to good foreigners of node {}".format(rightor, node))
                    
            
            if op == 'OR':
                if leftand and rightand:
                    pass 
                else:
                    if leftand:
                        goodforeignersofnodes[node].update(leftand)
                        print("Add leftand {} to good foreigners of node {}".format(leftand, node))
                    elif rightand:
                        goodforeignersofnodes[node].update(rightand)
                        print("Add rightand {} to good foreigners of node {}".format(rightand, node))
        print('\n')
    
    print("-------Good foreigners information-------")
    for node, goodforeigners in goodforeignersofnodes.items():
        if len(goodforeigners) > 0:
            print("Node {} has good foreigners {}".format(node, goodforeigners))
    
    print("-------To be filtered foreigners information-------")
    for node, tobefilteredforeigners in tobefilteredforeignersofnodes.items():
        if len(tobefilteredforeigners) > 0:
            print("Node {} has to be filtered foreigners {}".format(node, tobefilteredforeigners))

    return tobefilteredforeignersofnodes, goodforeignersofnodes 

def convergediamond(node, tobefilteredofnodes, rowsofsinknodes, formulas, index, aindex, extranodes, goodforeignersofnodes, countedrowsofnode):
    if node in tobefilteredofnodes:
        print(f"Node {node} has branching and needs to be filtered out, process it")
    else:
        print(f"Node {node} has no branching, do nothing")
        return
    
    sinks = tobefilteredofnodes[node] 
    for sink, dicts in sinks.items():
        if sink not in rowsofsinknodes:
            continue
        
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
                if node_ not in goodforeignersofnodes[node]:
                    allrows = filterrows('OR', node_, sink, index, aindex, extranodes, allrows)
            leftrows = copy.deepcopy(allrows)
            rightrows = copy.deepcopy(allrows) 
            
            for node_ in leftand:
                if node_ not in goodforeignersofnodes[node]:
                    leftrows = filterrows('AND', node_, sink, index, aindex, extranodes, leftrows)
            for node_ in rightand:
                if node_ not in goodforeignersofnodes[node]:
                    rightrows = filterrows('AND', node_, sink, index, aindex, extranodes, rightrows)

            countedrowsofnode[node].update(leftrows.union(rightrows))
                    
        elif sinkop == 'AND':
            # filter out all the AND operators
            # for OR op, if belong to only one side then no filter 
            for node_ in leftand.union(rightand):
                if node_ not in goodforeignersofnodes[node]:
                    allrows = filterrows('AND', node_, sink, index, aindex, extranodes, allrows)

            if leftor and rightor:
                print("Both leftor and rightor are not empty, need to filter")
                for node_ in leftor.union(rightor):
                    if node_ not in goodforeignersofnodes[node]:
                        allrows = filterrows('OR', node_, sink, index, aindex, extranodes, allrows)

            countedrowsofnode[node].update(allrows) 
        else:
            print("Dot not support operator {}".format(sinkop))
        
        

def propageteBottomUpWithBranchedInformation(net, simtable, index, aindex, outname, formulas, extranodes, tobefilteredofnodes, goodforeignersofnodes):
    curs = [outname] 
    countedrowsofnode = dict() # countedrowsofnode[node] = set of rows that make node count 
    rowsofsinknodes = dict() # rows that are about to be filter for a diamond 
    converged = dict() # save the nodes that are already converged for a node (key) to avoid doing it again 
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
                        
                    # get the rows that the root node is counted 

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
                        # leftrows = set()

                    rowsofthispair = index[form.left.val].intersection(index[form.right.val]).intersection(countrows)
                    
                elif op == 'AND':
                    leftrows = index[form.right.val]
                    # rows that right counted are rows that left = True 
                    rightrows = index[form.left.val]

                    rowsofthispair = aindex[form.left.val].intersection(aindex[form.right.val]).intersection(countrows)
     
                else:
                    print(f"Do not support operator {op}")
                    break
                    
                rowsofsinknodes[cur] = rowsofthispair # save the rows of the sink node

                
                #left first
                # intersect with rows that cur is count 
                leftrows = leftrows.intersection(countrows)
                print(f"After operator {op}, {form.left.val} count only rows {sorted(list(leftrows))}")
                if form.left.val not in countedrowsofnode:
                    countedrowsofnode[form.left.val] = set()
                countedrowsofnode[form.left.val].update(leftrows)

                # intersect with rows that cur is count 
                rightrows = rightrows.intersection(countrows)
                print(f"After operator {op}, {form.right.val} count only rows {sorted(list(rightrows))}")
                if form.right.val not in countedrowsofnode:
                    countedrowsofnode[form.right.val] = set()
                countedrowsofnode[form.right.val].update(rightrows)
                
                # do the convergence here if needed 
                if len(list(net.out_edges(form.left.val))) > 1: # left first
                    convergediamond(form.left.val, tobefilteredofnodes, rowsofsinknodes, formulas, index, aindex, extranodes, goodforeignersofnodes, countedrowsofnode)
                    
                    
                if len(list(net.out_edges(form.right.val))) > 1:
                    convergediamond(form.right.val, tobefilteredofnodes, rowsofsinknodes, formulas, index, aindex, extranodes, goodforeignersofnodes, countedrowsofnode)
                   
                   
                # add left and right to the next layers 
                if form.right.val not in nextlayer:
                    nextlayer.append(form.right.val)
                if form.left.val not in nextlayer:  
                    nextlayer.append(form.left.val)

            elif len(inedges) == 1:
                if form.val == "NOT":
                    singlemom = form.right.val 
                    print(f"{cur} === NOT {singlemom} inherits rows {sorted(list(countrows))}")
                else:
                    singlemom = form.val
                    print(f"{cur} === {singlemom} inherits rows {sorted(list(countrows))}")

                if singlemom not in countedrowsofnode:
                    countedrowsofnode[singlemom] = set()
                countedrowsofnode[singlemom].update(countrows)

                # check branching
                if len(list(net.out_edges(singlemom))) > 1:
                    convergediamond(singlemom, tobefilteredofnodes, rowsofsinknodes, formulas, index, aindex, extranodes, goodforeignersofnodes, countedrowsofnode)
                    
                # add this node to the next layer 
                if singlemom not in nextlayer:
                    nextlayer.append(singlemom)
            else:
                print("Reach node {} without in-comming edges".format(cur))
                continue 
        curs = nextlayer
    return countedrowsofnode 
    



# do everything with binary network (convert, get speciesnames, simulate...)          
def workwithBinaryNetwork(formulas, inputnames, outputnames, orispeciesnames, networkname, sortedinput, sortedinter, isko = False, iski = False, debug=False, extranodes=None, isprop=False):
    bistrformulas = toBinaryFormulas(formulas, True)
    binet, nodes_layer = convertBiBooleanFormulas2Network(bistrformulas, inputnames, orispeciesnames, "bi" + networkname, False, debug, extranodes) 
    
    biformulas = [] # this list of dictionary is for the simulation
    biformulasdict = dict() # this is for the expanding function 
    bispeciesnames = set() 
    # for term, bistrformula in bistrformulas.items():
    for formula in bistrformulas:
        term = formula['term']
        bistrformula = formula['formula']
        thisfor = {'left': term, "right": bistrformula}

        thisbiformula = dict()
        thisbiformula['term'] = term
        thisbiformula['formula'] = parseFormula(thisfor, debug) 
        biformulasdict[term] = thisbiformula['formula'] 

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
        showNetwork(binet, outname, bikoinshapss[outname], bikiinshapss[outname], None, None, 'binary.html')

    
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

                all_vars = sorted({var for inner in kioutputs for var in inner})
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
            tobefilteredforeignersofnodes, goodforeignersofnodes = preprocessBranching(binet, outname, biformulasdict)
            
            # rowsofnodes = propagateBottomUp(binet, table, index, aindex, outname, biformulasdict, extranodes)
            rowsofnodes = propageteBottomUpWithBranchedInformation(binet, table, index, aindex, outname, biformulasdict, extranodes, tobefilteredforeignersofnodes, goodforeignersofnodes)

            # # correct the rows for input nodes 
            # for input in inputnames:
            #     rowsofnodes[input] = koinrows[input].union(kiinrows[input])
    
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
        
        print("Number of nodes is {}".format(num))
        print("Error KO: ", errorko/num)
        print("Error KI: ", errorki/num)

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

        
        if korows and kirows:
            for input in inputnames:
                if input in rowsofnodes: 
                    gt = koinrows[input].union(kiinrows[input])
                    lack = gt.difference(rowsofnodes[input])
                    extra = rowsofnodes[input].difference(gt)
                    if len(lack) > 0 or len(extra) > 0:
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
