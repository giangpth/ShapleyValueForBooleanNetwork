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
           

def enrichfilter(op, nodetofilter, onenode, index, aindex, extranodes, countedchild, allrows, indirect=False):
    print("---FILTERING---")
    if op == 'OR':
        if indirect:
            print("INDRIRECT, do not filter")
            pass
        else:
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
    
                             
def processBranching(net, toconvergenode, carryon, countedrowsofnodes, rowsofpairs, resofpairs, index, aindex, formulas, extranodes, convergedpairs, simtable):
    if toconvergenode in carryon:
        des = nx.descendants(net, toconvergenode)
        desincluded = copy.deepcopy(des)
        desincluded.add(toconvergenode)
        print("BRANCHING: Node {} has branching having desendants {}".format(toconvergenode, desincluded))
        countedpairs = set()

        print("Checking convergence for the pairs {}".format(carryon[toconvergenode]))
        for pair in carryon[toconvergenode]:
            left, right = pair[0], pair[1] 
            if left in desincluded and right in desincluded:
                countedpairs.add(pair) 
                
        if len(countedpairs) > 0:
            print("Need convergence of the pairs {}".format(countedpairs))
        else:
            print("No pair to converge, continue")
            return
        for pair in countedpairs:
            print("PAIR {}: {}".format(pair, rowsofpairs[pair]))
            if toconvergenode not in convergedpairs:
                convergedpairs[toconvergenode] = set()
            if pair in convergedpairs[toconvergenode]:
                print("Already counted the pair {}, continue".format(pair))
                continue
            convergedpairs[toconvergenode].add(pair)
            # carryon[toconvergenode].remove(pair) # remove the pair that is already counted since any pair that is already converged to a node will not be carried on 
            ornode = resofpairs[pair] # get node resulted by the pair 
            countedchild = dict() # key are childs, value is the rows counted by the current node to this child 
            paths = nx.all_simple_paths(net, source=toconvergenode, target=ornode)

            us = set ()
            for path in paths:
                us.update(path)
            # print("To check set for node {} to {} is {}".format(toconvergenode, ornode, us))

            orrows = copy.deepcopy(rowsofpairs[pair]) # get rows counted by the pair, this time is only to pass to enrichfilter to test 

            tocheck = list(us) 
            # print("Initial to check list for node {} to {} is {}".format(toconvergenode, ornode, tocheck))
            while tocheck:
                onenode = tocheck.pop(0)
                try:
                    form = formulas[onenode]
                except:
                    print("Pass the input node {}".format(onenode))
                    continue
                op = form.val
                if op == 'OR' or op == 'AND':
                    if form.right.val in desincluded and form.left.val in desincluded:
                        if form.right.val not in us:
                            tocheck.append(form.right.val)
                            us.add(form.right.val)
                        if form.left.val not in us:
                            tocheck.append(form.left.val)
                            us.add(form.left.val)
                    elif form.right.val in desincluded and form.left.val not in desincluded:
                        nodepartner = form.left.val
                        orrows = enrichfilter(op, nodepartner, onenode, index, aindex, extranodes, countedchild, orrows)
                    elif form.right.val not in desincluded and form.left.val in desincluded:
                        nodepartner = form.right.val
                        orrows = enrichfilter(op, nodepartner, onenode, index, aindex, extranodes, countedchild, orrows)
                    else:
                        print("Both parents {} and {} of node {} are not in descendants".format(form.left.val, form.right.val, onenode))
                        
                else:
                    print("Unary operator, do not concern")
                    
            orrows = rowsofpairs[pair] # redeclare rows counted by the pair 
            print("All rows counted by the pair {} {}".format(pair, orrows))
            print("Need to intersect with the rows counted by the node {}: {} ".format(ornode, countedrowsofnodes[ornode]))
            orrows = orrows.intersection(countedrowsofnodes[ornode]) # visited[ornode] is the rows counted by the node resulted by the pair
            if countedchild:
                for child, rows in countedchild.items():
                    # print("Child {} has rows {}".format(child, rows))
                    orrows = orrows.intersection(rows)


            # from here, recheck rows that are filtered out since they may be counted
            rowsneedtocheck = rowsofpairs[pair].difference(orrows) 
            print("Rows left to check of pair {} are: \n {}".format(pair, rowsneedtocheck))
            for id in rowsneedtocheck:
                row = simtable[id] 
                if checkrow(net, row, toconvergenode, ornode, formulas):
                    orrows.add(id)

            # add this row to the value of the node  
            if toconvergenode not in countedrowsofnodes:
                countedrowsofnodes[toconvergenode] = set()
            countedrowsofnodes[toconvergenode].update(orrows)
        print("After BRANCHING, the row for node {} is {}".format(toconvergenode, countedrowsofnodes[toconvergenode]))
    else:
        print("No carryon pair to converge to node {}".format(toconvergenode))


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
                        
                                initrows = enrichfilter(op, form.left.val, node_, index, aindex, extranodes, node_, initrows)
                                if len(initrows) == 0:
                                    print("EMPTY ROWS, stop!")
                                    break
                            elif form.right.val not in desincluded and form.left.val in desincluded:
                                print('Parent {} of node {} is not in descendants'.format(form.right.val, node_))
                                if form.left.val not in checked:
                                    tocheck.append(form.left.val)
               
                                initrows = enrichfilter(op, form.right.val, node_, index, aindex, extranodes, node_, initrows)
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


# def getKnockoutOutput(formulas, inputstate, knockoutlist, isbi = False,  \
#                       maxStep=1000, debug=False, extranodes = None, isKnockin=False) -> dict:

# def knockoutSingleNode(formulas, inputstate, ):


def propagateBottomUp_v3(net, simtable, index, aindex, outname, formulas, extranodes):
    diamonds = diamondDigger(net, outname, formulas) 
    print("----------Diamond information:-----------")
    for node, diamond in diamonds.items():
        print("{}: {}".format(node, diamond))
    carryonofnodes = dict() 
    countedrowsofnodes = dict() 
    curs = [outname]
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

            # get rows that make current node cur counted 
            if cur not in countedrowsofnodes:
                print(f"Init all rows for current node {cur} since it inherits from noone")
                currows = set(range(len(simtable)))
                countedrowsofnodes[cur] = currows
            else: # count only rows that cur is counted
                currows = copy.deepcopy(countedrowsofnodes[cur]) 

            
            if cur in carryonofnodes:
                curry = carryonofnodes[cur]
            else:
                curry = dict()

            print("\nCounted row for current node {} is \n{}".format(cur, sorted(list(currows))))

            print("Carried on of current node {} is:".format(cur))
            for node_, carriedrows in curry.items():
                for side, rows in carriedrows.items():
                    print("{}-{}: {}".format(node_, side, sorted(list(rows))))

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
                if op == 'OR':
                    if not rightselfloop:
                        leftrows = aindex[form.right.val]
                    else:
                        leftrows = countedrowsofnodes[cur]

                    # rows that right counted are rows that left = False
                    if not leftselfloop:
                        rightrows = aindex[form.left.val]
                    else:
                        rightrows = countedrowsofnodes[cur]

                    carryrowofop = index[form.left.val].intersection(index[form.right.val]).intersection(currows)
                    

                elif op == 'AND':
                    leftrows = index[form.right.val]
                    # rows that right counted are rows that left = True
                    rightrows = index[form.left.val]

                    carryrowofop = aindex[form.left.val].intersection(aindex[form.right.val]).intersection(currows)
                else:
                    print("Do not support operator {}")
                
                # first assign the counted rows to left and right 
                if form.left.val not in countedrowsofnodes:
                    countedrowsofnodes[form.left.val] = set()
                countedrowsofnodes[form.left.val].update(leftrows.intersection(currows))
                if form.right.val not in countedrowsofnodes:
                    countedrowsofnodes[form.right.val] = set ()
                countedrowsofnodes[form.right.val].update(rightrows.intersection(currows))

                print("Rows assigned to node {} from operator {} are: \n{}".format(form.left.val, op, sorted(list(countedrowsofnodes[form.left.val]))))
                print("Rows assigned to node {} from operator {} are: \n{}".format(form.left.val, op, sorted(list(countedrowsofnodes[form.right.val])))) 

                '''
                # now process carry on 
                print("{} operator carried rows are: \n{}".format(op, sorted(list(carryrowofop))))
                if form.left.val not in carryonofnodes:
                    carryonofnodes[form.left.val] = dict()
                if form.right.val not in carryonofnodes:
                    carryonofnodes[form.right.val] = dict() 

                # first, both parties carry on rows of current operator, no filter 
                if cur not in carryonofnodes[form.left.val]:
                    carryonofnodes[form.left.val][cur] = dict()
                if cur not in carryonofnodes[form.right.val]:
                    carryonofnodes[form.right.val][cur] = dict()

                carryonofnodes[form.left.val][cur]['M'] = carryrowofop
                carryonofnodes[form.right.val][cur]['M'] = carryrowofop 

                # then now passing rows carried by current node cur, filter accordingly to the operator 
                print("Current node {} carries on:".format(cur))
                for node, sides in curry.items():
                    if node not in carryonofnodes[form.left.val]:
                        carryonofnodes[form.left.val][node] = dict()
                    if node not in carryonofnodes[form.right.val]:
                        carryonofnodes[form.right.val][node] = dict() 

                    for side, rows in sides.items():
                        print("{}-{}:{}".format(node, side, sorted(list(rows)))) 
                        lrows = rows.intersection(leftrows)
                        rrows = rows.intersection(rightrows)
                        mrows = rows.intersection(carryrowofop) 

                        if 'L' not in carryonofnodes[form.left.val][node]:
                            carryonofnodes[form.left.val][node]['L'] = set()
                        carryonofnodes[form.left.val][node]['L'].update(lrows) 
                        if 'M' not in carryonofnodes[form.left.val][node]:
                            carryonofnodes[form.left.val][node]['M'] = set()
                        carryonofnodes[form.left.val][node]['M'].update(mrows)

                        if 'R' not in carryonofnodes[form.right.val][node]:
                            carryonofnodes[form.right.val][node]['R'] = set()
                        carryonofnodes[form.right.val][node]['R'].update(rrows)
                        if 'M' not in carryonofnodes[form.right.val][node]:
                            carryonofnodes[form.right.val][node]['M'] = set()
                        carryonofnodes[form.right.val][node]['M'].update(mrows)

                # now merge, left first  
                if form.left.val in diamonds:
                    for node, siderows in carryonofnodes[form.left.val].items(): 
                        if node in diamonds[form.left.val]:
                            for side, rows in siderows.items():
                                if side == 'L':
                                    lrows = rows 
                                else:
                                    lrows = set ()
                                if side == 'R':
                                    rrows = rows
                                else:
                                    rrows = set ()
                                if side == 'M':
                                    mrows = rows
                                else:
                                    mrows = set()

                            rowstomerge = lrows.intersection(rrows)
                            rowstomerge = mrows.union(rowstomerge)
                            print("Need to merge rows below for node {} at node {}".format(node, form.left.val))
                            print(sorted(list(rowstomerge))) 
                            countedrowsofnodes[form.left.val].update(rowstomerge) 

                            # # now remove rows already merged from the carry on of the node 
                            # for side, rows in siderows.items(): 


                # now merge right 
                if form.right.val in diamonds:
                    for node, siderows in carryonofnodes[form.right.val].items(): 
                        if node in diamonds[form.right.val]:
                            for side, rows in siderows.items():
                                if side == 'L':
                                    lrows = rows 
                                else:
                                    lrows = set ()
                                if side == 'R':
                                    rrows = rows
                                else:
                                    rrows = set ()
                                if side == 'M':
                                    mrows = rows
                                else:
                                    mrows = set()

                            rowstomerge = lrows.intersection(rrows)
                            rowstomerge = mrows.union(rowstomerge)
                            print("Need to merge rows for node {} at node {}".format(node, form.right.val))
                            print(sorted(list(rowstomerge))) 
                            countedrowsofnodes[form.right.val].update(rowstomerge) 
                            
                '''
                
                if form.left.val not in nextlayer:
                    nextlayer.append(form.left.val)
                if form.right.val not in nextlayer: 
                    nextlayer.append(form.right.val)
                

            elif len(inedges) == 1:
                if form.val == 'NOT':
                    print(f"{cur} ==== NOT {form.right.val}")
                    singlemom = form.right.val
                else:
                    print(f"{cur} ==== {form.val}")
                    singlemom = form.val
                
                # first assign counted rows to singlemom 
                if singlemom not in countedrowsofnodes:
                    countedrowsofnodes[singlemom] = set()
                countedrowsofnodes[singlemom].update(currows)

                '''
                # now passing all the carry on of current node cur to its singlemom 
                if singlemom not in carryonofnodes:
                    carryonofnodes[singlemom] = dict()
                print("Current node {} carries on:".format(cur))
                for node, sides in curry.items():
                    if node not in carryonofnodes[singlemom]:
                        carryonofnodes[singlemom][node] = dict()

                    for side, rows in sides.items():
                        print("{}-{}:{}".format(node, side, sorted(list(rows)))) 
                        if side not in carryonofnodes[singlemom][node]:
                            carryonofnodes[singlemom][node][side] = set() 
                        carryonofnodes[singlemom][node][side].update(rows) 

                # now merge here 
                if singlemom in diamonds:
                    for node, siderows in carryonofnodes[singlemom].items(): 
                        if node in diamonds[singlemom]:
                            for side, rows in siderows.items():
                                if side == 'L':
                                    lrows = rows 
                                else:
                                    lrows = set ()
                                if side == 'R':
                                    rrows = rows
                                else:
                                    rrows = set ()
                                if side == 'M':
                                    mrows = rows
                                else:
                                    mrows = set()

                            rowstomerge = lrows.intersection(rrows)
                            rowstomerge = mrows.union(rowstomerge)
                            print("Need to merge rows for node {} at node {}".format(node, singlemom))
                            print(sorted(list(rowstomerge))) 
                            countedrowsofnodes[singlemom].update(rowstomerge) 
                '''

                if singlemom not in nextlayer:
                    nextlayer.append(singlemom)

            else:
                print(f"Encounter input node {cur}")
        curs = nextlayer
    return countedrowsofnodes




def propagateBottomUp_v2(net, simtable, index, aindex, outname, formulas, extranodes):
    carryonofnodes = dict() # carryonofnodes[node] = dict(), rows of pairs that are carried on to the next layer
    countedrowsofnodes = dict() # countedrowsofnodes[node] = set(), rows that make node count
    passingtimes = dict()
    passingby = dict()

    curs = [outname] # contain nodes at the same levels 
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

            # get rows that make cur counted 
            if cur not in countedrowsofnodes: # count all rows
                print(f"Init all rows for current node {cur} since it inherits from noone")
                countrows = set(range(len(simtable)))
                countedrowsofnodes[cur] = countrows
            else: # count only rows that cur is counted
                countrows = copy.deepcopy(countedrowsofnodes[cur])
            
            # get rows that carried on by curs (rows that need to be processed to merge)
            if cur in carryonofnodes:
                carryonofcur = carryonofnodes[cur] 
            else:
                carryonofcur = dict()
            
            if cur not in passingtimes:
                passingtimes[cur] = dict()
                
            print("\nCounted row for current node {} is \n{}".format(cur, sorted(list(countedrowsofnodes[cur]))))

            print("Carried on of current node {} is:".format(cur))
            for node_, carriedrows in carryonofcur.items():
                print("{}: {}".format(node_, sorted(list(carriedrows))))

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
                        leftrows = countedrowsofnodes[cur]

                    # rows that right counted are rows that left = False
                    if not leftselfloop:
                        rightrows = aindex[form.left.val]
                    else:
                        rightrows = countedrowsofnodes[cur]

                    carryrowofop = index[form.left.val].intersection(index[form.right.val]).intersection(countrows.union(carryonofcur))
                    
                elif op == 'AND':
                    leftrows = index[form.right.val]
                    # rows that right counted are rows that left = True
                    rightrows = index[form.left.val]

                    carryrowofop = aindex[form.left.val].intersection(aindex[form.right.val]).intersection(countrows.union(carryonofcur))
                else:
                    print("Do not support operator {}".format(op))

                print ("{} operator, carried rows are: \n{}".format(op, sorted(list(carryrowofop))))

                if form.left.val not in carryonofnodes:
                    carryonofnodes[form.left.val] = dict()

                if form.left.val not in passingtimes:
                    passingtimes[form.left.val] = dict()
                    # print(f"Init passing times for {form.left.val}")
                    # print(passingtimes[form.left.val])
                if form.left.val not in passingby:
                    passingby[form.left.val] = set()

                carryonofnodes[form.left.val][cur + '_LEFT_'] = carryrowofop # both parties carry on rows of current operator, no filter 
                passingtimes[form.left.val][cur] = 1
                passingby[form.left.val].add(cur)

                if form.right.val not in carryonofnodes:
                    carryonofnodes[form.right.val] = dict()

                if form.right.val not in passingtimes:
                    passingtimes[form.right.val] = dict() 
                    # print(f"Init passing times for {form.right.val}")
                    # print(passingtimes[form.right.val])
                if form.right.val not in passingby:
                    passingby[form.right.val] = set ()

                carryonofnodes[form.right.val][cur + "_RIGHT_"] = carryrowofop # both parties carry on rows of current operator, no filter
                passingtimes[form.right.val][cur] = 1 
                passingby[form.right.val].add(cur) 
                

                # now assign the counted rows to nodes
                leftrows = leftrows.intersection(countrows)
                if form.left.val not in countedrowsofnodes:
                    countedrowsofnodes[form.left.val] = set()
                countedrowsofnodes[form.left.val].update(leftrows)
                print("After operator {}, {} counts only rows: \n{}".format(op, form.left.val ,sorted(list(countedrowsofnodes[form.left.val]))))

                rightrows = rightrows.intersection(countrows)
                if form.right.val not in countedrowsofnodes:
                    countedrowsofnodes[form.right.val] = set()
                countedrowsofnodes[form.right.val].update(rightrows)
                print("After operator {}, {} counts only rows: \n{}".format(op, form.right.val ,sorted(list(countedrowsofnodes[form.right.val]))))

                # from here passing carried rows up 
                if cur in carryonofnodes:
                    for node_, carriedrows in carryonofnodes[cur].items():
                        # left first 
                        if node_ in carryonofnodes[form.left.val]:
                            if op == 'OR':
                                # leftcarryrows = 
                                carryonofnodes[form.left.val][node_].update(carriedrows.intersection(aindex[form.right.val]))

                            elif op == 'AND':
                                carryonofnodes[form.left.val][node_].update(carriedrows.intersection(index[form.right.val]))

                            if node_ in passingtimes[form.left.val]:
                                passingtimes[form.left.val][node_] += 1
                            else:
                                passingtimes[form.left.val][node_] = 1
                        else:
                            if op == 'OR':
                                carryonofnodes[form.left.val][node_] = (carriedrows.intersection(aindex[form.right.val]))
                            elif op == 'AND':
                                carryonofnodes[form.left.val][node_] = (carriedrows.intersection(index[form.right.val]))

                            passingtimes[form.left.val][node_] = 1

                        # now right 
                        if node_ in carryonofnodes[form.right.val]:
                            if op == 'OR':
                                carryonofnodes[form.right.val][node_].update(carriedrows.intersection(aindex[form.left.val]))
                            elif op == 'AND':
                                carryonofnodes[form.right.val][node_].update(carriedrows.intersection(index[form.left.val]))

                            if node_ in passingtimes[form.right.val]:
                                passingtimes[form.right.val][node_] += 1
                            else:
                                passingtimes[form.right.val][node_] = 1
                            
                        else:
                            if op == 'OR':
                                carryonofnodes[form.right.val][node_] = (carriedrows.intersection(aindex[form.left.val]))
                            elif op == 'AND':
                                carryonofnodes[form.right.val][node_] = (carriedrows.intersection(index[form.left.val]))

                            passingtimes[form.right.val][node_] = 1
                        print("Now passing time of {} is: {}".format(form.left.val, passingtimes[form.left.val]))
                        print("Now passing time of {} is: {}".format(form.right.val, passingtimes[form.right.val]))


                print("Passing by information:")
                print("{:20}: {}".format(form.left.val, passingby[form.left.val]))
                print("{:20}: {}".format(form.right.val, passingby[form.right.val]))


                # now merge here if possible 
                # left first 
                # if len(passingby[form.left.val]) > 1:
                #     for node_, times in passingtimes[form.left.val].items():
                #         print(f"Node {node_} passed by {form.left.val} {times} times")
                #         if times > 1:
                #             print("Merge rows carried for {} to node {}".format(node_, form.left.val))
                #             countedrowsofnodes[form.left.val].update(carryonofnodes[form.left.val][node_])
                #             print("Now rows of {} are: \n{}".format(form.left.val, sorted(list(countedrowsofnodes[form.left.val]))))
                # for node_, carriedrows in carryonofnodes:
                
                # now right 
                if len(passingby[form.right.val]) > 1:
                    for node_, times, in passingtimes[form.right.val].items():
                        print(f"Node {node_} passed by {form.right.val} {times} times")
                        if times > 1:
                            print("Merge rows carried for {} to node {}".format(node_, form.right.val))
                            countedrowsofnodes[form.right.val].update(carryonofnodes[form.right.val][node_])
                            print("Now rows of {} are: \n{}".format(form.right.val, sorted(list(countedrowsofnodes[form.right.val]))))


                print(f"After operator {op}, {form.left.val} carries on rows: ")
                for node_, carriedrows in carryonofnodes[form.left.val].items():
                    print("{:20}: {}".format(node_, sorted(list(carriedrows))))

                print(f"After operator {op}, {form.right.val} carries on rows: ")
                for node_, carriedrows in carryonofnodes[form.right.val].items():
                    print("{:20}: {}".format(node_, sorted(list(carriedrows)))) 


                # add parents to the next layer
                if form.left.val not in nextlayer:
                    nextlayer.append(form.left.val)
                if form.right.val not in nextlayer:
                    nextlayer.append(form.right.val)

            elif len(inedges) == 1:
                if form.val == 'NOT':
                    print(f"{cur} ==== NOT {form.right.val}")
                    singlemom = form.right.val
                else:
                    print(f"{cur} ==== {form.val}")
                    singlemom = form.val

                # first pass all the counted rows of current node ot its single mom since no filter needed 
                if singlemom not in countedrowsofnodes: 
                    countedrowsofnodes[singlemom] = dict()
                countedrowsofnodes[singlemom] = countrows
                print("Passing all the rows counted by current node {} to its single mom {}".format(cur, singlemom))
                print("After unary operator rows of {} are: \n{}".format(singlemom, sorted(list(countedrowsofnodes[singlemom]))))

                # now pass carryonrows of current node to its single parent 
                if singlemom not in carryonofnodes:
                    carryonofnodes[singlemom] = dict()
                
                if singlemom not in passingtimes:
                    passingtimes[singlemom] = dict()
                    # print(f"Init passing times for node {singlemom}")
                    # print(passingtimes[singlemom])
                if singlemom not in passingby:
                    passingby[singlemom] = set()

                passingby[singlemom].add(cur)

                if cur in carryonofnodes:
                    for node_, carriedrows in carryonofnodes[cur].items():
                        if node_ in carryonofnodes[singlemom]:
                            carryonofnodes[singlemom][node_].update(carriedrows)
                            passingtimes[singlemom][node_] += 1
                            print("Now passing time of {} is: {}".format(singlemom, passingtimes[singlemom]))
                        else:
                            carryonofnodes[singlemom][node_] = carriedrows
                            passingtimes[singlemom][node_] = 1
                            print("Now passing time of {} is: {}".format(singlemom, passingtimes[singlemom]))
                
                print(f"After unary operator {singlemom} carries on rows: ")
                for node_, carriedrows in carryonofnodes[singlemom].items():
                    print("{:20}: {}".format(node_, sorted(list(carriedrows))))

                print("Passing by information:")
                print("{:20}: {}".format(singlemom, passingby[singlemom]))
                
                # now do merging here 
                if len(passingby[singlemom]) > 1:
                    for node_, times, in passingtimes[singlemom].items():
                        if times > 1:
                            print("Merge rows carried from {} to node {}".format(node_, singlemom))
                            countedrowsofnodes[singlemom].update(carryonofnodes[singlemom][node_])
                            print("Now rows of {} are: \n{}".format(singlemom, sorted(list(countedrowsofnodes[singlemom]))))
                

                # add parents to the next layer
                if singlemom not in nextlayer:
                    nextlayer.append(singlemom)
            else:
                print("Reach node {} without in-comming edges".format(cur))
                continue
        
        print("Up until now counted rows for node to be process in the next layer are:")
        for node_ in nextlayer:
            print("{:20}: {}".format(node_, sorted(list(countedrowsofnodes[node_]))))

        curs = nextlayer
    print("Passing time of nodes are:")
    for node_, times in passingtimes.items():
        print("{:20}: {}".format(node_, times))
    return countedrowsofnodes 
                


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

                    androws = index[form.left.val].intersection(index[form.right.val])
                    pair = (form.left.val, form.right.val)
                    rowsofpairs[pair] = androws 
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
                    
                elif op == 'AND':
                    leftrows = index[form.right.val]
                    # rows that right counted are rows that left = True 
                    rightrows = index[form.left.val]

                    androws = aindex[form.left.val].intersection(aindex[form.right.val])
                    pair = (form.left.val, form.right.val)
                    rowsofpairs[pair] = androws
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
     
                else:
                    print(f"Do not support operator {op}")
                    break

                
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
                    processBranching_v2(net, form.left.val, carryonofnodes, countedrowsofnode, rowsofpairs, resofpairs, index, aindex, formulas, extranodes, converged, simtable)

                if len(list(net.out_edges(form.right.val))) > 1:
                    # print(list(net.out_edges(form.right.val)))
                    processBranching_v2(net, form.right.val, carryonofnodes, countedrowsofnode, rowsofpairs, resofpairs, index, aindex, formulas, extranodes, converged, simtable)
                
                # add left and right to the next layers 
                if form.right.val not in nextlayer:
                    nextlayer.append(form.right.val)
                if form.left.val not in nextlayer:  
                    nextlayer.append(form.left.val)

            elif len(inedges) == 1:
                if form.val == "NOT":
                    print(f"{cur} === NOT {form.right.val} inherits rows {sorted(list(countrows))}")
                    if cur in carryonofnodes:
                        carryonofnodes[form.right.val] = carryonofnodes[cur]

                    if form.right.val not in countedrowsofnode:
                        countedrowsofnode[form.right.val] = set()
                    countedrowsofnode[form.right.val].update(countrows)

                    # check branching
                    if len(list(net.out_edges(form.right.val))) > 1:
                        processBranching_v2(net, form.right.val, carryonofnodes, countedrowsofnode, rowsofpairs, resofpairs, index, aindex, formulas, extranodes, converged, simtable)
                   
                    # add this node to the next layer 
                    if form.right.val not in nextlayer:
                        nextlayer.append(form.right.val)
                else:
                    print(f"{cur} === {form.val} inherits rows {sorted(list(countrows))}") 
                    if cur in carryonofnodes:
                        carryonofnodes[form.val] = carryonofnodes[cur]
                    
                    if form.val not in countedrowsofnode:
                        countedrowsofnode[form.val] = set()
                    countedrowsofnode[form.val].update(countrows) 
                    
                    # check branching
                    if len(list(net.out_edges(form.val))) > 1:
                        processBranching_v2(net, form.val, carryonofnodes, countedrowsofnode, rowsofpairs, resofpairs, index, aindex, formulas, extranodes, converged, simtable)

                    # add new rows to the visited set 
                    if form.val not in nextlayer:
                        nextlayer.append(form.val)
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
            rowsofnodes = propagateBottomUp_v3(binet, table, index, aindex, outname, biformulasdict, extranodes)

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
