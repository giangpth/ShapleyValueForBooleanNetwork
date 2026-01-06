import copy 
import shap 
import numpy as np
import math 
import time
import random
import json

from pathlib import Path


from Shapley.visualization import showNetwork
from Shapley.utilities import parseArgumentsTest, parseArgumentsAnalysis, readfile,  getOrderedList, rowstovalues, toDecimal, genTableFromOutput, evaluation, print_ranked_keys, estimateErrorFromNodes, random_percentage_selection
from Shapley.shapleyCalculation import extractPhe, calKSV4Input, calKSV 
from Shapley.networkUtilities import diamondDigger, getOutput, getKnockoutOutput, simBinaryNetwork, findMasterDiamond, processMasterDiamonds, simulateOneNode, nodes_on_cycles_digraph, all_descendants
from Shapley.speciesHandler import getSpeciesName, getInputNames, genInput, initSpeciesStates
from Shapley.booleanFormulaHandler import parseFormula, toBinaryFormulas, getFormula
from Shapley.networkHandler import convertBiBooleanFormulas2Network, convertBooleanFormulas2Network, manipulateNetwork, rewireBinet
from benchmarker import parse_line

def BooleanShapleyAnalysis():
    parser = parseArgumentsAnalysis()
    args = vars(parser.parse_args())
    if len(args) < 2:
        parser.print_help()
    else:
        verbose = args['verbose']
        exp = args['expression']
        strout = args['output']
        mode = args['mode']
        outputname = strout

        print("The expression file is {}".format(exp))
        print("The interested output is: ")
        print(outputname)

        # now read the expression file and take the string formulas from this 
        lines = readfile(exp, verbose)

        speciesnames = getSpeciesName(lines, verbose) # get species names 
        print("List of all species:")
        print(speciesnames)
        
        if outputname not in speciesnames:
            print("The list of interested output is not a subset of species")
            return -1
        
        inputnames = getInputNames(lines, speciesnames, verbose)
        print("-----Input names-----")
        print(inputnames)
        internames = speciesnames.difference(inputnames).difference((outputname))
        print("----Intermediate nodes-----")
        print(internames)
        if verbose:
            print("----Getting boolean formulas-----")  
        strformulas = getFormula(lines, verbose) # list of formulas in strings
        formulas = dict()
        for strformula in strformulas:
            root = parseFormula(strformula, verbose)
            if verbose:
                print("Parsing formula for {}".format(strformula['left']))
                root.display()
                print("\n")
            formulas[strformula['left']] = root
        # get the network of the original model
        orinet, orifas = convertBooleanFormulas2Network(formulas, inputnames, \
            speciesnames, debug=verbose)
        print("Original network has {} nodes and {} edges".format(orinet.number_of_nodes(), orinet.number_of_edges()))
        if verbose:
            print("FEEBBACK ARCS IN THE ORIGINAL NETWORK:")
            print(orifas)

        sortedinput, sortedinter = getOrderedList(inputnames, internames, verbose)
        anet, aformulas, extranodes, nodes_positions = manipulateNetwork(orinet, inputnames, formulas, True, False, verbose)
        rowsofnodes, averightinrows, table, bidecimalpairs, protime = binarizeAndPropagate(formulas, aformulas, inputnames, outputname, \
                    speciesnames, "abinetwork", sortedinput, sortedinter, extranodes, mode, 15, totest = False, debug=verbose)
        print("Boolean Shapley Analysis completed")




def testEachModel(modpath, target, simprops):
    # now read the expression file and take the string formulas from this 
    lines = readfile(modpath, False)

    speciesnames = getSpeciesName(lines, False) # get species names 
    print("List of all species:")
    print(speciesnames)
    
    if target not in speciesnames:
        print("The list of interested output is not a subset of species")
        return -1
    
    inputnames = getInputNames(lines, speciesnames, False)
    print("-----Input names-----")
    print(inputnames)
    outname = target
    internames = speciesnames.difference(inputnames).difference((outname))
    print("----Intermediate nodes-----")
    print(internames)
    strformulas = getFormula(lines, False) # list of formulas in strings
    formulas = dict()
    for strformula in strformulas:
        root = parseFormula(strformula, False)
        formulas[strformula['left']] = root
    # get the network of the original model
    orinet, orifas = convertBooleanFormulas2Network(formulas, inputnames, \
        speciesnames, debug=False)
    print("Original network has {} nodes and {} edges".format(orinet.number_of_nodes(), orinet.number_of_edges()))

    # get the input states of the original network
    species = initSpeciesStates(speciesnames, False)
    inputstates = genInput(species, inputnames, False)

    onejob = dict()

    simrowsofnodes, oridecimalpairs, timesimori, timeori = workWithOriginalNetwork(orinet, inputnames, \
            inputstates, outname, internames, formulas, True, True, False)
    onejob['oritime'] = timeori

    onejob['propercentage'] = dict()
    for simprop in simprops:
        onejob['propercentage'][simprop] = dict()

        final_averightrows = 0
        final_protime = 0
        final_averightinrows = 0
        final_dncg = 0
        numsamples = 5
        for i in range(numsamples):
            sortedinput, sortedinter = getOrderedList(inputnames, internames, False)
            anet, aformulas, extranodes, nodes_positions = manipulateNetwork(orinet, inputnames, formulas, True, False, False)
            print(f"-----Now perform Boolean Shapley Analysis with simprop = {simprop} -----")
            proprows, averightinrows, table, bidecimalpairs, protime = binarizeAndPropagate(formulas, aformulas, inputnames, outname, \
                    speciesnames, "abinetwork", sortedinput, sortedinter, extranodes, 'Shapley', simprop, totest = True, debug=False)
            print(f"Boolean Shapley Analysis with {simprop} percentage of nodes to simulate completed")
            # print(f"Propagation time is: {protime}")
            averightrows, avg_dncg = evaluation(proprows, simrowsofnodes, table, outname)
            # print(f"Average number of correct rows: {averightrows}")
            final_averightrows += averightrows
            final_protime += protime
            final_averightinrows += averightinrows
            final_dncg += avg_dncg
        onejob['propercentage'][simprop]['protime'] = final_protime / numsamples
        onejob['propercentage'][simprop]['averightrows'] = final_averightrows / numsamples
        onejob['propercentage'][simprop]['averightinrows'] = final_averightinrows / numsamples
        onejob['propercentage'][simprop]['dncg'] = final_dncg / numsamples    
    return onejob

def percentageTest(inputfile, outputpath):
    # load csv file 
    txtpath = Path(inputfile)
    with txtpath.open("r", encoding="utf-8") as f:
        lines = f.readlines()

    jobs = []
    for ln in lines:
        # print(ln)
        parsed = parse_line(ln)
        if not parsed:
            continue
        model, targets = parsed
        for tgt in targets:
            jobs.append((model, tgt))
    print(jobs)
    perlist = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50]

    stat = dict()
    for job in jobs:
        print("\n\n========STARTING JOB FOR MODEL {} AND TARGET {} =========".format(job[0], job[1]))
        onejob =  testEachModel(job[0], job[1], perlist)
        stat[job] = onejob

    print("DONE ALL JOBS")
    print("START WRITING OUTPUT FILE")

    data_str_keys = {str(k): v for k, v in stat.items()}
    with open(outputpath, "w") as f:
        json.dump(data_str_keys, f, indent=4)



def BooleanShapleyTest():
    """
    Test function to perform Boolean Shapley Analysis and compare with the simulation results.
    """
    timeSimOriginalNetwork = 0 # only the time to simulate original network, lower bound 
    timePerformOriginalAnalysis = 0 # time to perform the Knockout/Knockin of the original network 
    timePropagation = 0 # time to propagate 
    timeHeading = 0 # time to read, parse formulas 
    
    bound = 0
    numinputs = 0
    numnnodes = 0
    numedges = 0

    parser = parseArgumentsTest()
    args= vars(parser.parse_args())
    if len(args) < 2:
        parser.print_help()
    else:
        time1 = time.time() 
        debug = args['debug']
        # isko = args['knockout']
        # iski = args['knockin']
        # isbi = args['binary']
        # isacyclic = args['acyclic']
        exp = args['expression']
        strout = args['output']
        mode = args['mode']
        outputname = strout

        print("The expression file is {}".format(exp))
        print("The interested output is: ")
        print(outputname)
        # if isko:
        #     print("Perform normal procedure for input and knockout procedure for intermediate nodes, mode: {}".format(mode))
        # else:
        #     print("Perform only normal procedure, mode: {}".format(mode))
        # now read the expression file and take the string formulas from this 
        lines = readfile(exp, debug)

        speciesnames = getSpeciesName(lines, debug) # get species names 
        print("List of all species:")
        print(speciesnames)
        numnnodes = len(speciesnames)

        if outputname not in speciesnames:
            print("The list of interested output is not a subset of species")
            return -1
        
        inputnames = getInputNames(lines, speciesnames, debug)
        print("-----Input names-----")
        print(inputnames)
        bound = len(inputnames) + 1
        numinputs = len(inputnames)


        internames = speciesnames.difference(inputnames).difference((outputname))
        print("----Intermediate nodes-----")
        print(internames) 


        if debug:
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
        
        numedges = orinet.number_of_edges()

        if debug:
            print("FEEBBACK ARCS IN THE ORIGINAL NETWORK:")
            print(orifas)

        
        # get the input states of the original network
        species = initSpeciesStates(speciesnames, debug)
        inputstates = genInput(species, inputnames, debug)
        
    
        time2 = time.time()
        timeHeading = time2 - time1 

        simrowsofnodes, oridecimalpairs, timesimori, timeori = workWithOriginalNetwork(orinet, inputnames, \
            inputstates, outputname, internames, formulas, True, mode=mode, debug=False) 

        timeSimOriginalNetwork += timesimori 
        timePerformOriginalAnalysis += timeori 
        
        '''
        # print ("-----Simulation output is:----------")
        # print(orioutputs)
        # for row in orioutputs:
        #     print(dict(sorted(row.items()))) 

        # for outputname in outputnames:
        #     workWithSHAP(list(inputnames), speciesnames, outputname, formulas, debug)
        '''


        print("----------WITH ACYCLIC FLAG-------------")
        # get the acyclic network 
        timeconver2acyclic = time.time()
        anet, aformulas, extranodes, nodes_positions = manipulateNetwork(orinet, inputnames, formulas, True, False, debug)
        timeconver2acyclic = time.time() - timeconver2acyclic

        showNetwork(anet, None, None, None, None, "acyclicnetwork.html")

        if debug:
            for term, formula in aformulas.items():
                print("{} = ".format(term))
                formula.display()

        # # now do the limiting procedure for each output node 
        # adecimalpairs, ashapss = workwithAcyclicNetwork(anet, \
        #     inputnames, internames, outputname, speciesnames, aformulas, extranodes, False, debug)
    
        
        # acount = 0
        # for oriinp, oriinter in oridecimalpairs.items():
        #     if oriinter != adecimalpairs[oriinp]:
        #         print(oriinp, oriinter, adecimalpairs[oriinp])
        #         acount += 1
        # print("\n----TOTAL NUMBER OF DIFFERENCES BETWEEN ORIGNAL NETWORK AND ACYCLIC NETWORK IS {}----\n".format(acount))
        

        # convert the acyclic network to binary network 
        print("-------------Now convert the acyclic network to the binary network---------------")
        sortedinput, sortedinter = getOrderedList(inputnames, internames, debug)

    
        prorowsofnodes, averightinrows, table,  bidecimalpairs, protime = binarizeAndPropagate(formulas, aformulas, inputnames, outputname, \
            speciesnames, "abinetwork", sortedinput, sortedinter, extranodes, mode, 15, totest = False, debug=debug)
        timePropagation += protime 
        # if 1:
        #     bicount = 0
        #     for ainp, ainter, in adecimalpairs.items():
        #         if ainter != bidecimalpairs[ainp]:
        #             print(ainp, ainter, bidecimalpairs[ainp])
        #             bicount += 1
        #     print("\n----TOTAL NUMBER OF DIFFERENCES BETWEEN ACYCLIC NETWORK AND BINARIZED NETWORK IS {}----\n".format(bicount))
        print(f"NUMBER OF NODES: {numnnodes}")
        print(f"NUMBER OF INPUTS: {numinputs}")
        print(f"NUMBER OF EDGES: {numedges}")
        print(f"TIME process heading: {timeHeading}")
        print(f"TIME simulate original network: {timeSimOriginalNetwork}")
        print(f"TIME perform original analysis process: {timePerformOriginalAnalysis}")
        print(f"TIME perform propagation: {timePropagation}")

        averightrows, ave_dncg = evaluation(prorowsofnodes, simrowsofnodes, table, outputname)
        

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



def workwithAcyclicNetwork(anet, inputnames, internames, outputname, speciesnames, aformulas, extranodes, iski, debug):
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
    blinkings = []
    for inputstate in cloneinputstates:
        output, blinking = getOutput(aformulas, inputstate, False, 1000, False, extranodes) 
        outputs.append(output)
        blinkings.append(blinking)
        inp, inter = toDecimal(output, sortedinput, sortedinter)
        decimalPairs[inp] = inter


    genphe = extractPhe(inputnames, outputname, outputs)
    koinshapss = calKSV4Input(genphe, inputnames, outputname)

    kiinshapss = None
    if iski: 
        kiinshapss = calKSV4Input(genphe, inputnames, outputname, True)

    
    print("----Input shapley value of acyclic network-----")
    for item in koinshapss.items():
        print(item)   
        print('\n')


    if kiinshapss:
        showNetwork(anet, koinshapss, kiinshapss, None, None, "acyclic.html")
    else:
        showNetwork(anet, koinshapss, kiinshapss, None, None, "acyclic.html") 

    return decimalPairs, koinshapss


def workWithOriginalNetwork(net, inputnames, inputstates, outputname, internames, \
                            formulas, totest, mode = 'Shapley', debug = False): 
    """
    Work with the original network to calculate Shapley values.
    Parameters:
        net: Original network
        inputnames: List of input node names
        speciesnames: List of all species names
        outputnames: List of output node names
        internames: List of intermediate node names
        formulas: Dictionary of boolean formulas for the original network
        totest: Boolean flag for simulating ko/ki to generate baseline for comparison or not
        mode: mode of calculation, can be 'Shapley' or 'Uniform', default is 'Shapley'
        debug: Boolean flag for debugging
    Returns:
        decimalPairs: Dictionary of input-output pairs in decimal representation
        koinputshapss: Dictionary of Shapley values for input nodes
    """
    koinrows, kiinrows, korows, kirows = None, None, None, None 
    time1 = time.time()
    cloneinputstates = copy.deepcopy(inputstates)

    sortedinput, sortedinter = getOrderedList(inputnames, internames, True)

    decimalPairs = dict()

    simoritimestart = time.time()
    outputs = []
    blingkings = []
    for inputstate in cloneinputstates:
        output, blinking = getOutput(formulas, inputstate, False, 1000, debug)        
        outputs.append(output) 
        blingkings.append(blinking) 
        inp, inter = toDecimal(output, sortedinput, sortedinter) 
        decimalPairs[inp] = inter 

    timesimori = time.time() - simoritimestart 

    simtable, index, aindex, bindex = genTableFromOutput(outputs, inputnames, blingkings, sortedinter, outputname, mode, False)

    genphe = extractPhe(inputnames, outputname, outputs) 
    koinputshapss, koinrows = calKSV4Input(genphe, inputnames, outputname, knockin = False, mode=mode, simtable=simtable, debug=False) 
    kiinputshapss, kiinrows = calKSV4Input(genphe, inputnames, outputname, knockin = True, mode=mode, simtable=simtable, debug=False) 

   
    print("======Time for SIMULATING ORIGINAL NETWORK : {} seconds===".format(timesimori))
    
    print("----Input knockout shapley value-----")
    for item in koinputshapss.items():
        print((item))
    
    print("----Input knockIN shapley value-----")
    for item in kiinputshapss.items():
        print((item))
   
    koshaps, kishaps = None, None
    rowsofnodes = None
    if totest:
        print("-------Now perform the knockout procedure to original network-------")
        # can use the inputstates 
        vs = {}
        for internode in internames:
            if internode not in outputname:
                print("Knockout {}".format(internode))
                clone2inputstates = copy.deepcopy(inputstates)
                koouputs = []
                koblinkings = []
                for inputstate in clone2inputstates:
                    output, koblinking = getKnockoutOutput(formulas, inputstate, [internode], False, 1000, False)
                    koouputs.append(output)
                    koblinkings.append(koblinking)
                kogenphe = extractPhe(inputnames, outputname, koouputs)
                vs[internode] = kogenphe 
        # for outname in outputnames:
        koshaps, korows = calKSV(genphe, vs, outputname, len(inputnames), mode=mode, simtable=simtable, inputnames=inputnames)
        print("----KNOCKOUT VALUE for output {}----".format(outputname))
        print(dict(sorted(koshaps.items())))
        print("\n")
        # print("counted rows for simluatated ko ")
        # print(korows)
    
        print("-------Now perform the KNOCKIN procedure to intermediate nodes-------")
        vski = {}
        for internode in internames:
            if internode not in outputname:
                print("Knockin {}".format(internode))
                clone3inputstates = copy.deepcopy(inputstates)
                kiouputs = []
                kiblinkings = []
                for inputstate in clone3inputstates:
                    output, kiblingking = getKnockoutOutput(formulas, inputstate, [internode],\
                                                False, 1000, False, None, True)
                    kiouputs.append(output)
                    kiblinkings.append(kiblingking)
                kigenphe = extractPhe(inputnames, outputname, kiouputs)
                vski[internode] = kigenphe 
        
        kishaps, kirows = calKSV(genphe, vski, outputname, len(inputnames), mode=mode, simtable=simtable, inputnames=inputnames)
        print("----KNOCKIN VALUE for output {}----".format(outputname))
        print(dict(sorted(kishaps.items())))
        print("\n")

        # print(f"koinrows {koinrows}")
        # print(f"kiinrows {kiinrows}")
        # print(f"korows {korows}")
        # print(f"kirows {kirows}")

        inputrows = {k: koinrows[k] | kiinrows[k] for k in koinrows}
        interrows = {k: korows[k] | kirows[k] for k in korows}


        rowsofnodes = inputrows | interrows

    timeori = time.time() - time1 
    print("======Time for FULLY SIMULATION ANALYSIS IS : {} seconds===".format(timeori))
    if 0: # koshaps and kishaps:
        showNetwork(net, koinputshapss, kiinputshapss, \
                koshaps, kishaps, "original.html")
    else:
        # showNetwork(net, koinputshapss[outname], kiinputshapss[outname], \
        #         None, None, "original.html")
        showNetwork(net, None, None, \
                None, None, "original.html")
     
    return rowsofnodes, decimalPairs, timesimori, timeori 
           

            
def propagateAcyclic(net, simtable, index, aindex, outname, formulas, masterDiamonds, inputnames, orderedBiNodes, debug=False):
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
    print("\n\n--------PROPAGATION WITH ACYCLIC NETWORK--------")
    curs = [outname]
    countedrowsofnode = dict() # countedrowsofnode[node] = set of rows that make node count
    rowsofsinknodes = dict() # rows that are about to be filter for a diamond 
    solddiamonds = set() # save the diamonds that are already converged for a node (key) to avoid doing it again
   
    while curs:
        if debug:
            print("\n\n---Processing layers of {}---".format(curs))
        nextlayer = []
        for cur in curs:
            if cur not in countedrowsofnode: # count all rows
                countrows = set(range(len(simtable)))
                countedrowsofnode[cur] = countrows
            else: # count only rows that cur is counted
                countrows = countedrowsofnode[cur]

            try:
                form = formulas[cur]
            except:
                if debug:
                    print("\nReach node {} without in-comming edges".format(cur))
                continue
            # first, check if cur is an extranode 
            

            # get incoming edge to cur node
            inedges = list(net.in_edges(cur))
            assert len(inedges) <= 2, print("Support only binary network")
            
            if debug:
                print("\nCounted row for current node {} is \n{}".format(cur, sorted(list(countrows))))
            if len(inedges) == 2: 
                if debug:
                    print(f"{cur} ==== {form.left.val} {form.val} {form.right.val}")
                
                # get operator
                op = form.val 

                # now start to process 
                if op == 'OR': 
                    # rows that left counted are rows that right = False 
                    leftrows = aindex[form.right.val]
                    # rows that right counted are rows that left = False 
                    rightrows = aindex[form.left.val]
                elif op == 'AND':
                    leftrows = index[form.right.val]
                    # rows that right counted are rows that left = True 
                    rightrows = index[form.left.val]
                else:
                    print(f"Do not support operator {op}")
                    break
                    
                # rowsofsinknodes[cur] = countedrowsofnode[cur] # save the rows of the sink node
                rowsofsinknodes[cur] = countrows
                
                if form.left.val not in inputnames:
                    if form.left.val not in countedrowsofnode:
                        countedrowsofnode[form.left.val] = set()

                    if form.left.val not in masterDiamonds:
                        # intersect with rows that cur is count 
                        leftrows = leftrows.intersection(countrows)
                        if debug:
                            print(f"After operator {op}, {form.left.val} count only rows: \n{sorted(list(leftrows))}")
                        countedrowsofnode[form.left.val].update(leftrows)
                    else: 
                        if debug:
                            print("PROCESS DIAMONDS for {}".format(form.left.val))
                        rowsfromdiamonds = processMasterDiamonds(net, simtable, form.left.val, masterDiamonds, rowsofsinknodes, index, aindex, formulas, solddiamonds, orderedBiNodes, debug)
                        # rowsfromdiamonds = processDiamond(net, simtable, form.left.val, masterDiamonds, diamonds, rowsofsinknodes, index, aindex, extranodes, formulas, solddiamonds)
                        countedrowsofnode[form.left.val].update(rowsfromdiamonds) 
               
                if form.right.val not in inputnames: # value for input is already calculated, no need to bother 
                    if form.right.val not in countedrowsofnode:
                        countedrowsofnode[form.right.val] = set()
                    
                    if form.right.val not in masterDiamonds:
                        # intersect with rows that cur is count 
                        rightrows = rightrows.intersection(countrows)
                        if debug:
                            print(f"After operator {op}, {form.right.val} count only rows: \n{sorted(list(rightrows))}")
                        countedrowsofnode[form.right.val].update(rightrows)
                    else:
                        if debug:
                            print("PROCESS DIAMONDS for {}".format(form.right.val))
                        rowsfromdiamonds = processMasterDiamonds(net, simtable, form.right.val, masterDiamonds, rowsofsinknodes, index, aindex, formulas, solddiamonds, orderedBiNodes, debug)
                        # rowsfromdiamonds = processDiamond(net, simtable, form.right.val, masterDiamonds, diamonds, rowsofsinknodes, index, aindex, extranodes, formulas, solddiamonds)
                        countedrowsofnode[form.right.val].update(rowsfromdiamonds)

                if form.left.val not in nextlayer: 
                    nextlayer.append(form.left.val)
                if form.right.val not in nextlayer:
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
                
                if singlemom not in countedrowsofnode:
                    countedrowsofnode[singlemom] = set()
                # check if current node cur is in a diamond of singlemom
                if singlemom not in masterDiamonds:
                    countedrowsofnode[singlemom].update(countrows)
                else:
                    if debug:
                        print("PROCESS DIAMONDS for {}".format(singlemom))
                    rowsfromdiamonds = processMasterDiamonds(net, simtable, singlemom, masterDiamonds, rowsofsinknodes, index, aindex, formulas, solddiamonds, orderedBiNodes, debug)
                    # rowsfromdiamonds = processDiamond(net, simtable, singlemom, masterDiamonds, diamonds, rowsofsinknodes, index, aindex, extranodes, formulas, solddiamonds)
                    countedrowsofnode[singlemom].update(rowsfromdiamonds)

                if singlemom not in nextlayer:
                    nextlayer.append(singlemom)
            if debug:
                print("\t\t\t---------------")
        curs = nextlayer
   
    return countedrowsofnode

def propagateCyclic(net, simtable, index, aindex, descendants, outname, formulas, simformulas, originalspecies, extranodes, masterDiamonds, inputnames, orderedBiNodes, simpercent = 15, debug=False):
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

    curs = [outname]
    countedrowsofnode = dict() # countedrowsofnode[node] = set of rows that make node count
    rowsofsinknodes = dict() # rows that are about to be filter for a diamond 
    solddiamonds = set() # save the diamonds that are already converged for a node (key) to avoid doing it again
    alreadysimulated = set() # save the nodes that are already simulated to avoid doing it again

    
    # potentiallylack = dict() # save the rows that are potentially lacky 
    # potentiallyextra = dict() # save the rows that are potentially extra 
    # encountered = set()
    while curs:
        if debug:
            print("\n\n---Processing layers of {}---".format(curs))
        nextlayer = []
        for cur in curs:
            if '_to_' in cur: # if cur is the node that is wired back from an extra node
                root, des = cur.split("_to_")[0], cur.split("_to_")[1]
                # print("\nCurrent node {} is an extra node, wire back to root node {}".format(cur, root))

                '''
                # inherit the potentially lack and extra from cur to root
                if root not in potentiallylack:
                    potentiallylack[root] = set()
                if cur in potentiallylack:
                    potentiallylack[root].update(potentiallylack[cur])
                if root not in potentiallyextra:
                    potentiallyextra[root] = set()
                if cur in potentiallyextra:
                    potentiallyextra[root].update(potentiallyextra[cur])
                '''
                if cur not in countedrowsofnode: # count all rows
                    countrows = set(range(len(simtable)))
                    countedrowsofnode[cur] = countrows  
                else: # count only rows that cur is counted

                    countrows = countedrowsofnode[cur]
                cur = root 
                if cur not in countedrowsofnode: # count all rows
                    countedrowsofnode[cur] = countrows
                else:
                    countedrowsofnode[cur].update(countrows)
            else:
                if cur not in countedrowsofnode: # count all rows
                    countrows = set(range(len(simtable)))
                    countedrowsofnode[cur] = countrows
                else: # count only rows that cur is counted
                    countrows = countedrowsofnode[cur]

            try:
                form = formulas[cur]
            except:
                if debug:
                    print("\nReach node {} without in-comming edges".format(cur))
                continue
            # first, check if cur is an extranode 
            

            # get incoming edge to cur node
            inedges = list(net.in_edges(cur))
            assert len(inedges) <= 2, print("Support only binary network")
            
            if debug:
                print("\nCounted row for current node {} is \n{}".format(cur, sorted(list(countrows))))
            if len(inedges) == 2: 
                if debug:
                    print(f"{cur} ==== {form.left.val} {form.val} {form.right.val}")
                
                '''
                leftdependent = False
                rightdependent = False
                doubledependent = False
                if form.left.val in descendants[form.right.val]:
                    if debug:
                        print(f"Left node {form.left.val} depends on right node {form.right.val}")
                    leftdependent = True
                if form.right.val in descendants[form.left.val]:
                    if debug:
                        print(f"Right node {form.right.val} depends on left node {form.left.val}")
                    if leftdependent:
                        doubledependent = True
                    else:
                        rightdependent = True

                # tocheckrows = countrows.union(potentiallylack.get(cur, set()))
                tocheckrows = set(range(len(simtable)))
                bothtrue, bothfalse, leftTrightF, leftFrightT = agreement(index, aindex, form.left.val, form.right.val, tocheckrows)
                if doubledependent:
                    if form.left.val not in potentiallylack:
                        potentiallylack[form.left.val] = set()
                    if form.right.val not in potentiallylack:
                        potentiallylack[form.right.val] = set()

                    if form.left.val not in potentiallyextra:
                        potentiallyextra[form.left.val] = set()
                    if form.right.val not in potentiallyextra:
                        potentiallyextra[form.right.val] = set()
                    if debug:
                        print("MUTALLY DEPENDENT nodes detected")
                    if form.val == 'OR': # if no cycles, for left takes right = F and for right takes left = F 
                        potentiallylack[form.left.val].update(bothtrue)
                        # potentiallylack[form.left.val].update(leftFrightT)
                        potentiallyextra[form.left.val].update(leftTrightF) 

                        potentiallylack[form.right.val].update(bothtrue)
                        # potentiallylack[form.right.val].update(leftTrightF)
                        potentiallyextra[form.right.val].update(leftFrightT)
                        
                    if form.val == 'AND': # if no cycles, for left takes right = T and for right takes left = T
                        potentiallylack[form.left.val].update(bothfalse)
                        potentiallyextra[form.left.val].update(leftFrightT)
                        
                        potentiallylack[form.right.val].update(bothfalse) 
                        potentiallyextra[form.right.val].update(leftTrightF)
                        
                elif rightdependent: # right node depends on left node
                    if form.left.val not in potentiallylack:
                        potentiallylack[form.left.val] = set()
                    if debug:
                        print("RIGHT node depends on LEFT node") 
                    potentiallylack[form.left.val].update(bothtrue)
                    if form.val == 'AND':
                        if form.left.val not in potentiallyextra:
                            potentiallyextra[form.left.val] = set()
                        potentiallyextra[form.left.val].update(leftFrightT)
                    
                        
                elif leftdependent: # left node depends on right node
                    if form.right.val not in potentiallylack:
                        potentiallylack[form.right.val] = set()
                    if debug:
                        print("LEFT node depends on RIGHT node")
                    potentiallylack[form.right.val].update(bothtrue)
                    if form.val == 'AND':
                        if form.right.val not in potentiallyextra:
                            potentiallyextra[form.right.val] = set()
                        potentiallyextra[form.right.val].update(leftTrightF)
                    
                # also update with the potential lack and extra from current node
                if cur in potentiallylack:
                    if form.left.val not in potentiallylack:
                        potentiallylack[form.left.val] = set()
                    if form.right.val not in potentiallylack:
                        potentiallylack[form.right.val] = set()
                    potentiallylack[form.left.val].update(potentiallylack[cur])
                    potentiallylack[form.right.val].update(potentiallylack[cur])    

                if cur in potentiallyextra:
                    if form.left.val not in potentiallyextra:
                        potentiallyextra[form.left.val] = set()
                    if form.right.val not in potentiallyextra:
                        potentiallyextra[form.right.val] = set()
                    potentiallyextra[form.left.val].update(potentiallyextra[cur])
                    potentiallyextra[form.right.val].update(potentiallyextra[cur])
                '''
                
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
                    
                    if form.left.val in selected and form.left.val not in alreadysimulated:
                        # check probablity to simulate this node 
                        alreadysimulated.add(form.left.val)
                        if debug:
                            print(f"Left node {form.left.val} is a node to SIMULATE, simulate it with all rows")
                        ressimrows = simulateOneNode(simtable, form.left.val, outname, simformulas, set(range(len(simtable))), extranodes, inputnames)
                        countedrowsofnode[form.left.val] = ressimrows

                    if form.left.val not in alreadysimulated:   
                        if form.left.val not in masterDiamonds:
                            # intersect with rows that cur is count 
                            leftrows = leftrows.intersection(countrows)
                            if debug:
                                print(f"After operator {op}, {form.left.val} count only rows: \n{sorted(list(leftrows))}")
                            countedrowsofnode[form.left.val].update(leftrows)
                        else: 
                            if debug:
                                print("PROCESS DIAMONDS for {}".format(form.left.val))
                            rowsfromdiamonds = processMasterDiamonds(net, simtable, form.left.val, masterDiamonds, rowsofsinknodes, index, aindex, formulas, solddiamonds, orderedBiNodes, debug)
                            # rowsfromdiamonds = processDiamond(net, simtable, form.left.val, masterDiamonds, diamonds, rowsofsinknodes, index, aindex, extranodes, formulas, solddiamonds)
                            countedrowsofnode[form.left.val].update(rowsfromdiamonds) 
               
                # if form.right.val not in inputnames: # value for input is already calculated, no need to bother 
                if 1: # propagate also to input nodes to estimate error 
                    if form.right.val not in countedrowsofnode:
                        countedrowsofnode[form.right.val] = set()
                    
                    if form.right.val in selected and form.right.val not in alreadysimulated:
                        # check probablity to simulate this node 
                        alreadysimulated.add(form.right.val)
                        if debug:
                            print(f"Right node {form.right.val} is a node to SIMULATE, simulate it with all rows")
                        ressimrows = simulateOneNode(simtable, form.right.val, outname, simformulas, set(range(len(simtable))), extranodes, inputnames)
                        countedrowsofnode[form.right.val] = ressimrows

                    if form.right.val not in alreadysimulated:
                        if form.right.val not in masterDiamonds:
                            # intersect with rows that cur is count 
                            rightrows = rightrows.intersection(countrows)
                            if debug:
                                print(f"After operator {op}, {form.right.val} count only rows: \n{sorted(list(rightrows))}")
                            countedrowsofnode[form.right.val].update(rightrows)
                        else:
                            if debug:
                                print("PROCESS DIAMONDS for {}".format(form.right.val))
                            rowsfromdiamonds = processMasterDiamonds(net, simtable, form.right.val, masterDiamonds, rowsofsinknodes, index, aindex, formulas, solddiamonds, orderedBiNodes, debug)
                            # rowsfromdiamonds = processDiamond(net, simtable, form.right.val, masterDiamonds, diamonds, rowsofsinknodes, index, aindex, extranodes, formulas, solddiamonds)
                            countedrowsofnode[form.right.val].update(rowsfromdiamonds)

                if form.left.val not in nextlayer and not leftroot: 
                    nextlayer.append(form.left.val)
                if form.right.val not in nextlayer and not rightroot:
                    nextlayer.append(form.right.val)

                if leftroot and (leftroot not in countedrowsofnode) and (leftroot not in nextlayer):
                    # print(f"Extranodes {form.left.val}, add {form.left.val} to nextlayer")
                    nextlayer.append(form.left.val)
                if rightroot and (rightroot not in countedrowsofnode) and (rightroot not in nextlayer):
                    # print(f"Extranodes {form.right.val}, add  {form.right.val} to nextlayer")
                    nextlayer.append(form.right.val)

                '''
                if debug:
                    print(f"Operator is {form.val}")
                    print(f"POTENTIALLY lacky rows for left node {form.left.val}: {sorted(list(potentiallylack.get(form.left.val, set()).difference(countedrowsofnode.get(form.left.val, set()))) )}")
                    print(f"POTENTIALLY lacky rows for right node {form.right.val}: {sorted(list(potentiallylack.get(form.right.val, set()).difference(countedrowsofnode.get(form.right.val, set()))))}")
                    print(f"POTENTIALLY extra rows for left node {form.left.val}: {sorted(list(countedrowsofnode.get(form.left.val, set()).intersection(potentiallyextra.get(form.left.val, set()))))}")
                    print(f"POTENTIALLY extra rows for right node {form.right.val}: {sorted(list(countedrowsofnode.get(form.right.val, set()).intersection(potentiallyextra.get(form.right.val, set()))))}")
                '''
            else:
                if form.val == 'NOT':
                    singlemom = form.right.val
                    if debug:
                        print(f"{cur} === NOT {singlemom}")
                else:
                    singlemom = form.val
                    if debug:
                        print(f"{cur} === {singlemom}")
                '''
                if singlemom not in potentiallylack:
                    potentiallylack[singlemom] = set()
                if singlemom not in potentiallyextra:
                    potentiallyextra[singlemom] = set() 
                # also update with the potential lack and extra from current node
                if cur in potentiallylack:
                    potentiallylack[singlemom].update(potentiallylack[cur])
                if cur in potentiallyextra:
                    potentiallyextra[singlemom].update(potentiallyextra[cur])
                '''
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
                        ressimrows = simulateOneNode(simtable, singlemom, outname, simformulas, set(range(len(simtable))), extranodes, inputnames)
                        countedrowsofnode[singlemom] = ressimrows
                    if singlemom not in alreadysimulated:
                        if singlemom not in countedrowsofnode:
                            countedrowsofnode[singlemom] = set()
                        # check if current node cur is in a diamond of singlemom
                        if singlemom not in masterDiamonds:
                            countedrowsofnode[singlemom].update(countrows)
                        else:
                            if debug:
                                print("PROCESS DIAMONDS for {}".format(singlemom))
                            rowsfromdiamonds = processMasterDiamonds(net, simtable, singlemom, masterDiamonds, rowsofsinknodes, index, aindex, formulas, solddiamonds, orderedBiNodes, debug)
                            # rowsfromdiamonds = processDiamond(net, simtable, singlemom, masterDiamonds, diamonds, rowsofsinknodes, index, aindex, extranodes, formulas, solddiamonds)
                            countedrowsofnode[singlemom].update(rowsfromdiamonds)

                if singlemom not in nextlayer and not singleroot:
                    nextlayer.append(singlemom)
                if singleroot and singleroot not in countedrowsofnode and singleroot not in nextlayer:
                    # print(f"Extranodes {singlemom}, add {singlemom} to nextlayer")
                    nextlayer.append(singlemom)
            if debug:
                print("\t\t\t---------------")
        curs = nextlayer
    print("SIMULATE {} nodes over {} nodes".format(len(alreadysimulated), len(originalspecies) - len(inputnames)- 1))
    return countedrowsofnode #, potentiallylack, potentiallyextra
'''
def propagate(net, simtable, index, aindex, outname, formulas, simformulas, extranodes, nodestosim, masterDiamonds, inputnames, orderedBiNodes, debug=False):
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
    print("\n\n--------PROPAGATION--------")
    curs = [outname]
    countedrowsofnode = dict() # countedrowsofnode[node] = set of rows that make node count
    rowsofsinknodes = dict() # rows that are about to be filter for a diamond 
    solddiamonds = set() # save the diamonds that are already converged for a node (key) to avoid doing it again
    alreadysimulated = set() # save the nodes that are already simulated to avoid doing it again
    while curs:
        if debug:
            print("\n\n---Processing layers of {}---".format(curs))
        nextlayer = []
        for cur in curs:
            try:
                form = formulas[cur]
            except:
                if debug:
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
            if debug:
                print("\nCounted row for current node {} is \n{}".format(cur, sorted(list(countrows))))
            if len(inedges) == 2: 
                if debug:
                    print(f"{cur} ==== {form.left.val} {form.val} {form.right.val}")
                # get operator
                op = form.val 
                # check if one of the node is extranode 
                leftselfloop, rightselfloop = False, False
                rootname = None
                if form.left.val in extranodes: 
                    if debug:
                        print(f"Left node {form.left.val} is an extra node")
                    # get the root node of the extra node 
                    rootname = form.left.val.split("_to_")[0]
                    desname = form.left.val.split("_to_")[1]
                    if rootname == desname:
                        leftselfloop = True

                if form.right.val in extranodes:
                    if debug:
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
                if form.left.val in nodestosim:
                    if form.left.val not in inputnames:
                        if form.left.val in alreadysimulated:
                            if debug:
                                print(f"Left node {form.left.val} is already simulated, skip it")
                        else:
                            alreadysimulated.add(form.left.val)
                            if debug:
                                print(f"Left node {form.left.val} is a node to simulate, simulate it with all rows of {cur}")
                            ressimrows = simulateOneNode(simtable, form.left.val, cur, simformulas, countedrowsofnode[cur], extranodes, inputnames)
                            countedrowsofnode[form.left.val] = ressimrows
                else:
                    if form.left.val not in inputnames:
                        if form.left.val not in countedrowsofnode:
                            countedrowsofnode[form.left.val] = set()
                        if form.left.val not in masterDiamonds:
                            # intersect with rows that cur is count 
                            leftrows = leftrows.intersection(countrows)
                            if debug:
                                print(f"After operator {op}, {form.left.val} count only rows: \n{sorted(list(leftrows))}")
                            countedrowsofnode[form.left.val].update(leftrows)
                        else: 
                            if debug:
                                print("PROCESS DIAMONDS for {}".format(form.left.val))
                            rowsfromdiamonds = processMasterDiamonds(net, simtable, form.left.val, masterDiamonds, rowsofsinknodes, index, aindex, extranodes, formulas, solddiamonds, orderedBiNodes, debug)
                            # rowsfromdiamonds = processDiamond(net, simtable, form.left.val, masterDiamonds, diamonds, rowsofsinknodes, index, aindex, extranodes, formulas, solddiamonds)
                            countedrowsofnode[form.left.val].update(rowsfromdiamonds) 
                
                # check right parent
                if form.right.val in nodestosim:
                    if form.right.val not in inputnames: # value for input is already calculated, no need to bother 
                        if form.right.val in alreadysimulated:
                            if debug:
                                print(f"Right node {form.right.val} is already simulated, skip it")
                        else:
                            alreadysimulated.add(form.right.val)
                            if debug:
                                print(f"Right node {form.right.val} is a node to simulate, simulate it with all rows of {cur}")
                            ressimrows = simulateOneNode(simtable, form.right.val, cur, simformulas, countedrowsofnode[cur], extranodes, inputnames)
                            countedrowsofnode[form.right.val] = ressimrows
                    
                else:
                    if form.right.val not in inputnames: # value for input is already calculated, no need to bother 
                        if form.right.val not in countedrowsofnode:
                            countedrowsofnode[form.right.val] = set()
                        if form.right.val not in masterDiamonds:
                            # intersect with rows that cur is count 
                            rightrows = rightrows.intersection(countrows)
                            if debug:
                                print(f"After operator {op}, {form.right.val} count only rows: \n{sorted(list(rightrows))}")
                            countedrowsofnode[form.right.val].update(rightrows)
                        else:
                            if debug:
                                print("PROCESS DIAMONDS for {}".format(form.right.val))
                            rowsfromdiamonds = processMasterDiamonds(net, simtable, form.right.val, masterDiamonds, rowsofsinknodes, index, aindex, extranodes, formulas, solddiamonds, orderedBiNodes, debug)
                            # rowsfromdiamonds = processDiamond(net, simtable, form.right.val, masterDiamonds, diamonds, rowsofsinknodes, index, aindex, extranodes, formulas, solddiamonds)
                            countedrowsofnode[form.right.val].update(rowsfromdiamonds)

                if form.left.val not in nextlayer:
                    nextlayer.append(form.left.val)
                if form.right.val not in nextlayer:
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
                
                
                if singlemom in nodestosim:
                    if singlemom not in inputnames:
                        if singlemom in alreadysimulated:
                            if debug:
                                print(f"Node {singlemom} is already simulated, skip it")
                        else:
                            alreadysimulated.add(singlemom)
                            if debug:
                                print(f"Node {singlemom} is a node to simulate, simulate it with all rows of {cur}")
                            ressimrows = simulateOneNode(simtable, singlemom, cur, simformulas, countedrowsofnode[cur], extranodes, inputnames)
                            countedrowsofnode[singlemom] = ressimrows
                else:
                    if singlemom not in inputnames: 
                        if singlemom not in countedrowsofnode:
                            countedrowsofnode[singlemom] = set()
                        # check if current node cur is in a diamond of singlemom
                        if singlemom not in masterDiamonds:
                            countedrowsofnode[singlemom].update(countrows)
                        else:
                            if debug:
                                print("PROCESS DIAMONDS for {}".format(singlemom))
                            rowsfromdiamonds = processMasterDiamonds(net, simtable, singlemom, masterDiamonds, rowsofsinknodes, index, aindex, extranodes, formulas, solddiamonds, orderedBiNodes, debug)
                            # rowsfromdiamonds = processDiamond(net, simtable, singlemom, masterDiamonds, diamonds, rowsofsinknodes, index, aindex, extranodes, formulas, solddiamonds)
                            countedrowsofnode[singlemom].update(rowsfromdiamonds)

                if singlemom not in nextlayer:
                    nextlayer.append(singlemom)
            if debug:
                print("\t\t\t---------------")
        curs = nextlayer
    return countedrowsofnode


def propagateWithDiamonds(net, simtable, index, aindex, outname, formulas, extranodes, masterDiamonds, debug=False):
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
    print("\n\n--------PROPAGATION--------")
    curs = [outname]
    countedrowsofnode = dict() # countedrowsofnode[node] = set of rows that make node count
    rowsofsinknodes = dict() # rows that are about to be filter for a diamond 
    solddiamonds = set() # save the diamonds that are already converged for a node (key) to avoid doing it again
    while curs:
        if debug:
            print("\n\n---Processing layers of {}---".format(curs))
        nextlayer = []
        for cur in curs:
            try:
                form = formulas[cur]
            except:
                if debug:
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
            if debug:
                print("\nCounted row for current node {} is \n{}".format(cur, sorted(list(countrows))))
            if len(inedges) == 2: 
                if debug:
                    print(f"{cur} ==== {form.left.val} {form.val} {form.right.val}")
                # get operator
                op = form.val 
                # check if one of the node is extranode 
                leftselfloop, rightselfloop = False, False
                rootname = None
                if form.left.val in extranodes:
                    if debug:
                        print(f"Left node {form.left.val} is an extra node")
                    # get the root node of the extra node 
                    rootname = form.left.val.split("_to_")[0]
                    desname = form.left.val.split("_to_")[1]
                    if rootname == desname:
                        leftselfloop = True

                if form.right.val in extranodes:
                    if debug:
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
                if form.left.val not in masterDiamonds:
                    # intersect with rows that cur is count 
                    leftrows = leftrows.intersection(countrows)
                    if debug:
                        print(f"After operator {op}, {form.left.val} count only rows: \n{sorted(list(leftrows))}")
                    countedrowsofnode[form.left.val].update(leftrows)
                else: 
                    if debug:
                        print("PROCESS DIAMONDS for {}".format(form.left.val))
                    rowsfromdiamonds = processMasterDiamonds(net, simtable, form.left.val, masterDiamonds, rowsofsinknodes, index, aindex, extranodes, formulas, solddiamonds)
                    # rowsfromdiamonds = processDiamond(net, simtable, form.left.val, masterDiamonds, diamonds, rowsofsinknodes, index, aindex, extranodes, formulas, solddiamonds)
                    countedrowsofnode[form.left.val].update(rowsfromdiamonds) 
                
                # check right parent
                if form.right.val not in countedrowsofnode:
                    countedrowsofnode[form.right.val] = set()
                if form.right.val not in masterDiamonds:
                    # intersect with rows that cur is count 
                    rightrows = rightrows.intersection(countrows)
                    if debug:
                        print(f"After operator {op}, {form.right.val} count only rows: \n{sorted(list(rightrows))}")
                    countedrowsofnode[form.right.val].update(rightrows)
                else:
                    if debug:
                        print("PROCESS DIAMONDS for {}".format(form.right.val))
                    rowsfromdiamonds = processMasterDiamonds(net, simtable, form.right.val, masterDiamonds, rowsofsinknodes, index, aindex, extranodes, formulas, solddiamonds)
                    # rowsfromdiamonds = processDiamond(net, simtable, form.right.val, masterDiamonds, diamonds, rowsofsinknodes, index, aindex, extranodes, formulas, solddiamonds)
                    countedrowsofnode[form.right.val].update(rowsfromdiamonds)

                if form.left.val not in nextlayer:
                    nextlayer.append(form.left.val)
                if form.right.val not in nextlayer:
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
                
                if singlemom not in countedrowsofnode:
                    countedrowsofnode[singlemom] = set()
                
                # check if current node cur is in a diamond of singlemom
                if singlemom not in masterDiamonds:
                    countedrowsofnode[singlemom].update(countrows)
                else:
                    if debug:
                        print("PROCESS DIAMONDS for {}".format(singlemom))
                    rowsfromdiamonds = processMasterDiamonds(net, simtable, singlemom, masterDiamonds, rowsofsinknodes, index, aindex, extranodes, formulas, solddiamonds)
                    # rowsfromdiamonds = processDiamond(net, simtable, singlemom, masterDiamonds, diamonds, rowsofsinknodes, index, aindex, extranodes, formulas, solddiamonds)
                    countedrowsofnode[singlemom].update(rowsfromdiamonds)

                if singlemom not in nextlayer:
                    nextlayer.append(singlemom)
            if debug:
                print("\t\t\t---------------")
        curs = nextlayer
    return countedrowsofnode
'''

# do everything with binary network (convert, get speciesnames, simulate...)          
def binarizeAndPropagate(oriformulas, aformulas, inputnames, outputname, orispeciesnames, networkname, sortedinput, sortedinter, extranodes, mode = 'Shapley', simpercent = 15, totest=True, debug=False):
    """
    Work with the binary network to calculate Shapley values and perform knockout/knockin analyses.
    Parameters:
        oriformulas: List of boolean formulas for the original network
        formulas: List of boolean formulas for the acyclic network derived from the original network 
        inputnames: List of input node names
        outputnames: List of output node names      
        orispeciesnames: Set of original species names
        networkname: Name of the network
        sortedinput: List of sorted input node names
        sortedinter: List of sorted intermediate node names
        extranodes: List of extra nodes added to the network 
        mode: Analysis mode, e.g., 'Shapley' or 'Uniform' 
        totest: Boolean flag for knockout and knockin test 
        debug: Boolean flag for debug mode
    Returns:
        None
    """
    print("\n\n-------WORKING WITH BINARY NETWORK-------")
    # print(aformulas)
    protime = 0
    time1 = time.time() 
    bistrformulas = toBinaryFormulas(aformulas, debug)
    # bistrformulas = toBinaryFormulas(oriformulas, debug)
    binet, nodes_positions = convertBiBooleanFormulas2Network(bistrformulas, inputnames, orispeciesnames, "bi" + networkname, False, debug, extranodes) 

    # unfinised work on handling cycles
    
    '''
    cycledbinet = rewireBinet(binet, extranodes) # rewire the binet to remove self-loop caused by extra nodes
    desc = all_descendants(cycledbinet, method='auto')
    '''
    # if True:
    #     for _node in desc:
    #         print(f"Node {_node} has descendants: \n {desc[_node]}")
    # showNetwork(cycledbinet, None, None, None, None, 'cyclebinet.html')
    # nodesincycles = nodes_on_cycles_digraph(cycledbinet)
    # nodetosimulate = set()
    # if nodesincycles:
    #     print("There are {} nodes in cycles after rewiring:".format(len(nodesincycles)))
    #     print(nodesincycles)
    #     relatednodes = set()
    #     for node in nodesincycles:
    #         # get list of node that point to this node
    #         preds = list(cycledbinet.predecessors(node))
    #         relatednodes.update(preds)
    #     print("Related nodes that point to nodes in cycles:")
    #     print(relatednodes)
    #     nodetosimulate = nodesincycles.union(relatednodes)
    #     print("Need to simulate {} nodes".format(len(nodetosimulate))) 
    protime += time.time() - time1 

    showNetwork(binet, None, None, None, None, 'binary.html')

    time2 = time.time()
    biformulas = [] # this list of dictionary is for the simulation, maintaining correctly the order for the simulation 
    biformulasdict = dict() # this is for the diamond digger 
    bispeciesnames = set() 
    # for term, bistrformula in bistrformulas.items():
    orderedBiNodes = dict()
    for i, formula in enumerate(bistrformulas):
        term = formula['term']
        orderedBiNodes[term] = i
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
        
    bispecies = initSpeciesStates(bispeciesnames, debug)

    biinternames = bispeciesnames.difference(inputnames).difference((outputname))
    if debug:
        print(f"----There are {len(biinternames)} intermediate nodes in binarized network:-----")
        print(biinternames)

    # simulate binary network, ONLY TO TEST THE CONSISTENCY WITH THE ORIGINAL NETWORK 
    # print("IN BINARY NETWORK")
    realbioutputs, blinkings, bidecimalpairs = simBinaryNetwork(biformulas, inputnames, bispeciesnames, sortedinput, sortedinter, False, 1000, extranodes)

    table, index, aindex, bindex = genTableFromOutput(realbioutputs, inputnames, blinkings, biinternames, outputname, mode, False) 

    if debug:
        print("BINDEX")
        for node, rows in dict(sorted(bindex.items())).items():
            print("Node {:20}: Rows {}".format(node, sorted(list(rows))))

    # print("-----Binary network simulation table-----")
    # for i, row in table.items():
    #     print("Row {:3}: {}".format(i, row['PROP']))

    intactbigenphe = extractPhe(inputnames, outputname, realbioutputs)
    if debug:
        print("--------Binary network maps between genotype and phenotype-------")
        print("Uncomment to see the binary genotype-phenotype mapping")
        # print(intactbigenphe)
    

    bikoinshapss, koinrows = calKSV4Input(intactbigenphe, inputnames, outputname, False, 'Shapley', table)

    bikiinshapss, kiinrows = calKSV4Input(intactbigenphe, inputnames, outputname, True, 'Shapley', table)
    print("-----Calculated value of input nodes-----")
    for out, item in dict(sorted(bikoinshapss.items())).items():
        # print(out, ":", item) 
        print("{:20} \t\t\t KO: {:10} | KI: {:10}".format(out, bikoinshapss[out], bikiinshapss[out]))

    inrows = {k: koinrows[k] | kiinrows[k] for k in koinrows}

    # show the vanilla binary network here  first 
    # for outname in outputnames:
    #     showNetwork(binet, bikoinshapss[outname], bikiinshapss[outname], None, None, 'binarywithvalue.html')

    protime += time.time() - time2 

    propko, propki = None, None    
        
    time3 = time.time() 
    print("=======PROPAGATE FOR OUTPUT {}=======".format(outputname)) 
    diamonds = diamondDigger(binet, outputname, biformulasdict, debug=False) 
    # finediamonds = refineDiamonds(binet, biformulasdict, diamonds)
    
    if debug:
        print("-----Diamonds found in the binary network-----")
        for begin, ends in diamonds.items():
            print(f"Node {begin} has diamonds: {ends}")

    masterDiamonds = findMasterDiamond(binet, diamonds) 
    if True:
        print("----Master diamonds-----")
        print(masterDiamonds)
    # rowsofnodes = propagateWithDiamonds(binet, table, index, aindex, outname, biformulasdict, extranodes, masterDiamonds)
    # need to take order of nodes to simulate the diamond and pass to the function
    # rowsofnodes = propagateBlinking(binet, table, index, aindex, bindex, outname, biformulasdict, biformulas, extranodes, nodetosimulate, masterDiamonds, inputnames, orderedBiNodes, debug=True)
    averightinrows = 1.0
    if extranodes:
        desc = None
        rowsofnodes = propagateCyclic(binet, table, index, aindex, \
                                    desc, outputname, biformulasdict, oriformulas, orispeciesnames, \
                                    extranodes, masterDiamonds, inputnames, \
                                    orderedBiNodes, simpercent=simpercent, debug=debug)
        print("IN CYCLIC NETWORK: COMPARE ROWS PROPAGTED TO INPUT NODES AND BASELINE TO ESTIMATE ERROR")
        averightinrows = estimateErrorFromNodes(inrows, rowsofnodes, inputnames)

    else:
        rowsofnodes = propagateAcyclic(binet, table, index, aindex, outputname, biformulasdict, masterDiamonds, inputnames, orderedBiNodes, debug=debug)
    
    # print(inrows)
    rowsofnodes = rowsofnodes | inrows 

    propko, propki = rowstovalues(rowsofnodes, table, outputname)
    protime += time.time() - time3

    print(f"AVE_CORRECT_ROWS_INPUT: {round(averightinrows, 4)}")
    print("=========================")

    print("PROPAGATED VALUES:")
    for node, value in dict(sorted(propko.items())).items():
        if node in orispeciesnames:
            propko[node] = round(propko[node], 4)
            propki[node] = round(propki[node], 4)
            print("{:20} \t\t\t KO: {:10} | KI: {:10}".format(node, propko[node], propki[node]))

    print("KO RANKING OF NODES:")
    print_ranked_keys(propko, oriformulas)
    print("=======")
    print("KI RANKING OF NODES:")
    print_ranked_keys(propki, orispeciesnames)
    print("=======")
    
    
    
    # # ==========================FROM HERE DOWN IS FOR TEST MODE ONLY===========================
    # # ONLY TO TEST THE CONSISTENCY WITH THE ORIGINAL NETWORK 
    # korows, kirows = None, None
    # averightrows = None
    # if totest:
    #     # now do the knockout procedure with the binary network, 
    #     print("-----Now do the knockout procedure with the binary network-----")

    #     # first generate all possible input states 
    #     inputstates = genInput(bispecies, inputnames, False)
    #     vsko = {} # to store genphe of knockout network 
    #     for internode in biinternames:
    #         if internode not in outputname:
    #             print("Knockout {}".format(internode))
    #             inputstatescopy = copy.deepcopy(inputstates)
    #             kooutputs = []
    #             koblinkings = []
    #             for inputstate in inputstatescopy:
    #                 output, koblinking = getKnockoutOutput(biformulas, inputstate, [internode], True, 1000, False, extranodes)
    #                 kooutputs.append(output)
    #                 koblinkings.append(koblinking)

    #             genphe = extractPhe(inputnames, outputname, kooutputs)
    #             vsko[internode] = genphe
    #     koshaps, korows = calKSV(intactbigenphe, vsko, outputname, len(inputnames), mode, table, inputnames)

    #     print("---- Binary KNOCKOUT VALUE for output {}----".format(outputname))
    #     print(dict(sorted(koshaps.items())))
    #     print("\n")

    #     print("-----Now do the knockIN procedure with the binary network-----")
    #     # now do the knockIN procedure with the binary network 
    #     vski = {}
    #     inputstates = genInput(bispecies, inputnames, False)
    #     for internode in biinternames: 
    #         if internode not in outputname:
    #             print("KnockIN {}".format(internode))
    #             inputstatescopy = copy.deepcopy(inputstates)
    #             kioutputs = [] 
    #             kiblinkings = []
    #             for inputstate in inputstatescopy: 
    #                 output, kiblinking = getKnockoutOutput(biformulas, inputstate, [internode], True, 1000, False, extranodes, isKnockin=True)
    #                 kioutputs.append(output)
    #                 kiblinkings.append(kiblinking)
    #             genphe = extractPhe(inputnames, outputname, kioutputs)
    #             vski[internode] = genphe
    #     kishaps, kirows = calKSV(intactbigenphe, vski, outputname, len(inputnames), mode, table, inputnames)

    #     print("---- Binary KNOCKIN VALUE for output {}----".format(outputname))
    #     print(dict(sorted(kishaps.items())))
    #     print("\n")
        
    #     if totest:
    #         averightrows = runmetric(bikoinshapss, bikiinshapss, koshaps, kishaps, outputname, \
    #                 propko, propki, binet, inputnames, korows, kirows, koinrows, kiinrows, \
    #                 rowsofnodes, orispeciesnames)
    return rowsofnodes, averightinrows, table, bidecimalpairs, protime

