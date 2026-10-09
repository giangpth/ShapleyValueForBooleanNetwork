import time
import timeit
import copy
import networkx as nx
from Shapley.booleanFormulaHandler import toBinaryFormulas
from Shapley.networkHandler import convertBiBooleanFormulas2Network, manipulateNetwork
from Shapley.utilities import estimateErrorFromInputs, evaluation1, evaluation2, genTableFromOutput, getOrderedList, parseArgumentsAnalysis, compareOutputs, reportResults
from Shapley.visualization import showNetwork
from Shapley.networkUtilities import diamondDigger, findMasterDiamond
from Shapley.shapleyCalculation import calSV4InputFromTable
from ._core import binarizeAndPropagate, readfile, getSpeciesName, \
    getInputNames, getFormula, parseFormula, convertBooleanFormulas2Network, \
    initSpeciesStates, genInput, simulateNetworkCombined, workWithOriginalNetwork, \
    propagateCyclicVersion1, propagateAcyclic

def Analysis():
    """
    Perform Shapley value analysis using the propagation method.
    Args:
        None
    Returns:
        None
    """
    timeHeading = 0 # time to read, parse formulas 
    timegenInput = 0
    timegenFullInput = 0
    timePreprocessing = 0
    timeSimOriginal = 0
    timeSimOriginalAndInfer = 0
    timeCalc4Input = 0
    timegentable = 0
    propapationTime = 0
    timeSimulatePerturbations = 0
    timeSimulatePerFull = 0

    avewronginputkorow = 0
    avewronginputkirow = 0
    
    bound = 0
    numinputs = 0
    numnnodes = 0
    numedges = 0

    parser = parseArgumentsAnalysis()
    args= vars(parser.parse_args())
    if len(args) < 2:
            parser.print_help()

    else:
        headingTimeBegin = timeit.default_timer()
        debug = args['verbose']
        exp = args['expression']
        strout = args['output']
        samplingStrategy = args['sampling']
        samplePercentage = args['samplePer']
        simPercentage = args['simPer']
        mode = args['mode']
        scheme = args['scheme']
        sysname = args['name']
        declog = args['decimal']
        test = args['test']
        outputname = strout

        
        

        print("========== STARTING BOOLEAN SHAPLEY ANALYSIS VIA PROPAGATION FOR {0: <46}=========".format(sysname))
        print("========== EXPRESSION FILE IS {0: <81}=========".format(exp))
        print("========== THE OUTPUT OF INTEREST IS: {0: <73}=========".format(outputname))
        print("========== SAMPLING STRATEGY IS: {0: <78}=========".format(samplingStrategy))
        print("========== SCHEME IS: {0: <89}=========".format(scheme))
        print("========== MODE IS: {0: <91}=========".format(mode))

        
        lines = readfile(exp, debug)

        speciesnames = getSpeciesName(lines, debug) # get species names 
        print("List of all species:".upper())
        print(sorted(list(speciesnames)))

        if outputname not in speciesnames:
            print("The list of interested output is not a subset of species")
        
        inputnames = getInputNames(lines, speciesnames, debug)
        print("-----Input names-----".upper())
        print(sorted(list(inputnames)))

        bound = len(inputnames) + 1
        numinputs = len(inputnames)


        internames = speciesnames.difference(inputnames).difference((outputname))

        print("----Intermediate nodes-----".upper())
        print(sorted(list(internames)))

        print("===============================================================")
        sortedinputs, sortedinters = getOrderedList(inputnames, internames, True)
        print("SORTED INPUT")
        print(sortedinputs)
        print("SORTED INTER")
        print(sortedinters)
        if debug:
            print("----Getting boolean formulas-----")
        strformulas = getFormula(lines, debug) # list of formulas in strings 
        formulas = dict() 
        for strformula in strformulas:
            root = parseFormula(strformula, debug)
            if debug:
                print("Parsing formula for {}".format(strformula['left']))
                root.display()
                print("\n")
            formulas[strformula['left']] = root 
            
        orinet, orifas = convertBooleanFormulas2Network(formulas, inputnames, \
                                                        speciesnames, debug=debug)
        
        numedges = orinet.number_of_edges()
        numnnodes = orinet.number_of_nodes()
        
        numscc = nx.number_strongly_connected_components(orinet)

   
        print("THERE ARE {} FEEDBACK ARCS IN THE ORIGINAL NETWORK:".format(len(orifas)))
        print(orifas)

        

        headingTimeEnd = timeit.default_timer()
        timeHeading = headingTimeEnd - headingTimeBegin 

        showNetwork(orinet, None, None, None, None, "ori" + sysname + ".html")
        


        timegenInput = timeit.default_timer()
        species = initSpeciesStates(speciesnames, debug)
        sampledinputstates, maprowid2inputs = genInput(species, sortedinputs, samplingStrategy, samplePercentage, debug)
        timegenInput = timeit.default_timer() - timegenInput

        
        print("-----Mapping rowid to inputs:-----")
        for mrowid, minputs in maprowid2inputs.items():
            print(f"{mrowid}: {minputs}")

        if test:
            # generate full set of input configurations for testing 
            timegenFullInput = timeit.default_timer()
            fullinputstates, testmaprowid2inputs = genInput(species, inputnames, 'all', 100, debug)
            timegenFullInput = timeit.default_timer() - timegenFullInput

            print("TEST: Mapping rowid to inputs:")
            for trowid, tinputs in testmaprowid2inputs.items():
                print(f"{trowid}: {tinputs}")
        

        # binarize if mode is propagation
        if mode == 'Propagation':
            timePreprocessing = timeit.default_timer()
            print("-------------Now convert the original network to the acyclic network---------------")
            # sortedinput, sortedinter = getOrderedList(inputnames, internames, debug)
            anet, aformulas, extranodes, nodes_positions = manipulateNetwork(orinet, inputnames, formulas, True, False, debug)

            # extranodes = dict()
            # for extranode in extranodes:
            #     extranodes[extranode] = False

            print("-------------Now convert the acyclic network to the binary network---------------")
            bistrformulas = toBinaryFormulas(aformulas, debug)
            binet, nodes_positions = convertBiBooleanFormulas2Network(bistrformulas, inputnames, speciesnames, 'binet' + sysname, acyclic = False, debug=debug, extranodes=extranodes)
            biformulas = [] # list of dictionaries of formulas in exact ordering to update value 
            biformulasdict = dict() # dictionary with keys are nodes for accessibility 
            bispeciesnames = set()
            orderedBinodes = dict()
            for i, formula in enumerate(bistrformulas):
                term = formula['term'] 
                orderedBinodes[term] = i
                bistrformula = formula['formula']
                thisfor = {'left': term, "right": bistrformula}
        
                thisbiformula = dict()
                thisbiformula['term'] = term
                thisbiformula['formula'] = parseFormula(thisfor, False) 
                biformulasdict[term] = thisbiformula['formula'] 
        
                biformulas.append(thisbiformula)
        
                bispeciesnames.add(term)
                coms = bistrformula.split()
                for com in coms:
                    if com != '(' and com != ')' and com != 'AND' and com != 'OR' and com != 'NOT' and com != '=':
                        if '-' in com: # this is for excluding '-' in the formulas which is confusing for the compiler 
                            com = com.replace('-','_')
                        bispeciesnames.add(com)

            biinternames = bispeciesnames.difference(inputnames).difference((outputname))

            timePreprocessing = timeit.default_timer() - timePreprocessing
            
            showNetwork(binet, None, None, None, None, 'binet' + sysname + '.html')
            if debug:
                for formula in biformulas:
                    print("{}=".format(formula['term']))
                    formula['formula'].display()
                    print("\n")
        else:
            biformulas = None
            extranodes = None
            biformulasdict = None
            bistrformulas = None
            biinternames = internames


        print('----------Now simulate the original network----------')

        if test: 
            testinputstatesclone = copy.deepcopy(fullinputstates)
            timeSimOriginal = timeit.default_timer()
            gtoutputs, gtblinkings, gtdecimalPairs = simulateNetworkCombined(testinputstatesclone, sortedinputs, sortedinters, formulas, None, None, debug)
            gttable, index, aindex, bindex  = genTableFromOutput(gtoutputs, inputnames, gtblinkings, sortedinters, outputname, scheme, debug)
            gtkorows, gtkirows, gtkovalues, gtkivalues = calSV4InputFromTable(gttable, sortedinputs, outputname, debug)
            timeSimOriginal = timeit.default_timer() - timeSimOriginal 
            if declog:
                for tindec, toutdec in gtdecimalPairs.items():
                    print(f"{tindec}: {toutdec}") 

            


        # simulate the original network with option to infer value for binarized nodes and generate table with passed scheme 
        timeSimOriginalAndInfer = timeit.default_timer()
        sampledinputstatesclone1 = copy.deepcopy(sampledinputstates)
        outputs, blinkings, decimalPairs = simulateNetworkCombined(sampledinputstatesclone1, sortedinputs, sortedinters, formulas, extranodes, biformulas, debug)
        timeSimOriginalAndInfer = timeit.default_timer() - timeSimOriginalAndInfer
        if declog:
            for indec, outdec in decimalPairs.items():
                print(f"{indec}: {outdec}") 

        timegentable = timeit.default_timer()
        print("----------Now generate the truth table from simulation outputs----------")
        table, index, aindex, bindex = genTableFromOutput(outputs, inputnames, blinkings, biinternames, outputname, scheme, debug) 
        timegentable = timeit.default_timer() - timegentable
        
        timeCalc4Input = timeit.default_timer()
        koinrows, kiinrows, koinvalues, kiinvalues = calSV4InputFromTable(table, sortedinputs, outputname, debug)
        timeCalc4Input = timeit.default_timer() - timeCalc4Input

        if debug:
            for row, state in table.items():
                print(f"{row}: {state}")
        

        

        # if mode == "Propagation":
        #     print(f"TIME simulate original network AND INFER values of extranodes: {timeSimOriginalAndInfer}")
        # elif mode == "Simulation":
        #     print(f"TIME simulate original network: {timeSimOriginalAndInfer}")

        
        korows = dict() # dict of list of rows associated to nodes for knockout values, keys are nodes 
        kirows = dict() # dict of list of rows associated to nodes for knockin values, keys are nodes 
        kovalues = dict()
        kivalues = dict() 
        
        if mode == "Propagation":            
            print(f"-----------Now propagate within binarized network for output {outputname}-------------")
            propapationTime = timeit.default_timer()
            diamonds = diamondDigger(binet, outputname, biformulasdict, debug)
            if debug:
                print("-----Diamonds found in the binary network-----")
                for begin, ends in diamonds.items():
                    print(f"Node {begin} has diamonds: {ends}")
            masterDiamonds = findMasterDiamond(binet, diamonds)
            if debug:
                print("----Master diamonds-----")
                print(masterDiamonds)

            
            rowsofnodes = propagateCyclicVersion1(binet, table, index, aindex, \
                                                    outputname, biformulasdict, formulas, \
                                                    speciesnames, extranodes, masterDiamonds, \
                                                    inputnames, orderedBinodes, simPercentage, debug)

            for node, rows in rowsofnodes.items():
                if node in inputnames or node in internames:
                    kovalues[node], kivalues[node] = 0.0, 0.0
                    korows[node], kirows[node] = set(), set()
                    for arowid in rows:
                        onerow = table[arowid]
                        if onerow[outputname]:
                            if onerow[node]:
                                kovalues[node] += onerow['PROP']
                                korows[node].add(arowid)
                            else:
                                kivalues[node] += onerow['PROP']
                                kirows[node].add(arowid)
                        else:
                            if onerow[node]:
                                kovalues[node] -= onerow['PROP']
                                korows[node].add(arowid)
                            else:
                                kivalues[node] -= onerow['PROP']
                                kirows[node].add(arowid)

                if node in extranodes:
                    try:
                        rootname = node.split("_to_")[0]
                        if rootname not in inputnames and rootname not in internames:
                            print("There is no node named {} as the root node of {}".format(rootname, node))
                        else:
                            if rootname not in korows:
                                korows[rootname], kirows[rootname] = set(), set()
                                kovalues[rootname], kivalues[rootname] = 0.0, 0.0 
                            for arowid in rows:
                                if arowid not in korows[rootname].union(kirows[rootname]):
                                    onerow = table[arowid]
                                    if onerow[outputname]:
                                        if onerow[rootname]:
                                            kovalues[rootname] += onerow['PROP']
                                            korows[rootname].add(arowid)
                                        else:
                                            kivalues[rootname] += onerow['PROP']
                                            kirows[rootname].add(arowid)
                                    else:
                                        if onerow[rootname]:
                                            kovalues[rootname] -= onerow['PROP']
                                            korows[rootname].add(arowid)
                                        else:
                                            kivalues[rootname] -= onerow['PROP']
                                            kirows[rootname].add(arowid)
                    except:
                        print(f"ERROR: Problem during post processing of node {rootname}")

            zeronodes = (inputnames.union(internames)).difference({*rowsofnodes})      
            for node in zeronodes:
                korows[node], kirows[node] = set(), set() 
                kovalues[node], kivalues[node] = 0.0, 0.0      

            # for node, testrows in korows.items():
            #     print("TESSSSSSTTTTTTT")
            #     print("Node {}".format(node))
            #     print("KO rows: {}".format(sorted(list(testrows))))
            
            print("-----Now estimate the error by comparing propagated values to simulation values of input nodes-----")
            print("For KO values")
            avewronginputkorow = estimateErrorFromInputs(koinrows, korows, sortedinputs)

            print("For KI values")
            avewronginputkirow = estimateErrorFromInputs(kiinrows, kirows, sortedinputs)

            # korows.update(koinrows)
            # kirows.update(kiinrows)
            # kovalues.update(koinvalues)
            # kivalues.update(kiinvalues)

            propapationTime = timeit.default_timer()- propapationTime

            print("=========PROPAGATION RESULTS=========")
            reportResults(kovalues, kivalues, sorted(list(inputnames)), sorted(list(internames)), outputname, korows, kirows)

        else: # mode == 'Simulation':  # case of Simulation 
            print("-------Now performing pertubation simulations---------")
            timeSimulatePerturbations = timeit.default_timer()
            for internode in sortedinters:
                if debug:
                    print(f"Knocing OUT {internode}")
                sampledinputstatesclone2 = copy.deepcopy(sampledinputstates)
                kooutputs, _, kodecpairs = simulateNetworkCombined(sampledinputstatesclone2, sortedinputs, sortedinters, formulas, None, None, debug, [internode], None)

                if declog:
                    for koindec, kooutdec in kodecpairs.items():
                        print(f"{koindec}: {kooutdec}")


                if debug: 
                    print(f"Knocing IN  {internode}")
                sampledinputstatesclone3 = copy.deepcopy(sampledinputstates)
                kioutputs, _, kidecpairs = simulateNetworkCombined(sampledinputstatesclone3, sortedinputs, sortedinters, formulas, None, None, debug, None, [internode])

                if declog:
                    for kiindec, kioutdec in kidecpairs.items():
                        print(f"{kiindec}: {kioutdec}")

                
                korow, kirow = compareOutputs(outputs, outputname, inputnames, kooutputs, kioutputs, debug)
                korows[internode] = sorted(list(korow))
                kirows[internode] = sorted(list(kirow))
                kovalue, kivalue = 0.0, 0.0 

                for rowid in korow:
                    row = table[rowid]
                    if row[outputname]:
                        kovalue += row['PROP']
                    else:
                        kovalue -= row['PROP']
                for rowid in kirow:
                    row = table[rowid]
                    if row[outputname]:
                        kivalue += row['PROP']
                    else:
                        kivalue -= row['PROP']
                kovalues[internode], kivalues[internode] = kovalue, kivalue

            korows.update(koinrows)
            kirows.update(kiinrows)
            kovalues.update(koinvalues)
            kivalues.update(kiinvalues)
            timeSimulatePerturbations = timeit.default_timer()- timeSimulatePerturbations
            

            print("=========SIMULATION RESULTS=========")
            reportResults(kovalues, kivalues, sorted(list(inputnames)), sorted(list(internames)), outputname, korows, kirows)
            

        if test:
            print('-------In Test mode, perform knocking out with full input configurations to check accuracy')
            timeSimulatePerFull = timeit.default_timer()
            gtinterkorows = dict() # dict of list of rows associated to nodes for knockout values, keys are nodes 
            gtinterkirows = dict() # dict of list of rows associated to nodes for knockin values, keys are nodes 
            gtinterkovalues = dict()
            gtinterkivalues = dict()

            for internode in sortedinters:
                if debug:
                    print(f"Knocing OUT {internode}")
                fullinputstateclone1 = copy.deepcopy(fullinputstates)
                gtkooutputs, _, _ = simulateNetworkCombined(fullinputstateclone1, sortedinputs, sortedinters, formulas, None, None, debug, [internode], None)

                if debug: 
                    print(f"Knocing IN  {internode}")
                fullinputstateclone2 = copy.deepcopy(fullinputstates)
                gtkioutputs, _, _ = simulateNetworkCombined(fullinputstateclone2, sortedinputs, sortedinters, formulas, None, None, debug, None, [internode])

                korow, kirow = compareOutputs(gtoutputs, outputname, inputnames, gtkooutputs, gtkioutputs, debug)
                gtinterkorows[internode] = sorted(list(korow))
                gtinterkirows[internode] = sorted(list(kirow))
                kovalue, kivalue = 0.0, 0.0 

                for inputid, row in gttable.items():
                    if inputid in korow:
                        if row[outputname]:
                            kovalue += row['PROP']
                        else:
                            kovalue -= row['PROP']
                    if inputid in kirow:
                        if row[outputname]:
                            kivalue += row['PROP']
                        else:
                            kivalue -= row['PROP']
                gtinterkovalues[internode], gtinterkivalues[internode] = kovalue, kivalue

            
            gtkorows.update(gtinterkorows)
            gtkirows.update(gtinterkirows)
            gtkovalues.update(gtinterkovalues)
            gtkivalues.update(gtinterkivalues) 
            timeSimulatePerFull = timeit.default_timer() - timeSimulatePerFull

            print("-----Report result for TEST-----")
            reportResults(gtkovalues, gtkivalues, sorted(list(inputnames)), sorted(list(internames)), outputname, gtkorows, gtkirows)

            
            print("-----Now compare the calculated rows with grouth truth-----")
            # for nodes in sorted(list) 
            evaluation2(korows, kirows, kovalues, kivalues, gtkorows, gtkirows, gtkovalues, gtkivalues, sortedinputs, outputname, debug)
        

    print("===============MODEL information===============")
    print(f"NUMBER OF INPUT: \t {numinputs}")
    print(f"NUMBER OF NODES: \t {numnnodes}")
    print(f"NUMBER OF EDGES: \t {numedges}")
    print(f"NUMBER OF SCCs: \t {numscc}")
    print(f"NUMBER OF aFAS: \t {len(orifas)}")
    print("===============TIME information===============")
    print(f"TIME process heading: {timeHeading}")
    print(f"TIME generate input states: {timegenInput}")
    print(f"TIME convert to BBN (Propagation mode): {timePreprocessing}")
    print(f"TIME to fully simulate original network (if Test): {timeSimOriginal}")
    print(f"TIME simulate network (with sampled method) AND INFER values of extranodes if Propagation: {timeSimOriginalAndInfer}")
    print(f"TIME generating table {timegentable}")
    print(f"TIME calculating Shapley values for input nodes: {timeCalc4Input}")
    print(f"TIME propagate (in Propagation mode): {propapationTime}") 
    print(f"TIME simulate perturbations (with sampling percentage) in Simulation mode: {timeSimulatePerturbations}") 
    print(f"TIME simulate perturbations full if Test: {timeSimulatePerFull}") 

    if mode == 'Propagation':
        print(f"TIME total for Propagation: {timeHeading + timegenInput + timePreprocessing + timeSimOriginalAndInfer + timegentable + timeCalc4Input + propapationTime}")
    if mode == 'Simulation':
        print(f"TIME total for Simulation: {timeHeading + timegenInput + timeSimOriginalAndInfer + timegentable + timeCalc4Input + timeSimulatePerturbations}")
    if test:
        print(f"TIME total for FULLY Simulation: {timeHeading + timegenFullInput + timeSimOriginal + timeSimulatePerFull}")

    print("==========ERROR ESTIMATION FROM INPUTS=========")
    print(f"Average error for KO values: {avewronginputkorow}")
    print(f"Average error for KI values: {avewronginputkirow}")





        
# def runBatchTest():
#     """
#     Run batch test for multiple expression files in a folder.
#     Args:
#         None

def OldVersionTest():
    """
    Test function to compare perfomance of standard Shapley value calculation and the propagation method.
    Args:
        None
    Returns:
        None
    """
    timeSimOriginalNetwork = 0 # only the time to simulate original network, lower bound 
    timePerformOriginalAnalysis = 0 # time to perform the Knockout/Knockin of the original network 
    timePropagation = 0 # time to propagate 
    timeHeading = 0 # time to read, parse formulas 
    
    bound = 0
    numinputs = 0
    numnnodes = 0
    numedges = 0

    
    parser = parseArgumentsAnalysis()
    args= vars(parser.parse_args())
    if len(args) < 2:
            parser.print_help()

    else:
        headingTimeBegin = timeit.default_timer()
        debug = args['verbose']
        exp = args['expression']
        strout = args['output']
        samplingStrategy = args['sampling']
        samplePercentage = args['samplePer']
        simPercentage = args['simPer']
        mode = args['mode']
        scheme = args['scheme']
        sysname = args['name']
        declog = args['decimal']
        outputname = strout

        print("==========STARTING BOOLEAN SHAPLEY ANALYSIS TEST==========")
        print(f"==========EXPRESSION FILE IS {exp}")
        print(f"==========THE OUTPUT OF INTEREST IS: {outputname}")
        print("===============================================================")
        
        lines = readfile(exp, debug)

        speciesnames = getSpeciesName(lines, debug) # get species names 
        print("List of all species:")
        print(sorted(list(speciesnames)))
        numnnodes = len(speciesnames)

        if outputname not in speciesnames:
            print("The list of interested output is not a subset of species")
        
        inputnames = getInputNames(lines, speciesnames, debug)
        print("-----Input names-----")
        print(sorted(list(inputnames)))
        bound = len(inputnames) + 1
        numinputs = len(inputnames)


        internames = speciesnames.difference(inputnames).difference((outputname))
        print("----Intermediate nodes-----")
        print(sorted(list(internames)))


        if debug:
            print("----Getting boolean formulas-----")
        strformulas = getFormula(lines, debug) # list of formulas in strings 
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
            print("THERE ARE {} FEEDBACK ARCS IN THE ORIGINAL NETWORK:".format(len(orifas)))
            print(orifas)

        
        # get the input states of the original network
        species = initSpeciesStates(speciesnames, debug)
        inputstates, maprowid2inputs = genInput(species, inputnames)

        for rowid, state in maprowid2inputs.items():
            print(f"{rowid}: {state}")
        

        headingTimeEnd = timeit.default_timer()
        timeHeading = headingTimeEnd - headingTimeBegin 

        simrowsofnodes, oridecimalpairs, timesimori, timeori = workWithOriginalNetwork(orinet, inputnames, \
            inputstates, outputname, internames, formulas, True, mode=mode, debug=False) 
        if declog:
            for indec, outdec in oridecimalpairs.items():
                print(f"{indec}:  {outdec}") 


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
        timeconver2acyclic = timeit.default_timer()
        anet, aformulas, extranodes, nodes_positions = manipulateNetwork(orinet, inputnames, formulas, True, False, debug)
        timeconver2acyclic = timeit.default_timer() - timeconver2acyclic

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
        speciesnames, "abinetwork", sortedinput, sortedinter, extranodes, mode, 0, totest = False, debug=debug)
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

    averightrows, ave_dncg = evaluation1(prorowsofnodes, simrowsofnodes, table, outputname)
    print(f"Average number of correct rows after binarization: {averightrows} ")
    print(f"Average DNCG after binarization: {ave_dncg} ")
    

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



if __name__ == "__main__":
    start = timeit.default_timer()
    Analysis()
    # OldVersionTest()
    # runBatchTest('realbenchmark.txt', 'cycle1_1_nolog.json')
    # testEachModel('data/realbenchmark/1.txt', 'Apoptosis')
    print("-----------TOTAL RUNNING TIME---------------")
    print(timeit.default_timer() - start)  