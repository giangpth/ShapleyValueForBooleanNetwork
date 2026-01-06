import math
from Shapley.exceptions import UnknownError
from Shapley.utilities import subsets

def extractPhe(inputnames, outputname, outputstates, debug=False):
    """
    Extract phenotype information from output states.
    Args:
        inputnames (set): A set of input species names.
        outputnames (set): A set of output species names.
        outputstates (list): A list of dictionaries, each representing an output state with species names as keys and their boolean states as values.
        debug (bool): If True, print debug information.
    Returns:
        dict: A dictionary mapping frozensets of input names to sets of output names.
    """
    genphe = dict() 
    statistic = dict()
    statisticgen = dict()

    # print("Pure input")
    for outputstate in outputstates:
        gentem = set() # will be used as key
        for name in inputnames:
            if outputstate[name]:
                gentem.add(name)
        gen = frozenset(gentem)

        phe = set() # will be used as value 
        # for name in outputnames:
        #     # print(name, outputstate[name])
        if outputstate[outputname]:
            phe.add(outputname)
        # print('')
        genphe[gen] = phe

        #from now on until the end of for loop is just for statistic purpose 
        frozenphe = frozenset(phe)
        if frozenphe not in statistic:
            statistic[frozenphe] = 1
        else:
            statistic[frozenphe] += 1

        frozengen = frozenset(gen)
        if frozengen not in statisticgen:
            statisticgen[frozengen] = 1
        else:
            statisticgen[frozengen] += 1
    return genphe


def calKSV4Input(genphe, inputnames, outputname, knockin=False, mode = 'Shapley', simtable = None,  debug=False):
    """
    Calculate Shapley values for input species based on their contributions to output species.
    Args:
        genphe (dict): A dictionary mapping frozensets of input names to sets of output
        inputnames (set): A set of input species names.
        outputnames (set): A set of output species names.
        knockin (bool): If True, perform knock-in analysis; if False, perform knock-out analysis.
        mode (str): Calculation mode: 'Shapley' or 'random'
        simtable (dict): A dictionary mapping row IDs to dictionaries of species states.
        debug (bool): If True, print debug information.
    Returns:
        dict: A dictionary mapping output names to dictionaries of input names and their Shapley values.
        dict: A dictionary mapping input names to sets of row IDs where they contributed.
    """
    countedrows = dict()
    ids = list([i for i in range(len(inputnames))])
    # print(ids)
    allsets = subsets(ids)

    # calculate for each output component
    shapss = dict()
    shaps = dict()  
    # for each input component
    for inname in inputnames: # the knockout one 
        countedrows[inname] = set()
        map = dict() 
        for i in range(len(inputnames)):
            map[i] = list(inputnames)[i]
        numrow = 0
        sum = 0
        for oneset in allsets:
            phenotem = set() 
            for e in oneset:
                phenotem.add(map[e])
            pheintact = frozenset(phenotem)
            if knockin:
                pheknockout = phenotem.union({inname})
                pheknockout = frozenset(pheknockout)
            else: 
                pheknockout = frozenset(phenotem - {inname})
            # print(pheintact, pheknockout)
            try: 
                v_intact = 0
                v_knockout = 0
                if outputname in genphe[pheintact]:
                    v_intact = 1
                if outputname in genphe[pheknockout]:
                    v_knockout = 1
                
                gain  = v_intact - v_knockout
                if mode == 'random':
                    weightedGain = gain * math.factorial(len(oneset)) * math.factorial(len(inputnames) - len(oneset))
                else:
                    weightedGain = gain * math.factorial(len(oneset)) * math.factorial(len(inputnames) - len(oneset))

                sum += weightedGain
                if weightedGain != 0:
                    numrow += 1
                    if simtable:
                        for id, line in simtable.items():
                            ok = True
                            for input in inputnames:
                                if input in pheintact:
                                    if line[input] == False:
                                        ok = False
                                        break
                                else:
                                    if line[input] == True:
                                        ok = False
                                        break
                            if ok:
                                countedrows[inname].add(id)

            except UnknownError:
                print("Error within calKSV4Input function")
                continue
        shap = round(sum/(math.factorial(len(inputnames))),4)
        shaps[inname] = shap
    if simtable:
        return shaps, countedrows
    else:
        return shaps

def calKSV(intactgenphe, knockoutgenphes, outputname, numinput, mode = 'Shapley', simtable=None, inputnames=None, debug=False):
    """
    Calculate Shapley values for intermediate species based on their contributions to output species.
    Args:
        intactgenphe (dict): A dictionary mapping frozensets of input names to sets of output names for the intact network.
        knockoutgenphes (dict): A dictionary mapping intermediate species names to dictionaries that map frozensets of input names to sets of output names for the network with that species knocked out.
        outputnames (set): A set of output species names.
        numinput (int): The number of input species.
        simtable (dict): A dictionary mapping row IDs to dictionaries of species states.
        inputnames (set): A set of input species names.
        debug (bool): If True, print debug information.
    Returns:
        dict: A dictionary mapping output names to dictionaries of intermediate species names and their Shapley values.
        dict: A dictionary mapping intermediate species names to sets of row IDs where they contributed.
    """
    countedrow = dict()

    shaps = {}
    # for each intermediate node 
    for internode, genphes in knockoutgenphes.items():
        countedrow[internode] = set()
        # print("Working with {}".format(internode))
        sum = 0
        numrow = 0
        for gen, phe in genphes.items():
            # rowid += 1
            # print(gen, phe)
            s = len(gen)
            gain = 0
            # try:
            intactphe = intactgenphe[gen]
            # print("Intact phe")
            # print(intactphe)
            pheyes = 0
            pheno = 0
            # try:
            if outputname in intactphe:
                pheyes = 1
            if outputname in phe:
                pheno = 1
            gain = pheyes - pheno
            if mode == 'Uniform':
                weightedgain = gain * math.factorial(s) * math.factorial(numinput - s)
            else:
                weightedgain = gain * math.factorial(s) * math.factorial(numinput - s)
            sum += weightedgain
            if weightedgain != 0 and simtable:
                numrow += 1
                for id, line in simtable.items():
                    ok = True
                    for input in inputnames:
                        if input in gen:
                            if not line[input]:
                                ok = False
                                break
                        if input not in gen:
                            if line[input]:
                                ok = False
                                break
                    if ok:
                        countedrow[internode].add(id)
        if debug:
            print("Number of rows for {} is {}".format(internode, numrow))
        shaps[internode] = round(sum/math.factorial(numinput),4)
    if simtable:
        return shaps, countedrow
    else:
        return shaps, None
