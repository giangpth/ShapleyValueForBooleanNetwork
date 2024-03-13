import argparse as ap
import itertools
from typing import Dict, Any
import hashlib
import json
import math
import copy 
from booleanFormulaParser import parseFormula, deleteSomeNodes, getResult


def dict_hash(dictionary: Dict[str, Any]) -> str:
    """MD5 hash of a dictionary."""
    dhash = hashlib.md5()
    # We need to sort arguments so {'a': 1, 'b': 2} is
    # the same as {'b': 2, 'a': 1}
    encoded = json.dumps(dictionary, sort_keys=True).encode()
    dhash.update(encoded)
    return dhash.hexdigest()


def parseArguments():
    parser = ap.ArgumentParser(description='Input the expression file and ')
    parser.add_argument('-e', '--expression', type=str, help="path to expression file", required=True)
    parser.add_argument('-o', '--output', type=str, help='list of interested components to examine')
    parser.add_argument('-d', '--debug', help="Print log to debug or not", \
                        action='store_true')
    return parser

def readfile(path, debug=False):
    try:
        f = open(path, 'r')
        lines = []
        for line in f:
            line = line.strip()
            if debug:
                print(line)
            lines.append(line)
        return lines


    except IOError:
        print("Error opening file")
        return -1
    
def getSpeciesName(lines, debug=False):
    species = set()
    for line in lines:
        coms = line.split()
        for com in coms:
            if com != '(' and com != ')' and com != 'AND' and com != 'OR' and com != 'NOT' and com != '=':
                species.add(com)
    if debug:
        print("----GET SPECIES----")
        print(species)
    return species

def genSpecies(speciesnames, debug=False):
    species = dict()
    for s in speciesnames:
        species[s] = None
    return species

def getInputNames(lines, speciesnames, debug=False):
    # input = set()
    dep = set() 
    for line in lines:
        coms = line.split('=') 
        # print(coms)
        for tem in coms[0].split():
            dep.add(tem)
    # print(len(dep))
    input = speciesnames - dep
    if debug:
        print("----GET INPUT----")
        print(dep)
        print(input)
    return input

def getFormula (lines, debug=False):
    allfor = [] 
    for line in lines:
        formula = dict()
        line = line.replace(" AND ", " and ")
        line = line.replace(" OR ", " or ")
        line = line.replace(" NOT ", " not ")
        sides = line.split('=')
        assert len(sides) == 2, print("There is none or more than 1 equal sign in the formula")
        formula['left'] = sides[0].strip()
        formula['right'] = sides[1].strip()
        allfor.append(formula)
    return allfor 

#species is a dictionary, inputnames is a list 
def genInput(species, inputnames, debug=False):
    numin = len(inputnames)
    coms = list(itertools.product([0, 1], repeat=numin))
    if debug:
        print("There are {} combinations for input".format(len(coms)))
        # print (coms)
    inputstates = []
    for com in coms:
        inp = copy.deepcopy(species)
        # print(inp)
        for id, name in enumerate(inputnames):
            try:
                # print(name, id)
                if com[id] == 1:
                    inp[name] = True
                elif com[id] == 0:
                    inp[name] = False
                else:
                    inp[name] = com[id]
            except KeyError:
                print("There is no key named {}".format(name))
        inputstates.append(inp)
    
    return inputstates

    # from now on to the end of the function is just to test 
    # teminp = copy.deepcopy(species)
    # for id, name in enumerate(inputnames):

    #     teminp[name] = True
    #         # try:
    #             # print(name, id)
    #         #     if com[id] == 1:
    #         #         inp[name] = True
    #         #     elif com[id] == 0:
    #         #         inp[name] = False
    #         #     else:
    #         #         inp[name] = com[id]
    #         # except KeyError:
    #         #     print("There is no key named {}".format(name))
    # return [teminp]


def sim1step(formulas, inputstate, debug=False):
    # run command to declare all the variables 
    for spe in inputstate.keys():
        # print(spe)
        command =  str(spe) + ' = ' + str(inputstate[spe])
        # print(command)
        exec(command)
    

    # run all the formulas with a copy of the variables 
    for formula in formulas:
        # formula['right'].display()
        # command = formula['left'] + '_tem_ = ' + formula['right'] 
        res = getResult(formula['right'], inputstate, debug)
        command = formula['left'] + '_tem_ = ' + str(res)
        exec(command)

    # STAT6_tem_ = None

    # assign new value to real variables taking from value of copy variables
    for formula in formulas:
        command = formula['left'] + ' = ' + formula['left'] + '_tem_'
        exec(command)
        # print(command)
    
    # put back to inputstate to return 
    for spe in list(inputstate):
        command = 'inputstate[\'' + str(spe) + '\'] = ' + str(spe)
        # print(command)
        exec(command)
    # exec("print(STAT6)")
    return inputstate

def merge2states(state1, state2, debug=False):
    for item in state1.keys():
        # print(state1[item], state2[item])
        state1[item] = state1[item] or state2[item]
    return state1

# inputstate is a list of all species with the predefine value for input species
# this function is the simulation process 
def getOutput(formulas, inputstate, debug=False, maxStep = 10000) -> dict: 
    oldstate = dict()
    numstep = 0
    for _ in range(maxStep):
        numstep += 1
        inputstate = sim1step(formulas, inputstate, debug)
        hash = dict_hash(inputstate)


        if hash not in oldstate:
            oldstate[hash] = numstep
        else:
            if debug:
                print("Number of iteration {} and first loop point is {}".format(numstep, oldstate[hash]))

            # merge all the state inside the loop 
            returnstate = copy.deepcopy(inputstate)
            for i in range(numstep-oldstate[hash] - 1):
                inputstate = sim1step(formulas, inputstate, debug)
                returnstate = merge2states(returnstate, inputstate)
            return returnstate 
    # print(type(inputstate))
    return inputstate 

def subsets(s):  
    if len(s) == 0:  
        return [[]]  
    x = subsets(s[:-1])  
    return x + [[s[-1]] + y for y in x] 

def extractPhe(inputnames, outputnames, outputstates, inputstates=None, debug=False):
    genphe = dict() 
    statistic = dict()
    statisticgen = dict()
    if not inputstates: # take inputstates from outputstates in case of pure input, which will not change over the simulation 
        print("Pure input")
        for outputstate in outputstates:
            gentem = set() # will be used as key
            for name in inputnames:
                if outputstate[name]:
                    gentem.add(name)
            gen = frozenset(gentem)

            phe = set() # will be used as value 
            for name in outputnames:
                # print(name, outputstate[name])
                if outputstate[name]:
                    phe.add(name)
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

        print (statistic)
        print(statisticgen)
    else:
        print("Impure input")
        assert len(inputstates) == len(outputstates), "Different length of input and output states"
        for id, instate in enumerate (inputstates):
            outstate = outputstates[id]

            gentem = set ()
            for name in inputnames:
                if instate[name]:
                    gentem.add(name)
            gen = frozenset(gentem)

            phe = set()
            for name in outputnames:
                if outstate[name]:
                    phe.add(name)
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

        print (statistic)
        print(statisticgen)

    return genphe


# calculate original shapley value based on simulation result 
def calOriShapleyInput(genphe, inputnames, outputnames, debug=False):
    ids = list([i for i in range(len(inputnames) - 1)])
    # print(ids)
    allsets = subsets(ids)

    # calculate for each output component
    shapss = dict()
    for outname in outputnames:
        shaps = dict()  
        # for each input component
        for inname in inputnames:
            rest = list(inputnames - {inname})
            # print(rest)
            map = dict() 
            for i in range(len(inputnames) - 1):
                map[i] = rest[i]

            sum = 0
            for oneset in allsets:
                phenotem = set() # subset without element i 
                for e in oneset:
                    phenotem.add(map[e])
                # print(phe)
                pheno = frozenset(phenotem)
                
                # then add also the considered element i 
                phenotem.add(inname) # subset with element i 

                pheyes = frozenset(phenotem)
                
                # calculate gain, gain = v(pheon) - v(pheyes)
                # v(S) = 1 if outname is presence, v(S) = 0 if outname is not presence 
                # v_pheno = genphe[pheno] 
                try:
                    if outname in genphe[pheno]: 
                        v_pheno = 1
                    else:
                        v_pheno = 0

                    if outname in genphe[pheyes]:
                        v_pheyes = 1
                    else:
                        v_pheyes = 0

                    gain = v_pheyes - v_pheno
                    weightgain = gain * math.factorial(len(oneset)) * math.factorial(len(inputnames) - len(oneset) - 1)
                    sum += weightgain
                except:
                    continue
            
            shap = sum/(math.factorial(len(inputnames))) # this one is when divided by N! 
            # print(shap)
            shaps[inname] = round(shap, 3)
        shapss[outname] = shaps
    return shapss
        

def propagate1step(formulas, state, debug=False):
    for spe in state.keys():
        command = str(spe) + ' = ' + str(state[spe])
        # print(command)
        exec(command)
        
    for formula in formulas:
        formula['right'] = formula['right'].replace(' or ', ' + ')
        formula['right'] = formula['right'].replace(' OR ', ' + ')
        formula['right'] = formula['right'].replace(' and ', ' + ')
        formula['right'] = formula['right'].replace(' AND ', ' + ')
        formula['right'] = formula['right'].replace(' NOT ', ' - ')
        formula['right'] = formula['right'].replace(' not ', ' - ')

        command = formula['left'] + '_tem_ = ' + formula['right']
        # print(command)
        exec(command)
        
    for formula in formulas:
        command = formula['left'] + ' = ' + formula['left'] + '_tem_'
        # print(command)
        exec(command)
        
    for spe in list(state):
        command = 'state[\'' + str(spe) + '\'] = ' + str(spe)
        # print(command)
        exec(command)
    return state


def propagateSV(SVs, species, formulas, debug=False, maxStep=1000):
    ress = []
    for sv in SVs:
        res = {}
        # print(sv)
        inputs = SVs[sv]
        # print(input)
        speciescopy = copy.deepcopy(species)
        for spe in speciescopy:
            if spe in inputs:
                speciescopy[spe] = SVs[sv][spe]
            else:
                speciescopy[spe] = 0
        # for input in inputs:
        #     # print(input)
        #     speciescopy[input] = SVs[sv][input]
        print(speciescopy)
        oldstate = dict()
        numstep = 0
        for _ in range(maxStep):
            numstep += 1
            speciescopy = propagate1step(formulas, speciescopy, debug)
            print(speciescopy)
            hash = dict_hash(speciescopy)
            if hash not in oldstate:
                oldstate[hash] = numstep
            else:
                if debug:
                    print("Number of iteration {} and first loop point is {}".format(numstep, oldstate[hash]))
                returnstate = copy.deepcopy(speciescopy)
                for i in range(numstep-oldstate[hash] - 1):
                    print(i)
                ress.append(returnstate)
                break
        

def main_one_step():
    parser = parseArguments()
    args= vars(parser.parse_args())
    if len(args) < 2:
        parser.print_help()
    else:
        debug = args['debug']
        exp = args['expression']
        strouts = args['output'].split()
        outputnames = set()
        for strout in strouts:
            outputnames.add(strout)
        if debug:
            print("The expression file is {}".format(exp))
            print("The interested output is: ")
            print(outputnames)

        lines = readfile(exp, debug)
        speciesnames = getSpeciesName(lines, debug) # set of species names 
        if not outputnames.issubset(speciesnames):
            print("The list of interested output is not a subset of species")
            return -1
        species = genSpecies(speciesnames, debug) # dictionary of species with state None 
        print("---SPECIES---")
        print(species)
        inputnames = getInputNames(lines, speciesnames, debug) # set of names of input nodes 
        inputstates = genInput(species, inputnames, debug) # list of all possible input combinations 
        
        # print("----Generated inputstates----")
        # print(inputstates)

        print("-----Input names-----")
        print(inputnames)

        internames = speciesnames.difference(inputnames).difference(outputnames)
        print("----Intermediate nodes-----")
        print(internames)

        print("----Getting boolean formulas-----")
        strformulas = getFormula(lines, debug)
        formulas = []
        for formula in strformulas:
            root = parseFormula(formula, debug)
            root.display() 
            thisfor = {'left': formula['left'], 'right': root}
            formulas.append(thisfor)
        
        print("----Simulate the network with all possible combinations of inputs-----")
        # print(all)
        outputs = [] # list of output states 
        
        inputstatescopy = copy.deepcopy(inputstates)

        for inputstate in inputstatescopy:
            # print (inputstate)
            output = getOutput(formulas, inputstate, debug, 1000)
            outputs.append(output)
            # print(output)
            # print('\n')

        # print("-----Output------")
        # for output in outputs:
        #     print(output)

        genphe = extractPhe(inputnames, outputnames, outputs, None)
        interphe = extractPhe(internames, outputnames, outputs)
        shaps = calOriShapleyInput(genphe, inputnames, outputnames)
        intershaps = calOriShapleyInput(interphe, internames, outputnames)
        print("-----Shape for input with input simulation----")
        print("----Input shapley value-----")
        for item in shaps.items():
            print(item)     

        print("----Intermediate shapley value-----")
        for item in intershaps.items():
            print(item) 
        
        print("------------------------------------------------------------------------------------------------------")
        print("\n")
        print("----- Now simulate with intermediate nodes as input-----")
        intercoms = genInput(species, internames, debug) 
        
        print("----New boolean formulas excluding input nodes-----")
        print("After excluding {}".format(inputnames))
        newformulas = []
        for formula in strformulas:
            root = parseFormula(formula, debug)
            deleteSomeNodes(root, inputnames)
            # print("After excluding {}".format(inputnames))
            if root and root.val: # dont add None formula
                thisfor = {'left': formula['left'], 'right': root}
                newformulas.append(thisfor)
                print("Added")
                root.display()


        interoutputs = [] # list of output states 
        intercomcopy1 = copy.deepcopy(intercoms)
        for intercom in intercomcopy1:
            # print (intercom)
            interoutput = getOutput(newformulas, intercom, debug, 1000)
            # print(interoutput)
            interoutputs.append(interoutput)
            # print('\n')

        # one step simulation
        onestepoutput = []
        intercomcopy2 = copy.deepcopy(intercoms)
        for intercom in intercomcopy2:
            out = sim1step(newformulas, intercom)
            onestepoutput.append(out)
        
        onlyinterphe = extractPhe(internames, outputnames, interoutputs, onestepoutput)
        # print("----Only interphe----")
        # for tem in onlyinterphe.items():
        #     print(tem)

        onlyintershapss = calOriShapleyInput(onlyinterphe, internames, outputnames)
        print("----Only intermediate shapley value-----")
        for item in onlyintershapss.items():
            print(item) 
        

        
def main_one_step_2():
    parser = parseArguments()
    args= vars(parser.parse_args())
    if len(args) < 2:
        parser.print_help()
    else:
        debug = args['debug']
        exp = args['expression']
        strouts = args['output'].split()
        outputnames = set()
        for strout in strouts:
            outputnames.add(strout)
        if debug:
            print("The expression file is {}".format(exp))
            print("The interested output is: ")
            print(outputnames)

        lines = readfile(exp, debug)
        speciesnames = getSpeciesName(lines, debug) # set of species names 
        if not outputnames.issubset(speciesnames):
            print("The list of interested output is not a subset of species")
            return -1
        species = genSpecies(speciesnames, debug) # dictionary of species with state None 
        print("---SPECIES---")
        print(species)
        inputnames = getInputNames(lines, speciesnames, debug) # set of names of input nodes 
        inputstates = genInput(species, inputnames, debug) # list of all possible input combinations 
        
        # print("----Generated inputstates----")
        # print(inputstates)

        print("-----Input names-----")
        print(inputnames)

        internames = speciesnames.difference(inputnames).difference(outputnames)
        print("----Intermediate nodes-----")
        print(internames)

        print("----Getting boolean formulas-----")
        strformulas = getFormula(lines, debug)
        formulas = []
        for formula in strformulas:
            root = parseFormula(formula, debug)
            root.display() 
            thisfor = {'left': formula['left'], 'right': root}
            formulas.append(thisfor)
        
        print("----Simulate the network with all possible combinations of inputs-----")
        # print(all)
        outputs = [] # list of output states 
        
        inputstatescopy = copy.deepcopy(inputstates)

        for inputstate in inputstatescopy:
            # print (inputstate)
            output = getOutput(formulas, inputstate, debug, 1000)
            outputs.append(output)
            # print(output)
            # print('\n')

        # print("-----Output------")
        # for output in outputs:
        #     print(output)
        inputstatescopy2 = copy.deepcopy(inputstates)
        outsonestep = []

        for inputstate in inputstatescopy2:
            out = sim1step(formulas, inputstate)
            outsonestep.append(out)
        

        genphe = extractPhe(inputnames, outputnames, outputs, None)
        interphe = extractPhe(internames, outputnames, outputs, outsonestep)
        shaps = calOriShapleyInput(genphe, inputnames, outputnames)
        intershaps = calOriShapleyInput(interphe, internames, outputnames)
        print("-----Shape for input with input simulation----")
        print("----Input shapley value-----")
        for item in shaps.items():
            print(item)     

        print("----Intermediate shapley value-----")
        for item in intershaps.items():
            print(item) 
        
        print("------------------------------------------------------------------------------------------------------")
        print("\n")
        print("----- Now simulate with intermediate nodes as input-----")
        intercoms = genInput(species, internames, debug) 
        
        print("----New boolean formulas excluding input nodes-----")
        print("After excluding {}".format(inputnames))
        newformulas = []
        for formula in strformulas:
            root = parseFormula(formula, debug)
            deleteSomeNodes(root, inputnames)
            # print("After excluding {}".format(inputnames))
            if root and root.val: # dont add None formula
                thisfor = {'left': formula['left'], 'right': root}
                newformulas.append(thisfor)
                print("Added")
                root.display()


        interoutputs = [] # list of output states 
        intercomcopy1 = copy.deepcopy(intercoms)
        for intercom in intercomcopy1:
            # print (intercom)
            interoutput = getOutput(newformulas, intercom, debug, 1000)
            # print(interoutput)
            interoutputs.append(interoutput)
            # print('\n')

        # one step simulation
        onestepoutput = []
        intercomcopy2 = copy.deepcopy(intercoms)
        for intercom in intercomcopy2:
            out = sim1step(newformulas, intercom)
            onestepoutput.append(out)
        
        onlyinterphe = extractPhe(internames, outputnames, interoutputs, onestepoutput)
        # print("----Only interphe----")
        # for tem in onlyinterphe.items():
        #     print(tem)

        onlyintershapss = calOriShapleyInput(onlyinterphe, internames, outputnames)
        print("----Only intermediate shapley value-----")
        for item in onlyintershapss.items():
            print(item) 
        



def main_first_attempt():
    parser = parseArguments()
    args= vars(parser.parse_args())
    if len(args) < 2:
        parser.print_help()
    else:
        debug = args['debug']
        exp = args['expression']
        strouts = args['output'].split()
        outputnames = set()
        for strout in strouts:
            outputnames.add(strout)
        if debug:
            print("The expression file is {}".format(exp))
            print("The interested output is: ")
            print(outputnames)

        lines = readfile(exp, debug)
        speciesnames = getSpeciesName(lines, debug) # set of species names 
        if not outputnames.issubset(speciesnames):
            print("The list of interested output is not a subset of species")
            return -1
        species = genSpecies(speciesnames, debug) # dictionary of species with state None 
        print("---SPECIES---")
        print(species)
        inputnames = getInputNames(lines, speciesnames, debug) # set of names of input nodes 
        inputstates = genInput(species, inputnames, debug) # list of all possible input combinations 

        print("-----Input names-----")
        print(inputnames)

        internames = speciesnames.difference(inputnames).difference(outputnames)
        print("----Intermediate nodes-----")
        print(internames)

        # print("-----Input states-----")
        # print(inputstates)

        print("----Getting boolean formulas-----")
        strformulas = getFormula(lines, debug)
        formulas = []
        for formula in strformulas:
            root = parseFormula(formula, debug)
            root.display() 
            thisfor = {'left': formula['left'], 'right': root}
            formulas.append(thisfor)
            
    
        
        print("----Simulate the network with all possible combinations of inputs-----")
        # print(all)
        outputs = [] # list of output states 
        for inputstate in inputstates:
            # print (inputstate)
            output = getOutput(formulas, inputstate, debug, 1000)
            outputs.append(output)
            # print(output)
            # print('\n')

        # print("-----Output------")
        # for output in outputs:
        #     print(output)

        genphe = extractPhe(inputnames, outputnames, outputs)
        interphe = extractPhe(internames, outputnames, outputs)

        # print("----Genphe----")
        # for tem in genphe.items():
        #     print(tem)

        # print("----Interphe----")
        # for tem in interphe.items():
        #     print(tem)

        shapss = calOriShapleyInput(genphe, inputnames, outputnames)
        intershapss = calOriShapleyInput(interphe, internames, outputnames)

        print("----Input shapley value-----")
        for item in shapss.items():
            print(item)     

        print("----Intermediate shapley value-----")
        for item in intershapss.items():
            print(item) 

        print("\n\n")
        print("----Simulate the network with all possible combinations of intermediate nodes----")
        intercoms = genInput(species, internames, debug)
        # print(intercom)
        
        print("----New boolean formulas excluding input nodes-----")
        print("After excluding {}".format(inputnames))
        newformulas = []
        for formula in strformulas:
            root = parseFormula(formula, debug)
            deleteSomeNodes(root, inputnames)
            # print("After excluding {}".format(inputnames))
            if root and root.val: # dont add None formula
                thisfor = {'left': formula['left'], 'right': root}
                newformulas.append(thisfor)
                print("Added")
                root.display()
            print("\n\n\n")

        interoutputs = [] # list of output states 
        for intercom in intercoms:
            # print (intercom)
            interoutput = getOutput(newformulas, intercom, debug, 1000)
            # print(interoutput)
            interoutputs.append(interoutput)
            # print('\n')

        # for tem in interoutputs:
        #     print(interoutputs)
        
        onlyinterphe = extractPhe(internames, outputnames, interoutputs)
        print("----Only interphe----")
        for tem in onlyinterphe.items():
            print(tem)

        onlyintershapss = calOriShapleyInput(onlyinterphe, internames, outputnames)
        print("----Only intermediate shapley value-----")
        for item in onlyintershapss.items():
            print(item) 

        

            
if __name__ == "__main__":
    main_one_step_2()

