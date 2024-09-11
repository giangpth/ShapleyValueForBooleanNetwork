import argparse as ap
import itertools
from typing import Dict, Any
import hashlib
import json
import math
import copy 
from booleanFormulaParser import parseFormula, getResult, toBinaryFormulas, convertBiBooleanFormulas2Network, convertBooleanFormulas2Network, limitGraphAfterNode
import networkx as nx
from pyvis.network import Network
from datetime import datetime
import timeit


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
    parser.add_argument('-k', '--knockout', help='Perform knockout or not', \
                        action='store_true')
    parser.add_argument('-b', '--binary', help='Work with binary network or not', \
                        action='store_true')
    parser.add_argument('-a', '--acyclic', help='Extract acyclic network with respect to output node', \
                        action='store_true')
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
                if '-' in com: # this is for excluding '-' in the formulas which is confusing for the compiler 
                    com = com.replace('-','_')
                species.add(com)
    if debug:
        print("----GET SPECIES----")
        print(species)
    return species

def genSpecies(speciesnames, debug=False):
    species = dict()
    for s in speciesnames:
        species[s] = False
    return species

def getInputNames(lines, speciesnames, debug=False):
    # input = set()
    dep = set() 
    for line in lines:
        coms = line.split('=') 
        # print(coms)
        for tem in coms[0].split():
            if '-' in tem: # this is for excluding '-' in the formulas which is confusing for the compiler 
                tem = tem.replace('-','_')
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

        # this is for excluding '-' in the formulas which is confusing for the compiler 
        left = sides[0].strip()
        right = sides[1].strip()

        if '-' in left:  # this is for excluding '-' in the formulas which is confusing for the compiler 
            left = left.replace('-', '_')

        if '-' in right:  # this is for excluding '-' in the formulas which is confusing for the compiler 
            right = right.replace('-', '_')

        formula['left'] = left
        formula['right'] = right
        # formula[left] = right 
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

# test the correctness of the function toBinaryNetwork      
def simBinaryNetwork(biformulas, inputnames, speciesnames, debug=False, maxstep=1000):
    print("----Simulate binary network----")
    species = genSpecies(speciesnames, debug) # take list the names of all the species and set them to be None 
    inputstates = genInput(species, inputnames, debug) # each inputstate contains all species, input species are set to a certain value 

    outputs = []
    for inputstate in inputstates:
        output = getOutput(biformulas, inputstate, True, 1000, debug) 
        outputs.append(output)

    return outputs 

    # run the simulation with the binary network 

def sim1step(formulas, inputstate, debug=False):
    for spe in inputstate.keys():
        command = str(spe) + ' = ' + str(inputstate[spe])
        exec (command)
    
    for term, biformula in formulas.items():
        # print(term, biformula)
        res = getResult(biformula, inputstate, debug)
        command = str(term) + '_tem_ = ' + str(res)
        exec(command)
    
    for term, biformula in formulas.items():
        command = str(term) + ' = ' + str(term) + '_tem_'
        exec (command)

    for spe in list(inputstate):
        command = 'inputstate[\'' + str(spe) + '\'] = ' + str(spe)
        # print (command)
        exec (command)
    return inputstate


    

def sim1bistep(biformulas, inputstate, debug=False):
    # run command to declare all the variables 
    for spe in inputstate.keys():
        # print(spe)
        command =  str(spe) + ' = ' + str(inputstate[spe])
        # print(command)
        exec(command)
    

    # run all the formulas with a copy of the variables 
    for formula in biformulas:
        # formula['right'].display()
        # command = formula['left'] + '_tem_ = ' + formula['right'] 
        res = getResult(formula['formula'], inputstate, debug)
        command = formula['term'] + '_tem_ = ' + str(res)
        exec(command)

    # assign new value to real variables taking from value of copy variables
    for formula in biformulas:
        command = formula['term'] + ' = ' + formula['term'] + '_tem_'
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
def getOutput(formulas, inputstate, isbi = False, maxStep = 10000, debug=False) -> dict: 
    oldstate = dict()
    numstep = 0

    # toshow is to test 
    # toshow = True
    # for spe, value in inputstate.items():
    #     if spe in ['TGFb', 'IFNg_e', 'TCR', 'IL6_e']:
    #         if not value:
    #             toshow = False
    #     else:
    #         if value:
    #             toshow = False 
    
    # toshow = False 

    for _ in range(maxStep):
        numstep += 1
        if isbi:
            inputstate = sim1bistep(formulas, inputstate, debug)
        else:
            inputstate = sim1step(formulas, inputstate, debug)
        hash = dict_hash(inputstate)


        if hash not in oldstate:
            oldstate[hash] = numstep
        else:
            # if debug:
            #     print("Number of iteration {} and first loop point is {}".format(numstep, oldstate[hash]))

            # merge all the state inside the loop 
            returnstate = copy.deepcopy(inputstate)
            for i in range(numstep-oldstate[hash] - 1):
                if isbi:
                    inputstate = sim1bistep(formulas, inputstate, debug)
                else:
                    inputstate = sim1step(formulas, inputstate, debug)
                returnstate = merge2states(returnstate, inputstate)
            # if toshow:
            #     print("Converge at step {}".format(numstep)) 
            #     print(returnstate)
            return returnstate 
    # print(type(inputstate))
    print("Cannot converge after {} steps".format(maxStep))
    return inputstate 

def getKnockoutOutput(formulas, inputstate, knockoutlist, isbi = False,  maxStep=1000, debug=False) -> dict:
    oldstate = dict()
    numstep = 0
    for _ in range(maxStep): 
        numstep += 1
        if isbi:
            inputstate = sim1bistep(formulas, inputstate, debug)
        else:
            inputstate = sim1step(formulas, inputstate, debug)
        # reset knockedout node = False 

        for node in knockoutlist:
            inputstate[node] = False

        hash = dict_hash(inputstate)
        
        if hash not in oldstate:
            oldstate[hash] = numstep
        else:
            if debug:
                print("Number of iteration {} and first loop point is {}".format(numstep, oldstate[hash]))

            # merge all the state inside the loop 
            returnstate = copy.deepcopy(inputstate)
            for i in range(numstep-oldstate[hash] - 1):
                if isbi:
                    inputstate = sim1bistep(formulas, inputstate, debug)
                else:
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

        # print (statistic)
        # print(statisticgen)
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


def calKnockoutShapleyForInput(genphe, inputnames, outputnames, v = False ,debug=False):
    ids = list([i for i in range(len(inputnames))])
    # print(ids)
    allsets = subsets(ids)

    # calculate for each output component
    shapss = dict()
    for outname in outputnames:
        shaps = dict()  
        # for each input component
        for inname in inputnames: # the knockout one 
            map = dict() 
            for i in range(len(inputnames)):
                map[i] = list(inputnames)[i]

            sum = 0
            for oneset in allsets:
                phenotem = set() 
                for e in oneset:
                    phenotem.add(map[e])
                pheintact = frozenset(phenotem)
                pheknockout = frozenset(phenotem - {inname})
                try: 
                    v_intact = 0
                    v_knockout = 0
                    if outname in genphe[pheintact]:
                        v_intact = 1
                    if outname in genphe[pheknockout]:
                        v_knockout = 1
                    
                    gain = v_intact - v_knockout
                    weightgain = gain * math.factorial(len(oneset)) * math.factorial(len(inputnames) - len(oneset))
                    sum += weightgain
                except:
                    continue
                
            shap = sum/(math.factorial(len(inputnames))) # this one is when divided by N! 
            # print(shap)
            shaps[inname] = round(shap, 3)
        shapss[outname] = shaps
    return shapss


def getBiOutputFromOriOutputs(orioutputs, biformulas, debug=False):
    print("Get psuedo binary output from ordinary output")
    bioutputs = []
    for output in orioutputs: 
        thisbioutput = copy.deepcopy(output)
        # print (output)
        for term, formula in biformulas.items():
            if term not in thisbioutput:
                value = getResult(formula, thisbioutput, debug) 
                thisbioutput[term] = value 
        bioutputs.append(thisbioutput)
            # print(term)
            # formula.display()  
    return bioutputs

# calculate original shapley value based on simulation result 
def calOriShapleyValue(genphe, inputnames, outputnames, v = False ,debug=False):
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
                
                # calculate gain, gain = v(pheyes) - v(pheno)
                # v(S) = 1 if outname is presence, v(S) = 0 if outname is not presence 
                # v_pheno = genphe[pheno] 
                try:
                    if not v:
                        if outname in genphe[pheno]: 
                            v_pheno = 1
                        else:
                            v_pheno = 0

                        if outname in genphe[pheyes]:
                            v_pheyes = 1
                        else:
                            v_pheyes = 0
                    else:
                        v_pheno = genphe[pheno]
                        v_pheyes = genphe[pheyes]

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
              
def genTableFromOutput(outputs, inputnames, debug=False):
    for line in outputs:
        # print(line)
        # each line is a dictionary with keys are name of species and value is true or false
        size = 0
        for input in inputnames:
            if line[input]:
                size += 1
        line['SIZE'] = size
        line['PROP'] = math.factorial(size)*math.factorial(len(inputnames) - size)/math.factorial(len(inputnames))
    if debug:
        print(outputs)
    return outputs


def main():
    parser = parseArguments()
    args= vars(parser.parse_args())
    if len(args) < 2:
        parser.print_help()
    else:
        debug = args['debug']
        isko = args['knockout']
        isbi = args['binary']
        isacyclic = args['acyclic']
        exp = args['expression']
        strouts = args['output'].split()
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
        
        # get the network of the original model
        orinet = convertBooleanFormulas2Network(formulas, inputnames, speciesnames, "network", debug)
        
        # the function returns the simulation output 
        orioutputs = workWithOriginalNetwork(inputnames, speciesnames, outputnames, internames, formulas, isko, debug)

        if isacyclic:
            # now do the limiting procedure for each output node 
            workwithLimitedNetwork(orinet, inputnames, internames, outputnames, speciesnames, False, formulas,debug)
        
        
        if isbi:
            print("------Now working with the binary network-------")
            workwithbinarynetwork(formulas, inputnames, outputnames, speciesnames, "binetwork", isko, debug)


def workwithLimitedNetwork(orinet, inputnames, internames, outputnames, speciesnames, isko, formulas, debug):
    print("Work with the minimum spanning tree extracted from the network with respect to output {}".format(outputnames))
    for outputname in outputnames:
        print("--------Now limit the network to output {}-------".format(outputname))
        temformulas = copy.deepcopy(formulas)
        limitGraphAfterNode(orinet, inputnames, outputname, temformulas, debug=True)
        workWithOriginalNetwork(inputnames, speciesnames, outputnames, internames, temformulas, isko, debug)


def workWithOriginalNetwork(inputnames, speciesnames, outputnames, internames, formulas, isko, debug):
    species = genSpecies(speciesnames, debug)
    inputstates = genInput(species, inputnames, debug)
    cloneinputstates = copy.deepcopy(inputstates)
    
    # now show the network 
    

    outputs = []
    for inputstate in cloneinputstates:
        output = getOutput(formulas, inputstate, False, 1000, debug)
        outputs.append(output)
    table = genTableFromOutput(outputs, inputnames)

    genphe = extractPhe(inputnames, outputnames, outputs)
    shapss = calKnockoutShapleyForInput(genphe, inputnames, outputnames)

    for outname in outputnames:
        outshap = 0
        for line in table:
            if line[outname]:
                outshap += line['PROP']
        shapss[outname][outname] = round(outshap,3)

    
    print("----Input shapley value-----")
    for item in shapss.items():
        print(item)   
        print('\n')
    
    if isko:
        print("-------Now perform the knockout procedure-------")
        
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
        for outname in outputnames:
            shaps = calknockoutShapleyValue(genphe, vs, outname, len(inputnames))
            print("----KNOCKOUT VALUE for output {}----".format(outname))
            print(shaps)
            print("\n")

    return outputs
     
def getknockoutlist(intercom, internames):
    knockouts = []
    for name in internames:
        if not intercom[name]:
            knockouts.append(name)
    return knockouts 

def calknockoutShapleyValue(intactgenphe, knockoutgenphes, output, numinput, debug=False):
    # print("Intact genphe")
    # print (intactgenphe)
    shaps = {}
    # for each intermediate node 
    for internode, genphes in knockoutgenphes.items():
        v = 0
        numrows = 0
        absgain = 0
        for gen, phe in genphes.items():
            # print(gen, phe)
            s = len(gen)
            gain = 0
            try:
                intactphe = intactgenphe[gen]
                # print("Intact phe")
                # print(intactphe)
                pheyes = 0
                pheno = 0
                try:
                    if output in intactphe:
                        pheyes = 1
                    if output in phe:
                        pheno = 1
                    gain = pheyes - pheno
                    if gain != 0:
                        numrows += 1
                    weightgain = gain * math.factorial(s) * math.factorial(numinput - s)
                    v += weightgain
                    absgain += abs(weightgain)
                except:
                    assert False, "Cannot find {} in the output set".format(output)
            except:
                assert False, "Cannot find the set {} for the intact network".format(gen)
        shaps[internode] = round(v/math.factorial(numinput),3)
        absgain = absgain/math.factorial(numinput)
        print("{} contributes to {} coalitions with the absolute value of {}".format(internode, numrows, absgain))
        if numrows != 0:
            print("The average is {} \n".format(round (absgain/numrows),3))
        else:
            print("")
    return shaps 


def getedges(left, root, net):
    if (root.right is None) and (root.left is None): # leave node  
        if root.parent and root.parent.val == "AND":
            net.add_edge (root.val, left, color='#008000', type='arrow', width=3) 
        elif root.parent and root.parent.val == "NOT":
            net.add_edge(root.val, left, color ="#FF0000", type='bar', width=3)
        else: 
            net.add_edge(root.val, left, color ='#808080', type='arrow', width=3)
        # print("Add edge {} {}".format(root.val, left))
    else:
        if root.right:
            getedges(left, root.right, net) 
        if root.left:
            getedges(left, root.left, net)
        


def graphicRepresentation(speciesnames, inputnames, outputnames, formulas, filename, shaps = None, debug=False):
    net = nx.DiGraph()
    for speciesname in speciesnames:
        if speciesname in outputnames:
            net.add_node(speciesname, label=speciesname, color='#FFFF00')
        elif speciesname in inputnames:
            net.add_node(speciesname, label=speciesname, color='#0000FF')
        else:
            net.add_node(speciesname, label=speciesname)
    for formula in formulas:
        root = formula['right'] 
        left = formula['left']

        getedges(left, root, net)

    nt = Network(directed=True)
    nt.from_nx(net)
    # nt.toggle_physics(False)
    nt.show(str(filename)+str(".html"), notebook=False)
    return net 

def findConstraints(net, output, formulas, debug = False):
    # convert formulas from list to dictionary
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
            print("Working with {}".format(cur))
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
                    print(com1, ":", filters[com1])
                    if cur in filters[com1]:
                        filters[com1].pop(cur)

                if com2 in filters:
                    print(com2, ":", filters[com2])
                    if cur in filters[com2]:
                        filters[com2].pop(cur)

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



# propagate function for binary tree 
def propagate(net, shaps, output, simtable, formulas, genphe, debug=False):
    
    filters = findConstraints(net, output, formulas, True)
    
    values = dict()
    for item in  shaps[output].items():
        # print(item)
        values[item[0]] = dict()
        values[item[0]]['lb'] = shaps[output][item[0]]
        values[item[0]]['ub'] = shaps[output][item[0]]
    # print(values)

    print("-----PROPAGATE------")
    # formulas should be ordered accordingly to layers and follow binary tree format 
    for formula in formulas:
        left = formula['left']
        if left != output:
            right = formula['right']

            neg = False 
            coms = []
            if right.val != "NOT":
                if right.right:
                    op = right.val
                    print("Operator is {}".format(op))
                else: 
                    coms.append(right.val)
            else:
                print ("Not supporting NOT operator at the moment")
                neg = True 
                return 

            # right.display()

            # get component 
            
            if not neg: # binary root but still can be unary components 
                if right.right: # binary
                    coms.append(right.right.val)
                    try:
                        coms.append(right.left.val)
                    except: 
                        print("Incomplete formula")
                else:
                    print ("There is no right component")
            
            if debug:
                print("Components are {}".format(coms))
            
            if len(coms) == 1: # the node to be propagated depends on only 1 node 
                outedges = net.out_edges([coms[0]])
                print (outedges)
                if len(outedges) == 1: 
                    print("Fully inherit")
                    try:
                        values[left] = values[coms[0]]
                    except:
                        values[left] = dict()
                        values[left]['ub'] = values[output]['ub']
                        values[left]['lb'] = 0
                else: # the parent node influences more than one node 
                        values[left] = dict()
                        values[left]['ub'] = values[output]['ub']
                        values[left]['lb'] = 0

            elif len(coms) == 2: # binary operator 
                outedges1 = net.out_edges(coms[0])
                outedges2 = net.out_edges(coms[1])
                print(outedges1)
                print(outedges2)
                if len(outedges1) == 1 and len(outedges2) == 1: # no branching from parent nodes 
                    if op == "AND": # case of AND 
                        if coms[0] in values:
                            values[left] = values[coms[0]]
                            if coms[1] not in values:
                                values[coms[1]] = values[coms[0]]
                        elif coms[1] in values:
                            values[left] = values[coms[1]]
                            values[coms[0]] = values[coms[1]]
                        else:
                            values[left] = dict()
                            values[left]['ub'] = values[output]['ub']
                            values[left]['lb'] = 0
                    elif op == "OR": 
                        values[left] = dict()
                        # print(values[output])

                        if coms[0] in values and coms[1] in values: # do normal OR propagate 
                            # this bound is super loose 
                            # values[left]['ub'] = values[output]['ub']
                            # values[left]['lb'] = values[coms[0]]['lb'] + values[coms[1]]['lb']  

                            # use constraints to compute value of left 
                            # firstly need to find the constraint for left 
                            constraints = filters[left]
                            leftshap = 0
                            for row in simtable:
                                if row[left]:
                                    satisfied = True
                                    for node, cons in constraints.items():
                                        if row[node] != cons:
                                            satisfied = False
                                    if satisfied:
                                        leftshap += row['PROP']
                            print("Pure OR operator, exact value of {} is {}".format(left, leftshap))
                            values[left]['ub'] = leftshap
                            values[left]['lb'] = leftshap

                        elif coms[0] in values: # only know value of one parent (coms[0])
                            values[left]['ub'] = values[output]['ub']
                            values[left]['lb'] = values[coms[0]]['lb'] 
                        elif coms[1] in values: # only know value of one parent (coms[1])
                            values[left]['ub'] = values[output]['ub']
                            values[left]['lb'] = values[coms[1]]['lb'] 
                        else: # no information about any parent 
                            values[left]['lb'] = 0
                            values[left]['ub'] = values[output]['ub']
                elif len(outedges1) + len(outedges2) == 3:
                    if op == "AND":
                        if len(outedges1) == 1:
                            values[left] = values[coms[0]]
                        else:
                            values[left] = values[coms[1]]
                    elif op == "OR":
                        values[left] = dict()
                        values[left]['ub'] = values[output]['ub']
                        if len(outedges1) == 1:
                            values[left]['lb'] = values[coms[0]]['lb']
                        else:
                            values[left]['lb'] = values[coms[1]]['lb'] 
                else:
                    values[left] = dict()
                    values[left]['lb'] = 0
                    if op == "AND":
                        values[left]['ub'] = min(values[coms[0]]['ub'], values[coms[0]]['ub'])
                    elif op == "OR":
                        values[left]['ub'] = values[output]['ub']
            else:
                print("Still hasn't support more than 2 branches")
    return values



# do everything with binary network (convert, get speciesnames, simulate...)          
def workwithbinarynetwork(formulas, inputnames, outputnames, orispeciesnames, networkname, isko = False, debug=False):
    bistrformulas = toBinaryFormulas(formulas, True)

    # to visualize the binary network
    binet = convertBiBooleanFormulas2Network(bistrformulas, inputnames, orispeciesnames, "bi" + networkname, debug) 
    
    biformulas = [] 
    bispeciesnames = set() 
    # for term, bistrformula in bistrformulas.items():
    for formula in bistrformulas:
        term = formula['term']
        bistrformula = formula['formula']
        thisfor = {'left': term, "right": bistrformula}

        thisbiformula = dict()
        thisbiformula['term'] = term
        thisbiformula['formula'] = parseFormula(thisfor, debug) 

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
        print("----Intermediate nodes-----")
        print(biinternames)
        print("----Intermediate nodes size is {}-----".format(len(biinternames)))
    # simulate binary network 
    realbioutputs = simBinaryNetwork(biformulas, inputnames, bispeciesnames, False)
    intactbigenphe = extractPhe(inputnames, outputnames, realbioutputs)
    if debug:
        print(intactbigenphe)
    
    bikoinshapss = calKnockoutShapleyForInput(intactbigenphe, inputnames, outputnames)
    print("-----Knockout Shapley value of input nodes in the binary network-----")
    for out, item in bikoinshapss.items():
        print(out, ":", item)

    
    if isko:
        # now do the knockout procedure with the binary network 
        print("-----Now do the knockout procedure with the binary network-----")

        # first generate all possible input states 
        inputstates = genInput(bispecies, inputnames, debug)
        vs = {} # to store genphe of knockout network 
        for internode in biinternames:
            if internode not in outputnames:
                print("Knockout {}".format(internode))
                inputstatescopy = copy.deepcopy(inputstates)
                kooutputs = []
                for inputstate in inputstatescopy:
                    output = getKnockoutOutput(biformulas, inputstate, [internode], True, 1000, debug)
                    kooutputs.append(output)
                    if debug:
                        print(output)
                
                genphe = extractPhe(inputnames, outputnames, kooutputs)
                vs[internode] = genphe

        for outname in outputnames:
            shaps = calknockoutShapleyValue(intactbigenphe, vs, outname, len(inputnames))
            print("---- Binary KNOCKOUT VALUE for output {}----".format(outname))
            print(shaps)
            print("\n")
               
            
if __name__ == "__main__":
    main()
