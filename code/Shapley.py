import argparse as ap
import itertools
from typing import Dict, Any
import hashlib
import json
import math
import copy 
from booleanFormulaHandler import parseFormula, getResult, toBinaryFormulas, expandFunction, propagateFromTarget, expandFormula
from networkHandler import convertBiBooleanFormulas2Network, convertBooleanFormulas2Network, manipulateNetwork, showNetwork
import networkx as nx
from pyvis.network import Network
from datetime import datetime 
import timeit
import random 
import shap 
import numpy as np
from collections import defaultdict

def rank_dict_values(d):
    # Sort items by value descending
    sorted_items = sorted(d.items(), key=lambda x: -x[1])
    
    ranks = {}
    current_rank = 1
    last_value = None

    for key, value in sorted_items:
        if value != last_value:
            last_value = value
            ranks[key] = current_rank
            current_rank += 1
        else:
            ranks[key] = current_rank - 1  # Same as previous rank
    
    return ranks

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
    parser.add_argument('-i', '--knockin', help='Perform knockin or not', \
                        action='store_true')
    parser.add_argument('-b', '--binary', help='Work with binary network or not', \
                        action='store_true')
    parser.add_argument('-a', '--acyclic', help='Extract acyclic network with respect to output node', \
                        action='store_true')
    parser.add_argument('-p', '--propagate', help='Run propagation for binary network or not', \
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

# species is a dictionary with keys are species, values are the state of the species 
# this function convert a state of the Boolean network to a decimal number 
def getOrderedList(inputnames, internames, debug=False):
    sortedinput = sorted(list(inputnames))
    if debug:
        print('Sorted input names')
        print(sortedinput)
        print("There are {} inputs with the range of input: {} {}".format(len(sortedinput),0,2**(len(sortedinput)) - 1))
    sortedinter = sorted(list(internames))
    if debug:
        print('Sorted internames:')
        print(sortedinter)
        print("There are {} intermediate nodes with the range of inter: {} {}".format(len(sortedinter), 0,2**(len(sortedinter)) - 1))
    return sortedinput, sortedinter 

def toDecimal(state, sortedinput, sortedinter, debug=False):
    try:
        inputbitstr = ''
        for input in sortedinput:
            if state[input]:
                inputbitstr += '1'
            else:
                inputbitstr += '0'
        try:
            interbitstr = ''
            for inter in sortedinter:
                if state[inter]:
                    interbitstr += '1'
                else:
                    interbitstr += '0'
            
            intinputs = int(inputbitstr, 2)
            intinters = int(interbitstr, 2)
            return intinputs, intinters 
        except:
            print("Cannot find intermediate node {} in the state".format(inter))
    except:
        print("Cannot find input {} in the state".format(input))


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

# test the correctness of the function toBinaryNetwork      
def simBinaryNetwork(biformulas, inputnames, speciesnames, sortedinput, sortedinter, debug=False, maxstep=1000, extranodes=None):
    print("----Simulate binary network----")
    species = genSpecies(speciesnames, debug) # take list the names of all the species and set them to be None 
    inputstates = genInput(species, inputnames, debug) # each inputstate contains all species, input species are set to a certain value 

    # if debug:
        # print("Integer version of state of Binary network")
    decimalPairs = dict() 
    if extranodes:
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


    

def sim1bistep(biformulas, inputstate, debug=False, knockoutlist=None, isKnockin=False):
    # run command to declare all the variables 
    for spe in inputstate.keys():
        # print(spe)
        command =  str(spe) + ' = ' + str(inputstate[spe])
        # print(command)
        exec(command)
    
    if knockoutlist:
        for knockoutnode in knockoutlist:
            # print("Setting value of {} to False".format(knockoutnode))
            if isKnockin:
                inputstate[knockoutnode] = True
            else: 
                inputstate[knockoutnode] = False 

    # run all the formulas with a copy of the variables 
    for formula in biformulas:
        # formula['right'].display()
        # command = formula['left'] + '_tem_ = ' + formula['right'] 
        res = getResult(formula['formula'], inputstate, debug)

        if "_XTR_" not in formula['term']: # for original node, run all the formulas then assign the value back after this 
            command = formula['term'] + '_tem_ = ' + str(res)
            # print(command)
            exec(command)
        else: # but for extra node added to binarize the network, assign value immediately, also update the value of '_XTR_' node in inputstate 
            command = formula['term'] + ' = ' + str(res) 
            # print(command)
            exec(command)

            updateinputstatecommand = 'inputstate[\'' + formula['term'] + '\'] = ' + str(res)
            # print(updateinputstatecommand)
            exec(updateinputstatecommand)
            if knockoutlist:
                if formula['term'] in knockoutlist:
                    # print("Dont update value of {}, set to False/True".format(formula['term']))
                    # set the value of this node to be False/True since it is knocked out/in 
                    if isKnockin:
                        inputstate[formula['term']] = True
                    else:
                        inputstate[formula['term']] = False 
            

    # assign new value to real variables taking from value of copy variables
    for formula in biformulas:
        if "_XTR_" not in formula['term']: # now assign the value for the original node 
            command = formula['term'] + ' = ' + formula['term'] + '_tem_'
            exec(command)
            # print(command)  
    
    # put back to inputstate to return 
    for spe in list(inputstate):
        command = 'inputstate[\'' + str(spe) + '\'] = ' + str(spe)
        # print(command)
        exec(command)

    return inputstate

def merge2states(state1, state2, debug=False):
    for item in state1.keys():
        # print(state1[item], state2[item])
        state1[item] = state1[item] or state2[item]
    return state1

# inputstate is a list of all species with the predefine value for input species
# this function is the simulation process 
def getOutput(formulas, inputstate, isbi = False, maxStep = 10000, debug=False, extranodes=None) -> dict: 
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
                    except:
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
                    except:
                        print("Cannot find the root of the extra node {}".format(extranode)) 
                        return

            inputstate = sim1step(formulas, inputstate, debug)

            # for extranode in extranodes:
                # print("Extranode {} get value meanwhile rootnode have value {}")


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
                    if extranodes:
                        for extranode in extranodes:
                            try:
                                rootname = extranode.split("_to_")[0] 
                                if rootname not in inputstate:
                                    print("There is no node named {} as the root node of {}".format(rootname, extranode))
                                inputstate[extranode] = inputstate[rootname]
                            except:
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
                            except:
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

def getKnockoutOutput(formulas, inputstate, knockoutlist, isbi = False,  \
                      maxStep=1000, debug=False, extranodes = None, isKnockin=False) -> dict:
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
                    except:
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
                    except:
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
                            except:
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
                            except:
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

def subsets(s):  
    if len(s) == 0:  
        return [[]]  
    x = subsets(s[:-1])  
    return x + [[s[-1]] + y for y in x] 

def extractPhe(inputnames, outputnames, outputstates, debug=False):
    genphe = dict() 
    statistic = dict()
    statisticgen = dict()

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
    return genphe


def calKSV4Input(genphe, inputnames, outputnames, knockin=False, simtable = None, debug=False):
    countedrows = dict()
    ids = list([i for i in range(len(inputnames))])
    # print(ids)
    allsets = subsets(ids)

    # calculate for each output component
    shapss = dict()
    for outname in outputnames:
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
                    if outname in genphe[pheintact]:
                        v_intact = 1
                    if outname in genphe[pheknockout]:
                        v_knockout = 1
                    
                    gain  = v_intact - v_knockout
                    weightedGain = gain * math.factorial(len(oneset)) * math.factorial(len(inputnames) - len(oneset))
                    sum += weightedGain
                    if weightedGain != 0:
                        numrow += 1
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

                except:
                    continue
            shap = round(sum/(math.factorial(len(inputnames))),4)
            shaps[inname] = shap
        shapss[outname] = shaps
    if simtable:
        return shapss, countedrows
    else:
        return shapss

              
def genTableFromOutput(simoutputs, inputnames, sortedinput, sortedinter, outputnames, debug=False):
    index = dict() # for each node, save the IDs of row that the node is TRUE
    aindex = dict() # for each node, save the IDs of row that the node is TRUE
    dictresult = dict()
    for inter in sortedinter:
        index[inter] = set()
        aindex[inter] = set()
    for input in inputnames:
        index[input] = set()
        aindex[input] = set()
    for outputname in outputnames:
        index[outputname] = set() 
        aindex[outputname] = set() 
    
    # print("Integer verson of output:")
    for id, line in enumerate(simoutputs):
        # inp, inter = toDecimal(line, sortedinput, sortedinter)
        # print(inp, inter)
        # print(line)
        # each line is a dictionary with keys are name of species and value is true or false
        
        size = 0
        for input in inputnames:
            if line[input]:
                size += 1
        line['SIZE'] = size
        line['PROP'] = round(math.factorial(size)*math.factorial(len(inputnames) - size)/math.factorial(len(inputnames)),4)
        for inter in sortedinter:
            if line[inter]:
                index[inter].add(id)
            else:
                aindex[inter].add(id)
        for input in inputnames:
            if line[input]:
                index[input].add(id)
            else:
                aindex[input].add(id)
        for outputname in outputnames: 
            if line[outputname]:
                index[outputname].add(id)
            else:
                aindex[outputname].add(id)

        dictresult[id] = line
    return dictresult, index, aindex 


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



# a wrapper for SHAP, taking a dataset of samples and computes the output of the model for those samples
class BNmodel:
    def __init__(self, map, species, formulas, outputname):
        self.map = map # list of input in order 
        self.species = species 
        self.formulas = formulas 
        self.outputname = outputname
    
    def predict(self, inputs): # input is 2 dimension numpy array 
        # print(inputs.shape)
        outputs = np.zeros((inputs.shape[0],1))
        # outputs = []
        # print(outputs.shape)

        # print(len(outputs))
        # print(len(inputs[0]))
        for i in range(inputs.shape[0]):
            # print(inputs[i])
            inputstate = copy.deepcopy(self.species)
            for j in range(inputs.shape[1]):
                if inputs[i][j] == 1:
                    inputstate[self.map[j]] = True
                if inputs[i][j] == 0:
                    inputstate[self.map[j]] = False
            # print(inputstate)
            output = getOutput(self.formulas, inputstate, False, 1000, False)
            # print(output)
            if output[self.outputname]:
                outputs[i][0] = 1
            else:
                outputs[i][0] = 0
            # if output[self.outputname]:
            #     outputs.append(1)
            # else:
            #     outputs.append(0)
        return outputs

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


def genRandomInputForSHAP(numfeature, numsample):
    inputs = []
    for i in range(numsample):
        randint = random.randint(0, pow(2, numfeature) - 1) 
        code = '{0:0' + str(numfeature) + 'b}'
        bistr = code.format(randint)
        # print(bistr)
        # oneinput = np.array(list(bistr))
        inputs.append(bistr)
    
    binary_array = np.array([list(map(int, list(s))) for s in inputs])
    # print(binary_array.shape)

    return binary_array
    # return inputs 

def genAllInputForSHAP(numfeature):
    inputs = []
    for i in range(pow(2,numfeature)):
        code = '{0:0' + str(numfeature) + 'b}'
        bistr = code.format(i)
        # print(bistr)
        # oneinput = np.array(list(bistr))
        inputs.append(bistr)
    
    binary_array = np.array([list(map(int, list(s))) for s in inputs])
    # print(binary_array.shape)

    return binary_array
    # return inputs 

def identity_masker(data):
    return data

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
     

def calKSV(intactgenphe, knockoutgenphes, outputnames, numinput, simtable=None, inputnames=None, debug=False):
    countedrow = dict()
    shapss = dict()
    for output in outputnames:
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
                if output in intactphe:
                    pheyes = 1
                if output in phe:
                    pheno = 1
                gain = pheyes - pheno

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
                            # if internode == 'CI_p' or internode == 'ci':
                            #     print(knockoutgenphes[internode])
                # except:
                    # assert False, "Cannot find {} in the output set".format(output)
                # except:
                    # assert False, "Cannot find the set {} for the intact network".format(gen)
            print("Number of rows for {} is {}".format(internode, numrow))
            shaps[internode] = round(sum/math.factorial(numinput),4)
        shapss[output] = shaps
    if simtable:
        return shapss, countedrow
    else:
        return shapss, None

        
def findConstraints(net, output, formulas, debug = False):
    print("\nFinding constraints for output {}".format(output))
    # convert formulas from list to dictionary
    # formulas pass here is a list of string
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
            # print("Working with {}".format(cur))
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
                    # print(com1, ":", filters[com1])
                    if cur in filters[com1]:
                        filters[com1].pop(cur)

                if com2 in filters:
                    # print(com2, ":", filters[com2])
                    if cur in filters[com2]:
                        filters[com2].pop(cur)
                
                # if fordict[cur] is a string, convert to Node 
                if isinstance(fordict[cur], str):
                    temdict = dict()
                    temdict['right'] = fordict[cur]
                    temdict['left'] = cur
                    fordict[cur] = parseFormula(temdict, debug) 

                            
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
    else:
        print("Do not support {} operator".format(op))
    
    return allrows
    
                                      
def processBranching(net, toconvergenode, carryon, countedrowsofnodes, rowsofpairs, resofpairs, index, aindex, formulas, extranodes, convergedpairs):
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
            print(f"Rows of pair {pair} that really account to converge node {toconvergenode} are {orrows}")
            # add this row to the value of the node  
            if toconvergenode not in countedrowsofnodes:
                countedrowsofnodes[toconvergenode] = set()
            countedrowsofnodes[toconvergenode].update(orrows)
        print("After BRANCHING, the row for node {} is {}".format(toconvergenode, countedrowsofnodes[toconvergenode]))
    else:
        print("No carryon pair to converge to node {}".format(toconvergenode))
    
def processBranching_v3(net, toconvergenode, carryon, countedrowsofnodes, rowsofpairs, resofpairs, index, aindex, formulas, extranodes, convergedpairs):
    if toconvergenode in carryon:
        des = nx.descendants(net, toconvergenode)
        desincluded = copy.deepcopy(des)
        desincluded.add(toconvergenode)
        print("---BRANCHING---\nNode {} has branching having desendants {}".format(toconvergenode, desincluded))
        countedpairs = set()

        print("Checking convergence for the pairs {}".format(carryon[toconvergenode]))
        for pair in carryon[toconvergenode]:
            if pair[0] in des and pair[1] in des:
                countedpairs.add(pair) 

        if len(countedpairs) > 0:
            print("PAIRS TO CONVERGED: \n{}".format(countedpairs))
        else:
            print("NO PAIR TO CONVERGE, continue")
            return

        for pair in countedpairs:
            if toconvergenode not in convergedpairs:
                convergedpairs[toconvergenode] = set()
            if pair in convergedpairs[toconvergenode]:
                print("Already counted the pair {}, continue".format(pair))
                continue
            convergedpairs[toconvergenode].add(pair)
            print("ROWS carried by pair {}: \n {}".format(pair, sorted(list(rowsofpairs[pair]))))
            # now get resulted node and operator
            resnode = resofpairs[pair] 
            rootform = formulas[resnode]
            rootop = rootform.val
            tocheck = [resnode] 
            alreadychecked = set()
            countedrows = rowsofpairs[pair] # get rows counted by the pair, this time is only to pass to enrichfilter to test
            while tocheck:
                node_ = tocheck.pop(0)
                if node_ in alreadychecked:
                    continue
                alreadychecked.add(node_)
                if node_ == toconvergenode:
                    continue

                try:
                    form = formulas[node_]
                except:
                    print("Pass the input node {}".format(node_))
                    continue

                op = form.val
                if op == 'OR' or op == 'AND': # binary operator 
                    left = form.left.val
                    right = form.right.val
                    if left in desincluded and right in desincluded:
                        print("BOTH parents {} and {} of node {} are in descendants".format(left, right, node_)) 
                        # if operator is or, both parents need to be checked 
                        if op == 'OR':
                            if left not in alreadychecked:
                                tocheck.append(left)
                            if right not in alreadychecked:
                                tocheck.append(right)
                        # if operator is and, only one parent involve in a next AND operator need to be checked
                        if op == 'AND':
                            try:
                                leftform = formulas[left]
                                leftop = leftform.val
                                if leftop == 'AND':
                                    if left not in alreadychecked:
                                        tocheck.append(left)
                            except:
                                print("Pass the input node {}".format(left))
                                continue
                            try: 
                                rightform = formulas[right] 
                                rightop = rightform.val
                                if rightop == 'AND':
                                    if right not in alreadychecked:
                                        tocheck.append(right)
                            except:
                                print("Pass the input node {}".format(right)) 

                    elif left in desincluded and right not in desincluded:
                        print("Parent {} of node {} is not in descendants".format(right, node_))
                        if op == 'OR':
                            countedrows = enrichfilter(op, right, node_, index, aindex, extranodes, countedrows, countedrows) 
                        # if left no
                    elif left not in desincluded and right in desincluded:
                        print("Parent {} of node {} is not in descendants".format(left, node_))
                    else:
                        print("BOTH parents {} and {} of node {} are NOT in descendants".format(left, right, node_))
                else:
                    print("Unary operator")
                    print("Add node {} to the list of tocheck nodes".format(node_))
                    if node_ not in alreadychecked:
                        tocheck.append(node_)
                    continue



def  processBranching_v2(net, toconvergenode, carryon, countedrowsofnodes, rowsofpairs, resofpairs, index, aindex, formulas, extranodes, convergedpairs):
    if toconvergenode in carryon:
        des = nx.descendants(net, toconvergenode)
        desincluded = copy.deepcopy(des)
        desincluded.add(toconvergenode)
        print("BRANCHING: Node {} has branching having desendants {}".format(toconvergenode, des))
        countedpairs = set()

        print("Checking convergence for the pairs {}".format(carryon[toconvergenode]))
        for pair in carryon[toconvergenode]:
            if pair[0] in des and pair[1] in des:
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
            convergedpairs[toconvergenode].add(pair)
            rowsoffactors = dict()
            print("Rows carried by pair {} are {}".format(pair, sorted(list(rowsofpairs[pair]))))
            for factor in pair:
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
                                # tofilter.add((form.left.val, op))
                                if form.right.val not in checked:
                                    tocheck.append(form.right.val)
                                # if node_ not in path:
                                #     indirect = True
                                # else:
                                #     indirect = False
                                initrows = enrichfilter(op, form.left.val, node_, index, aindex, extranodes, node_, initrows)
                                if len(initrows) == 0:
                                    print("EMPTY ROWS, stop!")
                                    break
                            elif form.right.val not in desincluded and form.left.val in desincluded:
                                print('Parent {} of node {} is not in descendants'.format(form.right.val, node_))
                                if form.left.val not in checked:
                                    tocheck.append(form.left.val)
                                # if node_ not in path:
                                #     indirect = True
                                # else:
                                #     indirect = False
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
            print("Rows counted by the pair {} are {}".format(pair, sorted(list(rowstoconverge))))
            countedrowsofnodes[toconvergenode].update(rowstoconverge)
            print("--------")
    else:
        print("Carry nothing to converge")

def rowstovalues(rowdict, simtable, outname):
    kovalues = dict()
    kivalues = dict()
    for node, rowids in rowdict.items():
        kovalues[node] = 0.0
        kivalues[node] = 0.0
        for id in rowids:
            row = simtable[id]
            if row[outname]:
                if row[node]:
                    kovalues[node] += row['PROP']
                else:
                    kivalues[node] += row['PROP']
            else:
                if row[node]:
                    kovalues[node] -= row['PROP']
                else:
                    kivalues[node] -= row ['PROP']
    return kovalues, kivalues


def propagateBottomUp(net, simtable, index, aindex, outname, formulas, extranodes):
    curs = [outname] 
    countedrowsofnode = dict() # visited[node] = set(), rows that make node count 
    rowsofpairs = dict() # save the loss information of the OR operator 
    carryonofnodes = dict() # carryonofnodes[node] = set(), pairs that are carried on to the next layer
    resofpairs = dict() # the node that is the result of pairs 
    converged = dict() # save the pairs that are already converged for a node (key) to avoid doing it again 
    descendantsofnodes = dict() # save the descendants of a node to check if the node is already counted
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
            print("\nCounted row for current node {} is {}".format(cur, sorted(list(countrows))))

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
                    processBranching_v2(net, form.left.val, carryonofnodes, countedrowsofnode, rowsofpairs, resofpairs, index, aindex, formulas, extranodes, converged)

                if len(list(net.out_edges(form.right.val))) > 1:
                    # print(list(net.out_edges(form.right.val)))
                    processBranching_v2(net, form.right.val, carryonofnodes, countedrowsofnode, rowsofpairs, resofpairs, index, aindex, formulas, extranodes, converged)
                
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
                        processBranching_v2(net, form.right.val, carryonofnodes, countedrowsofnode, rowsofpairs, resofpairs, index, aindex, formulas, extranodes, converged)
                   
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
                        processBranching_v2(net, form.val, carryonofnodes, countedrowsofnode, rowsofpairs, resofpairs, index, aindex, formulas, extranodes, converged)

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
            rowsofnodes = propagateBottomUp(binet, table, index, aindex, outname, biformulasdict, extranodes)
            # correct the rows for input nodes 
            for input in inputnames:
                rowsofnodes[input] = koinrows[input].union(kiinrows[input])
            # propko, propki = rowstovalues(rowsofnodes, table, outname)
        
            ########################################
            ### now perform the top down propagation
            # for the top down need to calculate the rows for extranodes with no input
            # now calculate the knockout and knockin shapley value for the EXTRANODES
            # inputstates = genInput(bispecies, inputnames, False)
            # vsexko = {} # to store genphe of knockout network
            # vsexki = {} # to store genphe of knockin network
            # for enode in extranodes:
            #     print("Knockout {}".format(enode))
            #     inputstatescopy = copy.deepcopy(inputstates)
            #     kooutputs = []
            #     for inputstate in inputstatescopy:
            #         output = getKnockoutOutput(biformulas, inputstate, [enode], True, 1000, False, extranodes)
            #         kooutputs.append(output)

            #     genphe = extractPhe(inputnames, outputnames, kooutputs)
            #     vsexko[enode] = genphe

            #     print("KnockIN {}".format(enode))
            #     inputstatescopy2 = copy.deepcopy(inputstates)
            #     kioutputs = []
            #     for inputstate in inputstatescopy2:
            #         output = getKnockoutOutput(biformulas, inputstate, [enode], True, 1000, False, extranodes, isKnockin=True)
            #         kioutputs.append(output)

            #     genphe = extractPhe(inputnames, outputnames, kioutputs)
            #     vsexki[enode] = genphe
            # koexshaps, koexrows = calKSV(intactbigenphe, vsexko, outputnames, len(inputnames), table, inputnames)
            # kiexshaps, kiexrows = calKSV(intactbigenphe, vsexki, outputnames, len(inputnames), table, inputnames)
            # # merge ko and ki rows for inputs 
            # inrows = dict() 
            # for input in inputnames:
            #     if input in koinrows and input in kiinrows:
            #         inrows[input] = koinrows[input].union(kiinrows[input])
            #     elif input in koinrows:
            #         print("NOT FOUND {} in knockin list".format(input))
            #         inrows[input] = koinrows[input]
            #     elif input in kiinrows:
            #         print("NOT FOUND {} in knockout list".format(input))
            #         inrows[input] = kiinrows[input]
            #     else:
            #         print("NOT FOUND {} in both lists".format(input))
            #         inrows[input] = set()
            # for enode in extranodes:
            #     if enode in koexrows and enode in kiexrows:
            #         inrows[enode] = koexrows[enode].union(kiexrows[enode])
            #     elif enode in koexrows:
            #         print("NOT FOUND {} in knockin list".format(enode))
            #         inrows[enode] = koexrows[enode]
            #     elif enode in kiexrows:
            #         print("NOT FOUND {} in knockout list".format(enode))
            #         inrows[enode] = kiexrows[enode]
            #     else:
            #         print("NOT FOUND {} in both lists".format(enode))
            
            # rowsofnodes = topDown(binet, biformulasdict, rowsofnodesBU, inrows, nodes_layer)
            # ### end of top down propagation 
            # #############################################
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

        
        # if korows and kirows:
        #     for input in inputnames:
        #         if input in rowsofnodes: 
        #             gt = koinrows[input].union(kiinrows[input])
        #             lack = gt.difference(rowsofnodes[input])
        #             extra = rowsofnodes[input].difference(gt)
        #             if len(lack) > 0 or len(extra) > 0:
        #                 print(input)
        #                 print("Lack  {:5}: {}".format(len(lack), sorted(list(lack))))
        #                 print("Extra {:5}: {}".format(len(extra), sorted(list(extra))))
        #                 print("Real  {:5}: {}".format(len(gt), sorted(list(gt))))
        #                 print("Prop  {:5}: {}\n".format(len(rowsofnodes[input]), sorted(list(rowsofnodes[input]))))

        #     for inter, row in korows.items():
        #         if inter in rowsofnodes:
        #             gt = row.union(kirows[inter])
        #             lack = gt.difference(rowsofnodes[inter])
        #             extra = rowsofnodes[inter].difference(gt)
        #             # if len(lack) > 0 or len(extra) > 0:
        #             if True:
        #                 print(inter)
        #                 print("Lack  {:5}: {}".format(len(lack), sorted(list(lack))))
        #                 print("Extra {:5}: {}".format(len(extra), sorted(list(extra))))
        #                 print("Real  {:5}: {}".format(len(gt), sorted(list(gt))))
        #                 print("Prop  {:5}: {}\n".format(len(rowsofnodes[inter]), sorted(list(rowsofnodes[inter]))))

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
