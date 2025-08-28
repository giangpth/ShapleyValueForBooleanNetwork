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

def checkrow(net, rowraw, toconvergenode, targetnode, biformulas):
    # do BFS 
    row = copy.deepcopy(rowraw)
    bfslayers = dict(enumerate(nx.bfs_layers(net, sources=toconvergenode)))
    intactnode = row[toconvergenode]
    intactout = row[targetnode] 
    row[toconvergenode] = not intactnode 
    for layer, nodes in bfslayers.items():
        for node in nodes:
            if node == toconvergenode:
                continue
            try: 
                nodeform = biformulas[node] 
            except:
                print(f"Cannot find formula for node {node}")
            nodeafter = getResult(nodeform, row) 
            row[node] = nodeafter
            if node == targetnode:
                if nodeafter == intactout:
                    return False
                else:
                    return True 
    return False

def diamondDigger(net, outname, formulas):
    # first round, locate all the diamonds and their endpoints, starting from output node 
    curs = [outname]
    carryonofnodes = dict() 

    while curs:
        print(f'\n\n Processing layer of {curs}')
        nextlayer = []
        for cur in curs: 
            try:
                form = formulas[cur]
            except:
                print(f"Encounter node {cur} without incoming edges")
            inedges = list(net.in_edges(cur))
            assert len(inedges) <= 2, print("Support only binary network")

            if len(inedges) == 2: # care only about binary operator 
                left = form.left.val 
                right = form.right.val 
                # first, carryon  current node cur 
                if left not in carryonofnodes:
                    carryonofnodes[left] = dict() 
                
                if cur not in carryonofnodes[left]:
                    carryonofnodes[left][cur] = []
                carryonofnodes[left][cur].append('L')


                if right not in carryonofnodes:
                    carryonofnodes[right] = dict() 

                if cur not in carryonofnodes[right]:
                    carryonofnodes[right][cur] = []
                carryonofnodes[right][cur].append('R')

                if cur in carryonofnodes:
                    print(f"Now passing carry on of current node {cur} to its parents {left} and {right}")
                    # distribute the carryon of current node to its parents 
                    for node, carryon in carryonofnodes[cur].items():
                        # carryonofnodes[left][node] = copy.deepcopy(carryonofnodes[cur][node]) 
                        # print(node, carryon)
                        if node not in carryonofnodes[left]:
                            carryonofnodes[left][node] = []
                        else:
                            print("Node {} already in carry on of {}".format(node, left))
                        
                        if node not in carryonofnodes[right]:
                            carryonofnodes[right][node] = []
                        else:
                            print("Node {} already in carry on of {}".format(node, right)) 

                        for side in carryon:
                            carryonofnodes[left][node].append(side + 'L')
                            carryonofnodes[right][node].append(side + 'R')
                        
                        # carryonofnodes[right][node].extend(carryon)
                        # carryonofnodes[right][node].append('R')
                        # for side in carryonofnodes[left][node]: 
                        #     # print(side)
                        #     carryonofnodes[left][node].append('L')

                        # carryonofnodes[right][node] = copy.deepcopy(carryonofnodes[cur][node]) 
                        # for side in carryonofnodes[right][node]: 
                        #     # print(side)
                        #     carryonofnodes[right][node].append('R')

                if left not in nextlayer:
                    nextlayer.append(left)
                if right not in nextlayer:
                    nextlayer.append(right)
                
                print("Carryon of {:20}: {}".format(left,carryonofnodes[left]))
                print("Carryon of {:20}: {}".format(right,carryonofnodes[right]))
            elif len(inedges) == 1:
                singlemon = inedges[0][0] 

                #just pass all the carry on of current node cur to its single mom 
                if singlemon not in carryonofnodes:
                    carryonofnodes[singlemon] = dict()

                if cur in carryonofnodes:
                    print(f"Now passing carry on of current node {cur} to its singlemon {singlemon}")
                    for node, carryon in carryonofnodes[cur].items():
                        if node not in carryonofnodes[singlemon]:
                            carryonofnodes[singlemon][node] = []
                        else:
                            print(f"Node {node} already in carry on of node {singlemon}")
                        
                        for side in carryon:
                            carryonofnodes[singlemon][node].append(side)
                        

                if singlemon not in nextlayer:
                    nextlayer.append(singlemon)
                print("Carryon of {:20}: {}".format(singlemon,carryonofnodes[singlemon]))
            else:
                print(f"Node {cur} has no input, pass") 
        curs = nextlayer

    diamonds = dict()
    print ("-----------------Now check DIAMONDS from the resulted carryons--------------------------")
    for node, carryon in carryonofnodes.items():
        print(f"Carryon of {node}")
        if len(net.out_edges(node)) > 1:
            for one, onelist in carryon.items():
                print("{:20}: {}".format(one, onelist))
                L, R = False, False
                for astr in onelist:
                        tostop = False 
                        side = astr[0] 
                        if side == 'L':
                            L = True 
                        if side == 'R':
                            R = True
                if L and R:
                    print("\tNode {} should inherit the carryon rows of node {}".format(node, one))
                    if node not in diamonds:
                        diamonds[node] = set()
                    diamonds[node].add(one)

    return diamonds


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