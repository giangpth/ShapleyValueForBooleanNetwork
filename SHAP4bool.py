import argparse as ap
from ast import parse
from datetime import date
from logging import exception
import os
from datetime import date
import itertools
from typing import Dict, Any
import hashlib
import json
import math
import shap
import numpy as np 


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

# species is a dictionary, inputnames is a list 
def genInput(species, inputnames, debug=False):
    numin = len(inputnames)
    coms = list(itertools.product([0, 1], repeat=numin))
    if debug:
        print("There are {} combinations for input".format(len(coms)))
        # print (coms)
    inputstates = []
    for com in coms:
        inp = species.copy()
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
    # teminp = species.copy()
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
        command = formula['left'] + '_tem_ = ' + formula['right'] 
        exec(command)


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
    return inputstate

def merge2states(state1, state2, debug=False):
    for item in state1.keys():
        # print(state1[item], state2[item])
        state1[item] = state1[item] or state2[item]
    return state1

# inputstate is a list of all species with the predefine value for input species
# this function is the simulation process 
def getOutput(formulas, inputstate, debug=False, maxStep = 10000) -> dict: 
    # print(inputstate)

    # for spe in inputstate.keys():
    #     # print(spe)
    #     command =  str(spe) + ' = ' + str(inputstate[spe])
    #     # print(command)
    #     exec(command)
    
    oldstate = dict()
    numstep = 0
    for _ in range(maxStep):
        numstep += 1
        inputstate = sim1step(formulas, inputstate, debug)
        hash = dict_hash(inputstate)
        # print("Step {}".format(numstep))
        # print(inputstate)

        if hash not in oldstate:
            oldstate[hash] = numstep
        else:
            if debug:
                print("Number of iteration {} and first loop point is {}".format(numstep, oldstate[hash]))

            # merge all the state inside the loop 
            returnstate = inputstate.copy()
            for i in range(numstep-oldstate[hash] - 1):
                inputstate = sim1step(formulas, inputstate, debug)
                returnstate = merge2states(returnstate, inputstate)
            return returnstate 
    # print(type(inputstate))
    return inputstate 

def getVanilaOutput(genphe, inputs, outputnames, inputnames):
    outputs = []
    for one in inputs:
        gene = set()
        for id, name in enumerate(inputnames):
            if one[id] == 1:
                gene.add(name)
        frozengene = frozenset(gene)
        phe = genphe[frozengene]
        output = []
        for id, name in enumerate(outputnames):
            if name in phe:
                output.append(1)
            else:
                output.append(0)
        outputs.append(output)
    return np.asarray(outputs)



class boolnetwork:
    def __init__(self, speciesnames, inputnames, outputnames, formulas, maxsteps = 1000, debug=False):
        self.formulas = formulas
        self.maxsteps = maxsteps
        self.speciesnames = speciesnames
        self.inputnames = inputnames
        self.species = genSpecies(speciesnames, debug)
        self.outputnames = outputnames 
        print(inputnames)
        print(outputnames) 


    # input is just a simple list or vector of {0, 1} with respect to a predefined order of genes 
    # input will be converted to input states (with species and their presence/absence)
    def simulate(self, inputs):
        outputstates = []
        for input in inputs:
            inp = self.species.copy()
            for id, name in enumerate(self.inputnames):
                try:
                    if input[id] == 1:
                        inp[name] = True
                    elif input[id] == 0:
                        inp[name] = False
                    else:
                        inp[name] = input[id]
                except KeyError:
                    print("There is no key named {}".format(name))
            out = getOutput(self.formulas, inp, maxStep=self.maxsteps)
            outputstates.append(out)
        genphe = extractPhe(self.inputnames, self.outputnames, outputstates)
        # print (genphe)
        outputs = getVanilaOutput(genphe, inputs, self.outputnames, self.inputnames)
        return outputs

def subsets(s):  
    if len(s) == 0:  
        return [[]]  
    x = subsets(s[:-1])  
    return x + [[s[-1]] + y for y in x] 

def extractPhe(inputnames, outputnames, outputstates, debug=False):
    genphe = dict() 
    statistic = dict()
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
    print (statistic)

    return genphe

def calOriShapley(genphe, inputnames, outputnames, debug=False):
    ids = list([i for i in range(len(inputnames) - 1)])
    # print(ids)
    allsets = subsets(ids)
    # print(len(allsets))
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
            # print (map)
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
            shap = sum/(math.factorial(len(inputnames)))
            # print(shap)
            shaps[inname] = round(shap, 3)
        # print(shaps)
        shapss[outname] = shaps
    return shapss

def main():
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
        speciesnames = getSpeciesName(lines, debug)
        if not outputnames.issubset(speciesnames):
            print("The list of interested output is not a subset of species")
            return -1
        species = genSpecies(speciesnames, debug)
        print("---SPECIES---")
        print(species)
        inputnames = getInputNames(lines, speciesnames, debug)
        inputstates = genInput(species, inputnames, debug)
        all = getFormula(lines, debug)

        boolnet = boolnetwork(speciesnames=speciesnames, 
                              inputnames=inputnames, 
                              outputnames=outputnames, 
                              formulas=all)
        inputs = np.asarray(list(itertools.product([0, 1], repeat=len(inputnames))))
        
        # inputs = [[0, 0, 0, 1, 1], [0, 0, 1, 0, 1]]
        # outputs = boolnet.simulate(input)
        # print(outputs)
        explainer = shap.Explainer(boolnet.simulate, inputs)
        shap_values = explainer(inputs)
        # print(shap_values)
        # shap.plots.bar(shap_values)


        # outputs = []
        # for inputstate in inputstates:
        #     # print (inputstate)
        #     output = getOutput(all, inputstate, debug, 1000)
        #     outputs.append(output)
        #     # print(output)
        #     # print()

        # genphe = extractPhe(inputnames, outputnames, outputs)
        # for tem in genphe.items():
        #     print(tem)

        # shapss = calOriShapley(genphe, inputnames, outputnames)
        # # print(shapss)
        # for item in shapss.items():
        #     print(item)        

            

if __name__ == "__main__":
    main()