import itertools
import copy
from random import random, sample, randint

def getSpeciesName(lines, debug=False):
    """
    Extract unique species names from a list of lines.
    Args:
        lines (list): A list of strings, each representing a line from the input.
        debug (bool): If True, print the extracted species names.
    Returns:
        set: A set of unique species names.
    """
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

def initSpeciesStates(speciesnames, debug=False):
    """Initialize species states to False.
    Args:
        speciesnames (set): A set of species names.
        debug (bool): If True, print the initialized species states.
    Returns:
        dict: A dictionary with species names as keys and False as values.
    """
    speciesState = dict()
    for s in speciesnames:
        speciesState[s] = False
    return speciesState 

def getInputNames(lines, speciesnames, debug=False):
    """
    Extract input species names from a list of lines. 
    Input species are those that do not appear on the left side of any equation.

    Args:
        lines (list): A list of strings, each representing a line from the input.
        speciesnames (set): A set of all species names.
        debug (bool): If True, print the extracted input species names.
    Returns:
        set: A set of input species names.
    """
    dep = set() 
    for line in lines:
        coms = line.split('=') # split into left and right side
        for tem in coms[0].split():
            if '-' in tem: # this is for excluding '-' in the formulas which is confusing for the compiler 
                tem = tem.replace('-','_')
            dep.add(tem)
    input = speciesnames - dep
    if debug:
        print("----GET INPUT----")
        print(dep)
        print(input)
    return input

#species is a dictionary, inputnames is a list 
def genInput(species, inputnames, samplingStrategy='all', samplePercentage=60, debug=False):
    """
    Generate all possible input states for the given species and input names.
    Args:
        species (dict): A dictionary with species names as keys and their boolean states as values.
        inputnames (set): A set of input species names.
        samplingStrategy (str): Can be: "all" (default), "uniform", "antithetic" or "stratified" .
        samplePercentage (int): The percentage of combinations to sample if using a sampling strategy different from "all". Default is 60.
        debug (bool): If True, print debug information.
    Returns:
        list: A list of dictionaries, each representing a unique input state with all species.
    """
    print("Generating input states using strategy: {}".format(samplingStrategy))
    if samplingStrategy != 'all':
        if samplePercentage <= 0 or samplePercentage > 100:
            raise ValueError("IN GENERATE INPUT STATES: samplePercentage must be between 0 and 100.")
        print("Sample percentage: {}%".format(samplePercentage))
    
    numin = len(inputnames) 
    inputnamelist = sorted(list(inputnames)) # sort to have a fixed order of input names

    print(f"{numin} inputs causing in total {2**numin} combinations")
    print("Input list:")
    print(inputnamelist)

    if samplingStrategy == 'all':
        coms = list(itertools.product([0, 1], repeat=numin))

    elif samplingStrategy == 'uniform':
        numSamples = int((samplePercentage/100) * (2**numin))
        random_numbers = sample(range(2**numin), numSamples)

        coms = [list(map(int, bin(num)[2:].zfill(numin))) for num in random_numbers]
    elif samplingStrategy == 'antithetic':
        numSamples = int((samplePercentage/100) * (2**numin))
        half_samples = numSamples // 2

        random_numbers = sample(range(2**numin), half_samples)
        antithetic_numbers = [(2**numin - 1) - num for num in random_numbers]

        combined_numbers = random_numbers + antithetic_numbers
        coms = [list(map(int, bin(num)[2:].zfill(numin))) for num in combined_numbers]
    elif samplingStrategy == 'stratified':
        numSamples = int((samplePercentage/100) * (2**numin))
        strata_size = 2**numin // numSamples

        random_numbers = []
        for i in range(numSamples):
            start = i * strata_size
            end = start + strata_size - 1
            random_num = randint(start, end)
            random_numbers.append(random_num)

        coms = [list(map(int, bin(num)[2:].zfill(numin))) for num in random_numbers]
    else:
        raise ValueError("IN GENERATE INPUT STATES: Invalid sampling strategy. Choose from 'all', 'uniform', 'antithetic', or 'stratified'.")

    coms = sorted(coms, key=lambda x: int(''.join(map(str, x)), 2))  # Sort combinations by their decimal value

    inputstates = []
    mappingrow2inputs = dict()
    for comid, com in enumerate(coms):
        value = int(''.join(map(str, com)), 2)
        inp = copy.deepcopy(species)
        mappingrow2inputs[value] = [] 
        for id, name in enumerate(inputnamelist):
            try:
                # print(name, id)
                if com[id] == 1:
                    inp[name] = True
                    mappingrow2inputs[value].append(name)
                elif com[id] == 0:
                    inp[name] = False
                else:
                    inp[name] = com[id]
            except KeyError:
                print("IN GENERATE INPUT STATES: There is no key named {}".format(name))
        inputstates.append(inp)
    print("Number of samples in total: {}".format(len(inputstates)))

    # UNCOMMNENT THE FOLLOWING LINES TO PRINT EACH GENERATED INPUT STATE
    # for i, state in enumerate(inputstates):
    #     print(f"Sample {i+1}: {state}")

    return inputstates, mappingrow2inputs