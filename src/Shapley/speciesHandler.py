import itertools
import copy

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
def genInput(species, inputnames, debug=False):
    """
    Generate all possible input states for the given species and input names.
    Args:
        species (dict): A dictionary with species names as keys and their boolean states as values.
        inputnames (set): A set of input species names.
        debug (bool): If True, print debug information.
    Returns:
        list: A list of dictionaries, each representing a unique input state with all species.
    """
    numin = len(inputnames) 
    inputnamelist = sorted(list(inputnames)) # sort to have a fixed order of input names
    print("INPUT NAMES LIST:")
    print(inputnamelist)
    coms = list(itertools.product([0, 1], repeat=numin))
    if debug:
        print("There are {} combinations for input".format(len(coms)))
        # print (coms)
    inputstates = []
    for com in coms:
        inp = copy.deepcopy(species)
        # print(inp)
        for id, name in enumerate(inputnamelist):
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