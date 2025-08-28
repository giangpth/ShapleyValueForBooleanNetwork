import copy 
from collections import deque
import random 

from Shapley.Node import Node
from Shapley.visualization import showNetwork 

def sim1step(formulas, inputstate, debug=False):
    """
    Perform one simulation step for a Boolean network.
    Args:
        formulas (dict): A dictionary with species names as keys and their boolean formulas as values.
        inputstate (dict): A dictionary with species names as keys and their current boolean states as values.
        debug (bool): If True, print debug information.
    Returns:
        dict: Updated inputstate after one simulation step.
    """
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
    """
    Perform one simulation step for a binary Boolean network.
    Different with sim1step, this function handle the extranodes added to binarize the network 
    Args:
        biformulas (list): A list of dictionaries, each containing 'term' and 'formula' keys for the binary formulas.
        inputstate (dict): A dictionary with species names as keys and their current boolean states as values.
        debug (bool): If True, print debug information.
        knockoutlist (list): A list of species names to be knocked out (set to False) or knocked in (set to True).
        isKnockin (bool): If True, species in knockoutlist are set to True; if False, they are set to False.
    Returns:
        dict: Updated inputstate after one simulation step.
    """
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


# get result of a boolean formula given the value of inputs as dictionary  
def getResult(curnode, inputs, debug=False):
    """
    Recursively evaluate the boolean formula represented by the tree rooted at `curnode`.   
    Parameters:
        curnode (Node): The current node in the boolean formula tree.
        inputs (dict): A dictionary mapping variable names to their boolean values.
        debug (bool): If True, print debug information.
    Returns:
        bool: The result of evaluating the boolean formula.
    """
    if not curnode:
        return None
    if not curnode.left and not curnode.right: # leaf node 
        try:
            return inputs[curnode.val]
        except:
            print("Cannot find {} in input list, return None".format(curnode.val))
            if debug:
                print(inputs)
            return None
    else:
        if not curnode.left: # not operator with no left 
            if curnode.val.upper() != "NOT":
                curnode.display()
            assert curnode.val.upper() == "NOT", "Uncomplete binary operator"
            # print("Not {}".format(cur.right.val))
            return (not getResult(curnode.right, inputs, debug))

        else: # may be nested, may be nuclear operator 
            if not curnode.right:
                curnode.display()
            assert curnode.right, "Uncomplete binary operator"
            if curnode.val.upper() == "AND":
                return (getResult(curnode.left, inputs, debug) and getResult(curnode.right, inputs, debug))
            elif curnode.val.upper() == "OR":
                return (getResult(curnode.left, inputs, debug) or getResult(curnode.right, inputs, debug))
            else:
                assert False, "Non support operator"

    
def find_parens(s):
    """
    Find matching parentheses in a list of tokens.
    Parameters:
        s (list): A list of tokens (strings).
    Returns:
        dict: A dictionary where keys are the indices of opening parentheses and values are the indices of corresponding closing parentheses.
    Raises:
        IndexError: If there are unmatched parentheses.
    """
    close = {}
    level = {}
    pstack = []

    for i, c in enumerate(s):
        if c == '(':
            pstack.append(i)
        elif c == ')':
            if len(pstack) == 0:
                raise IndexError("No matching closing parens at: " + str(i))
            cur = pstack.pop()
            close[cur] = i
            try:
                level[len(pstack)].append((cur, close[cur]))
            except:
                level[len(pstack)] = []
                level[len(pstack)].append((cur, close[cur]))

    if len(pstack) > 0:
        raise IndexError("No matching opening parens at: " + str(pstack.pop()))
    # print(close)
    return level

# def parseBinaryOperator(toprocess, debug=False):
#     assert len(toprocess) >= 1, "Empty parathesis {}".format(toprocess)
#     assert len(toprocess) <= 3, "More than 3 terms in binary operator {}".format(toprocess)
#     if len(toprocess) == 1: # identical assignment 
#         cur = Node(toprocess[0])  
#         if debug:
#             cur.display()
#     if len(toprocess) == 2: # NOT operator 
#         assert toprocess[0].upper() == "NOT", "Uncomplete binary operator {}".format(toprocess) 
#         cur = Node("NOT")
#         cur.insert_right(Node(toprocess[1]))
#         if debug:
#             cur.display()
#     if len(toprocess) == 3: # OR or AND operator 
#         assert toprocess[1].upper() in ["AND", "OR"], "Uncomplete binary operator {}".format(toprocess)
#         cur = Node(toprocess[1].upper())
#         cur.insert_left(Node(toprocess[0]))
#         cur.insert_right(Node(toprocess[2]))
#         if debug:
#             cur.display()
#     return cur 

# parse a nuclear term with no parathese, return a node representing this term                    
def parseNuclearTerm(toprocess, debug=False): 
    """
    Parse a nuclear term (a boolean expression without parentheses) into a binary tree representation.
    Parameters:
        toprocess (list): A list of tokens representing the boolean expression.
        debug (bool): If True, print debug information.
    Returns:
        Node: The root node of the binary tree representing the boolean expression.
    """
    assert len(toprocess) >= 1, "Empty parathesis"
    # print(toprocess)
    andcount = 0
    orcount = 0
    for word in toprocess:
        if type(word) == str:
            if word.upper() == "AND":
                andcount +=1
            if word.upper() == "OR":
                orcount += 1

    if andcount + orcount >= 2:
        if debug:
            print("More than 1 binary operator, dont allow NOT")
    if andcount >= 1 and orcount >= 1:
        assert False, "Non-homogenous binary operators {}".format(toprocess)

    biopelist = ["AND", "OR"]

    if len(toprocess) > 1: # not only 1 term but a nuclear operator 

        # in case of NOT A, create node with only right term
        if type(toprocess[0]) == str and toprocess[0].upper() == "NOT":
            assert len(toprocess) == 2, "Confusing term found: {}, , add parentheses to clarify scope of NOT operator".format((toprocess))
            node = Node("NOT")
            if type(toprocess[1]) == Node:
                node.insert_right(toprocess[1])
            else:
                node.insert_right(Node(toprocess[1]))
            return node 
        else:
            
            # case of A and/or not B ( and/or B ...)
            if type(toprocess[0]) == str:
                assert toprocess[0].upper() not in biopelist, "Stand alone binary operator found {}".format((toprocess))

            
            
            node = Node(toprocess[1].upper())

            if type(toprocess[0]) == Node:
                node.insert_left((toprocess[0]))
            else:
                node.insert_left(Node(toprocess[0]))

            assert len(toprocess) >= 3, "Uncomplete binary operator found {}".format((toprocess))
            
            # case of "A and/or not B"
            if type(toprocess[2]) == str and toprocess[2].upper() == "NOT":
                assert len(toprocess) == 4, "Confusing term found: {}, add parentheses to clarify scope of NOT operator".format((toprocess))
                morenode = Node("NOT") 

                if type(toprocess[3]) == Node:
                    morenode.insert_right((toprocess[3]))
                else:
                    morenode.insert_right(Node(toprocess[3]))

                node.insert_right(morenode)
                return node 
            else:
                if type(toprocess[2]) == Node:
                    node.insert_right((toprocess[2]))
                else:
                    node.insert_right(Node(toprocess[2]))
                
                idx = 3
                oldnode = node
                while idx < len(toprocess): # còn nước còn tát 
                    newnode = Node(toprocess[idx].upper())
                    newnode.insert_left(oldnode)
                    try:
                        if type(toprocess[idx+1]) == Node:
                            newnode.insert_right((toprocess[idx+1]))
                        else:
                            newnode.insert_right(Node(toprocess[idx+1]))
                    except:
                        assert False, "Uncomplete binary operator".format(toprocess)
                    oldnode = newnode
                    idx += 2
                return oldnode
                 
                    
            
    else: # case of only 1 term 
        if type(toprocess[0]) == Node:
            return toprocess[0]
        
        node = Node(toprocess[0])
        return node     

    return Node(None)            

def getFormula (lines, debug=False):
    """
    Parse boolean formulas from a list of lines.
    Args:
        lines (list): A list of strings, each representing a line from the input.
        debug (bool): If True, print debug information.
    Returns:
        list: A list of dictionaries, each containing 'left' and 'right' keys for the formulas.
    """
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

def testshow():
    print("----Test show----")
    node1 = Node("A")
    node2 = Node("B")
    node3 = Node("AND")
    node3.insert_left(node1)
    node3.insert_right(node2)

    node4 = Node("C")
    node5 = Node("D")
    node6 = Node("AND") 
    node6.insert_left(node4)
    node6.insert_right(node5)

    node7 = Node("OR") 
    node7.insert_left(node3)
    node7.insert_right(node6)

    node7.display()


def checkrow(net, rowraw, toconvergenode, targetnode, biformulas):
    """
    Check if knocking out the toconvergenode in the given row can change the state of the targetnode.
    Args:
        net (nx.DiGraph): The directed graph representing the network.
        rowraw (dict): A dictionary representing the state of all nodes in the network.
        toconvergenode (str): The name of the node to be knocked out.
        targetnode (str): The name of the target node to check for state change.
        biformulas (dict): A dictionary with species names as keys and their boolean formulas as values.
    Returns:
        bool: True if knocking out toconvergenode changes the state of targetnode, False otherwise.
    """
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

# function to refine formulas to exclude all the input nodes 
def parseFormula(formula, debug=False):
    """
    Parse a boolean formula represented as a string into a binary tree representation.
    Parameters:
        formula (dict): A dictionary with keys 'left' (the variable being defined) and 'right' (the boolean expression as a string).
        debug (bool): If True, print debug information.
    Returns:
        Node: The root node of the binary tree representing the boolean expression.
    """
    words = formula['right'].split() # split by space

    level = find_parens(words) # dictionary with keys are the opening parathesis, value of each key is the corresponding closing parathese

    while level:
        myKeys = list(level.keys())
        maxlev = max(myKeys)
        
        # print("Maximum level is {}".format(maxlev))

        terms = level[maxlev] # get all the terms at the smallest scope (maximum level) and convert it to nodes 
        # also change the list 'words' to replace old things with new nodes 

        offset = 0
        for term in terms: # for now term is the opening and closing index of the parathesis 
            # print(term)
            toprocess = words[term[0] + 1 - offset: term[1] - offset] # dont take parathesis 
            # print(toprocess) # only in the form "(not) A and/or (not) B"
            curnode = parseNuclearTerm(toprocess, debug)
            # curnode.display()

            # now replace the toprocess with the new node 
            words = words[:term[0] - offset] + [curnode] + words[term[1] + 1 - offset:] # replace everything inside paratheses (including parathesis)
            offset += len(toprocess) + 2 - 1
            # print("New list")
            # print(words)

        level = find_parens(words)
    
    # no more para, only nuclear term 
    root = parseNuclearTerm(words, debug)
    return root 

def recShrinkFormula2Node(node, nodesyntacs, biformulas):
    """
    Recursively shrink a boolean formula tree by replacing sub-expressions with temporary nodes.
    Parameters:
        node (Node): The current node in the boolean formula tree.
        nodesyntacs (dict): A dictionary mapping string representations of sub-expressions to their corresponding temporary node names.
        biformulas (dict): A dictionary mapping temporary node names to their corresponding boolean expressions as strings.
    Returns:
        str: The name of the temporary node representing the sub-expression rooted at `node`.
    """
    # print("Recurisve call")
    # node.display()
    assert node, print("Dont take None input")
    if node.left and node.right: # binary operator
        # add new node 
        if node.left.right: # leaf node 
            left = recShrinkFormula2Node(node.left, nodesyntacs, biformulas)
            # print(1)
        else:
            left = node.left.val 
            # print(2)
        # print(left)

        if node.right.right:
            right = recShrinkFormula2Node(node.right, nodesyntacs, biformulas) 
            # print(3)
        else:
            right = node.right.val 
            # print(4)
        # print(right)

        newnode1 = left + node.val + right 
        newnode2 = right + node.val + left # for communitative properties of boolean formulas 
        if newnode1 not in nodesyntacs and newnode2 not in nodesyntacs: 
        # if 1: # to test cycle 
            nodename = "_XTR_" + str(len(nodesyntacs) + 1) + '_'
            nodesyntacs[newnode1] = nodename 
            biformulas[nodename] = (left + " " + node.val + " " + right).strip() 
            # print(5)
            return nodename
        else: 
            try:
                return nodesyntacs[newnode1]
            except:
                print("No node with syntac {}, try {}".format(newnode1, newnode2))
            try: 
                return nodesyntacs[newnode2]
            except:
                print("No node with both syntac {} and {}".format(newnode1, newnode2))
                print("Check for error")
                exit()
    else:
        if not node.left: # uniary operator 
            if node.right:
                # print(7)
                if node.right.right:
                    # print(8)
                    right = recShrinkFormula2Node(node.right, nodesyntacs, biformulas)
                else:
                    # print(9)
                    right = node.right.val
            else: 
                # print(10)
                right = ""
            print(right)
            newnode = node.val + right 
            if newnode not in nodesyntacs: 
            # if 1: # to test cycle 
                nodename = "_XTR_" + str(len(nodesyntacs) + 1) + '_'
                nodesyntacs[newnode] = nodename 
                biformulas[nodename] = (node.val + " " + right).strip() 
                # print(11)
                return nodename 
            else: 
                # print(12)
                return nodesyntacs[newnode]
        else:
            # print(13)
            return node.val 


# return a dictionary of binary formulars 
# with keys are species (or temporary node) 
# and values are the corresponding biformulas 
def toBinaryFormulas(formulas, debug=False):
    """
    Convert a set of boolean formulas into an equivalent set of binary formulas using temporary nodes.
    Parameters:
        formulas (dict): A dictionary where keys are variable names and values are their corresponding boolean expressions as strings.
        debug (bool): If True, print debug information.
    Returns:
        list: A list of dictionaries, each containing 'term' (the variable name) and 'formula' (the binary boolean expression as a string).
    """
    print("----THIS ONLY WORKS FOR ACYCLIC NETWORK------")
    nodesyntacs = dict() 
    biformulas = dict() 
    todel = set() 
    
    for term, formula in formulas.items():
        # left = formula['left']
        # node = formula['right']
        left = term 
        node = formula
        if debug:
            print("{} = ".format(term))
            node.display()
        tem = recShrinkFormula2Node(node, nodesyntacs, biformulas)
        print("TEM: {} = {}".format(tem, biformulas[tem]))
        biformulas[left] = biformulas[tem]
        todel.add(tem)
        if debug:
            print()
    
    notdel = set()
    for term, biformula in biformulas.items():
        for tem in todel:
            if tem in biformula:
                if debug:
                    print("{} is needed for {}, do not del".format(tem, biformula))
                notdel.add(tem)
    
    
    for tem in todel:
        if tem not in notdel:
            del biformulas[tem]
    
    # # now convert string formula to node formula
    # for term, formula in biformulas.items():
    #     temfor = {'left': term, 'right': formula}
    #     biformulas[term] = parseFormula(temfor, debug)

    # if debug:
    #     print("Equivalent binary boolean network")
    #     for left, right in biformulas.items():
    #         print(left, "= ", right)
            # right.display()
    
    # rearrange biformulas so that all the extra added nodes are updated 
    # before the original nodes are update 
    # order is tem nodes that take the value from pure original nodes first 
    # then tem nodes that take value from other tem nodes plus odinary nodes
    # then tem nodes that take value from only tem nodes 
    # scan all the binary formulas to assign a level to sort it after this 

    nodelevel = dict() # firstly assign default level -1 for all the term
    for term, formular in biformulas.items():
        nodelevel[term] = -1
    
    tobeassigned = []

    tosort = []
    # at first scan the biformulas once to assign layer to node that can be assign 
    # (orignal nodes have layer 100 and extra nodes which depend on only original nodes have layer 0 
    for term, formula in biformulas.items():
        form = dict()
        if '_XTR_' in term:
            count = formula.count('_XTR_')
            if count == 0: # there is no extra node in the formulas 
                form['level'] = 0 # highest level 
                form['term'] = term 
                form['formula'] = formula
                tosort.append(form)
                nodelevel[term] = 0
            else:
                tobeassigned.append(term)
        else: 
            form['level'] = 100 # least level of the original node 
            form['term'] = term 
            form['formula'] = formula
            tosort.append(form)
            nodelevel[term] = 100

    while tobeassigned: # while to be assigned is not empty
        toprocess = tobeassigned.pop(0)
        # print("Trying to assign level to {}".format(toprocess))
        formula = biformulas[toprocess]

        coms = formula.split()
        numcoms = len(coms)
        if numcoms >= 3 : # this is a binary operator 
            first = coms[0] 
            second = coms[2] 
            if first in nodelevel and second in nodelevel: # none of them are input 
                if nodelevel[first] == -1 or nodelevel[second] == -1: # one of them hasnt had a level yet
                    tobeassigned.append(toprocess)
                    # print("1: cannot assign level to {}, put back in the list".format(toprocess))
                    continue
                else: # we can assigned some level to this 
                    form = dict() 
                    if max(nodelevel[first], nodelevel[second]) < 100:
                        form['level'] = 1 + max(nodelevel[first], nodelevel[second]) # least level of the original node 
                        nodelevel[toprocess] = 1 + max(nodelevel[first], nodelevel[second])
                    else:
                        form['level'] = 1 + min(nodelevel[first], nodelevel[second]) 
                        nodelevel[toprocess] = 1 + min(nodelevel[first], nodelevel[second])

                    form['term'] = toprocess 
                    form['formula'] = formula
                    tosort.append(form)

                    
            else:
                if first not in nodelevel: # first is input, second need to be extra node
                    if second in nodelevel and nodelevel[second] == -1:
                        tobeassigned.append(toprocess)
                        # print("2: cannot assign level to {}, put back in the list".format(toprocess))
                        continue
                    else:
                        form = dict() 
                        form['level'] = 1 + nodelevel[second] # least level of the original node 
                        form['term'] = toprocess 
                        form['formula'] = formula
                        tosort.append(form)

                        nodelevel[toprocess] = 1 + nodelevel[second] 
                else: # first is in nodelevel, second need to be the input 
                    if nodelevel[first] == -1:
                        tobeassigned.append(toprocess)
                        # print("3: cannot assign level to {}, put back in the list".format(toprocess))
                        continue
                    else:
                        form = dict() 
                        form['level'] = 1 + nodelevel[first] # least level of the original node 
                        form['term'] = toprocess 
                        form['formula'] = formula
                        tosort.append(form)

                        nodelevel[toprocess] = 1 + nodelevel[first] 


        elif numcoms == 2: # this is an uniary operator
            uni = coms[1] 
            if nodelevel[uni] == -1: 
                tobeassigned.append(toprocess)
                # print("4: cannot assign level to {}, put back in the list".format(toprocess))
                continue
            else:
                form = dict() 
                form['level'] = 1 + nodelevel[uni] # least level of the original node 
                form['term'] = toprocess 
                form['formula'] = formula
                tosort.append(form)

                nodelevel[toprocess] = 1 + nodelevel[uni]
            
        elif numcoms == 1: # this is an identical 
            uni = coms[0]
            if nodelevel[uni] == -1: 
                tobeassigned.append(toprocess)
                # print("5: cannot assign level to {}, put back in the list".format(toprocess))
                continue
            else:
                form = dict() 
                form['level'] = 1 + nodelevel[uni] # least level of the original node 
                form['term'] = toprocess 
                form['formula'] = formula
                tosort.append(form)

                nodelevel[toprocess] = 1 + nodelevel[uni]
        else:
            print("There is some problem with formula {}".format(formula))

    # now short the list tosort according to level 
    toreturn = sorted(tosort, key=lambda d: d['level'])
    # print(toreturn)
    if debug:
        for form in toreturn:
            print(form)

    print("Sorted equivalent binary formulas are:")
    for form in toreturn:
        print("{} = {}".format(form['term'], form['formula']))

    return toreturn 

def unrollCurNode(binet, cur, parents, outedges, biformulasdictandorder, node_positions, edgedatas, unrolldict, debug=False):
    print("Unrolling node {}".format(cur)) 
    print("{} has {} parents: {} and {} out-edges: {}".format(cur, len(parents), parents,  len(outedges), outedges))
    # get the formulas of the current node 
    try:
        oriformula = biformulasdictandorder[cur][0] # get the formula of the current node
        order = biformulasdictandorder[cur][1] # get the order of the current node
    except KeyError:
        print("Cannot find formula of node {}, maybe it has been removed already".format(cur))
        return 
    oriformula.display()

    try:
        curpos = node_positions[cur] 
    except:
        print("Cannot find position of node {}, use default (0, 0)".format(cur))
        curpos = (0, 0) 
    
    for j, edge in enumerate(outedges): 
        # new node 
        newnode = cur + '_' + str(j) 
        unrolldict[newnode] = cur # add the new node to the unrolldict, so that we can find it later

        newformula = dict() 
        newformula['term'] = newnode 
        newformula['formula'] = oriformula  
        
        # now handle the result node 
        resnode = edge[1]   
        
        try:
            childform = copy.deepcopy(biformulasdictandorder[resnode][0]) # get the formula of the result node 
            # print(biformulasdictandorder[resnode])
        except KeyError:
            print("Cannot find formula of result node {}, maybe it has been removed already".format(resnode))
            return

        # print("Working with child node: {} = ".format(resnode))
        # childform.display() 

        
        # modify the formula to replace the current node with the new node
        if childform.val == 'OR' or childform.val == 'AND':
            if childform.left.val == cur:
                childform.left.val = newnode
            elif childform.right.val == cur:
                childform.right.val = newnode
        else:
            if childform.val == 'NOT':
                if childform.right.val == cur:
                    childform.right.val = newnode
                else:
                    print("Cannot find current node {} in child formula".format(cur))
                    childform.display()
            else:
                if childform.val == cur:
                    childform.val = newnode
                else:
                    print("Cannot find current node {} in child formula".format(cur))
                    childform.display()

        print("New formula:")
        childform.display()

        # now modify the formula of resnode to include the new node 
        biformulasdictandorder[resnode] = (childform, biformulasdictandorder[resnode][1]) # add the new node and its formula to the dictionary
        # now add the new node to the biformulasdictandorder
        biformulasdictandorder[newnode] = (oriformula, order)

        # print("Add new node {} with formula:".format(newnode))
        # biformulasdictandorder[newnode][0].display()

        # also need to rewire the network 
        # first connect parents to the new node and then remove the edge from parents to current node
        xpos = curpos[0] + random.randint(-20, 20)
        ypos = curpos[1] + random.randint(-20, 20)
        binet.add_node(newnode, x=xpos, y=ypos) # add new node with position
        node_positions[newnode] = (xpos, ypos) # update the position of the new node


        
        for parent in parents:
            edgedata = edgedatas[(parent, cur)] # get the edge data of the current node to the parent
            # print("Edge data from {} to {} is {}".format(parent, cur, edgedata))
            print("Add edge: {} {}".format(parent, newnode))
            binet.add_edge(parent, newnode, color=edgedata['color'], type='arrow', width=3) # keep the color of the edge 
            edgedatas[(parent, newnode)] = edgedata # update the edge data for the new edge
            try:
                binet.remove_edge(parent, cur) # remove the edge from parent to current node
                print("Remove edge: {} {}".format(parent, cur))
            except:
                print("Cannot remove edge from {} to {}, maybe it has been removed already".format(parent, cur))
                pass
            
        # now add the edge from the new node to the resnode also remove the edge from cur to resnode
        otheredgedata = edgedatas[(cur, resnode)] # get the edge data of the current node to the result node
        # print("Edge data from {} to {} is {}".format(cur, resnode, otheredgedata))
        binet.add_edge(newnode, resnode, color=otheredgedata['color'], type='arrow', width=3) 
        print("Add edge: {} {}".format(newnode, resnode))
        edgedatas[(newnode, resnode)] = otheredgedata # update the edge data for the new edge
        try:
            binet.remove_edge(cur, resnode) # remove the edge from current node to result node
            print("Remove edge: {} {}".format(cur, resnode))
        except:
            print("Cannot remove edge from {} to {}, maybe it has been removed already".format(cur, resnode))
            pass
        print('-----')
    return 

def unrollBinaryNetwork(binet, outname, biformulasdictandorder, node_positions, debug=False):
    """
    Unroll a binary network by duplicating nodes with multiple outgoing edges.
    Parameters:
        binet (networkx.DiGraph): The binary network to be unrolled.
        outname (str): The name of the output node.
        biformulasdictandorder (dict): A dictionary where keys are node names and values are tuples of (formula, order).
        node_positions (dict): A dictionary mapping node names to their positions (x, y).
        debug (bool): If True, print debug information.
    Returns:
        tuple: A tuple containing:
            - unrolldict (dict): A dictionary mapping new node names to their original node names.
            - biformulasdictandorder (dict): The updated dictionary of node formulas after unrolling.
    """
    # starting from output  
    # biformulasdict is just to quickly access the formula of a node
    # biformulaslist is to keep the order of the formulas
    # First of all create the index mapping for the term and its formula in the list first 
    unrolldict = dict()
    edgedatas = dict()
    unrolllist = set() 
    for edge in binet.edges(data=True):
        edgedatas[(edge[0], edge[1])] = binet.get_edge_data(edge[0], edge[1]) # store the edge data for later use


    # print(biformulaslist)
    print("--------UNROLLING--------")
    curs = [outname] 
    while curs:
        print("\n\n---Processing layers of {}---".format(curs))
        nextlayer = []
        for cur in curs:
            try:
                form = biformulasdictandorder[cur][0]
            except: 
                print("\nReach node {} without in-comming edges".format(cur))
                continue
            inedges = list(binet.in_edges(cur))
            assert len(inedges) <= 2, print("Support only binary network, but get {} in-edges for node {}".format(len(inedges), cur))
            if len(inedges) == 2: 
                print(f"{cur} ==== {form.left.val} {form.val} {form.right.val}")

                if form.left.val not in nextlayer:
                    nextlayer.append(form.left.val)
                if form.right.val not in nextlayer:
                    nextlayer.append(form.right.val)
                
                outedges = list(binet.out_edges(cur)) 
                if len(outedges) > 1:
                    unrolllist.add(cur) # add this node to the unroll list
                    unrollCurNode(binet, cur, [form.left.val, form.right.val], outedges, biformulasdictandorder, node_positions, edgedatas, unrolldict, debug)
                            

            elif len(inedges) == 1:
                if form.val == 'NOT':
                    singlemom = form.right.val
                    print(f"{cur} === NOT {singlemom}")
                else:
                    singlemom = form.val
                    print(f"{cur} === {singlemom}")
                if singlemom not in nextlayer:
                    nextlayer.append(singlemom)
                
                outedges = list(binet.out_edges(cur)) 
                if len(outedges) > 1:
                    unrolllist.add(cur) # add this node to the unroll list
                    unrollCurNode(binet, cur, [singlemom], outedges, biformulasdictandorder, node_positions, edgedatas, unrolldict, debug)
            # now unroll this current node if it has branches 
            if cur == outname:
                print("Output node, no need to unrolling")
                continue
        showNetwork(binet, None, None, None, None, 'unroll.html')
        curs = nextlayer


    # now delete all the nodes that are not in the unroll list
    for node in unrolllist:
        binet.remove_node(node) # remove the node from the network
        del biformulasdictandorder[node] # remove the node from the biformulasdictandorder

    # now visualize the network
    print("Unrolling finished, now visualize the network")
    showNetwork(binet, None, None, None, None, 'unroll.html')
    
    return unrolldict, biformulasdictandorder 
    

def expandFormula(biformulas, target, debug=False):
    """
    Expand a target boolean formula by replacing sub-expressions with their definitions from a set of binary formulas.
    Parameters:
        biformulas (dict): A dictionary where keys are variable names and values are their corresponding boolean expressions as strings.
        target (str): The variable name of the target formula to be expanded.
        debug (bool): If True, print debug information.
    Returns:
        Node: The root node of the binary tree representing the expanded boolean expression.
    """
    targetfunctioncoms = biformulas[target].split(' ') 
    assert len(targetfunctioncoms) <= 3, "Too many component in binary operator {}".format(targetfunctioncoms) 
    assert len(targetfunctioncoms) >= 1, "Too few component in binary operator {}".format(targetfunctioncoms) 
    # root = Node(targ)
    leaves = deque()
    if len(targetfunctioncoms) == 3: # binary operator
        if debug:
            print("Binary operator") 
        assert targetfunctioncoms[1].upper() in ["AND", "OR"], "Strange binary operator {}".format(targetfunctioncoms[1])
        root = Node(targetfunctioncoms[1].upper())
        root.insert_left(Node(targetfunctioncoms[0]))
        root.insert_right(Node(targetfunctioncoms[2]))
        # queue.append(root) # add the operator   
        leaves.append(root.left)
        leaves.append(root.right)

    elif len(targetfunctioncoms) == 2: 
        if debug:
            print("NOT operator")
        assert targetfunctioncoms[0].upper() == "NOT", "Expect NOT operator but get {}".format(targetfunctioncoms)
        root = Node(targetfunctioncoms[0].upper())
        # queue.append(root) # add the operator
        root.insert_right(Node(targetfunctioncoms[1])) 
        leaves.append(root.right) 
    else: # only 1 term 
        if debug:
            print("Only 1 term, identical operator")
        root = Node(targetfunctioncoms[0])
        leaves.append(root) # add the operator 
    while leaves: # while there is still some leaf node
        cur = leaves.popleft() 
        if cur.val in biformulas: 
            if debug:
                print("Expand node {} with {}".format(cur.val, biformulas[cur.val]))
            thiscoms = biformulas[cur.val].split(' ')
            assert len(thiscoms) <= 3, "Too many component in binary operator {}".format(thiscoms) 
            assert len(thiscoms) >= 1, "Too few component in binary operator {}".format(thiscoms) 
            if len(thiscoms) == 3: # binary operator
                assert thiscoms[1].upper() in ["AND", "OR"], "Strange binary operator {}".format(thiscoms[1])
                cur.oldVal = cur.val 
                cur.val = thiscoms[1].upper()
                left = Node(thiscoms[0])
                left.oldVal = left.val 
                right = Node(thiscoms[2]) 
                right.oldVal = right.val 
                cur.insert_left(left)
                cur.insert_right(right) 
                leaves.append(cur.left)
                leaves.append(cur.right)
            elif len(thiscoms) == 2:    
                assert thiscoms[0].upper() == "NOT", "Expect NOT operator but get {}".format(thiscoms)
                cur.oldVal = cur.val
                cur.val = thiscoms[0].upper()
                right = Node(thiscoms[1])
                right.oldVal = right.val 
                cur.insert_right(right) 
                leaves.append(cur.right)
            else: # only 1 term
                cur.oldVal = cur.val
                cur.val = thiscoms[0]
                leaves.append(cur) # add the operator
            
    return root

def expandFunction(biformulas, target, debug=False):
    """
    Expand a target boolean formula by replacing sub-expressions with their definitions from a set of binary formulas.
    Parameters:
        biformulas (dict): A dictionary where keys are variable names and values are their corresponding boolean expressions as strings.
        target (str): The variable name of the target formula to be expanded.
        debug (bool): If True, print debug information.
    Returns:
        Node: The root node of the binary tree representing the expanded boolean expression.
    """
    targetfunction = biformulas[target].split(' ') 
    change = True
    while change:
        oldtarget = copy.deepcopy(targetfunction) 
        for i, term in enumerate(targetfunction):
            if term in biformulas:
                change = True
                toreplace = biformulas[term].split(' ') 
                toreplace.insert(len(toreplace), ')') 
                for j in reversed(toreplace):
                    targetfunction.insert(i+1, j) 
                targetfunction[i] = '('

        if oldtarget == targetfunction:
            change = False
    
    targetstrformula = dict()
    targetstrformula['left'] = target
    targetstrformula['right'] = " ".join(targetfunction)
    if debug:
        print(targetstrformula['right'])
    targetformula = parseFormula(targetstrformula, debug) 
    # targetformula.display()
    return targetformula 


def setLevelTargetFunction(targetfunction, debug=False):
    """
    Set the level of each node in the binary tree representing a boolean function using BFS traversal.
    Parameters:
        targetfunction (Node): The root node of the binary tree representing the boolean function.
        debug (bool): If True, print debug information.
    Returns:
        list: A list of levels of nodes in the order they were visited during BFS traversal.
    """
    # ******* BFS ********
    if not targetfunction:
        return []

    targetfunction.level = 0 # Root node is at level 0
    queue = deque([targetfunction])  # Initialize queue with root node
    bfs_result = []  # List to store traversal order

    while queue:
        node = queue.popleft()  # Dequeue the front node
        node.level = node.parent.level + 1 if node.parent else 0 # Set level of current node
        bfs_result.append(node.level)  # Process node 
        if debug:
            print(node.level, node.val)
        # Enqueue left child if exists
        if node.left:
            queue.append(node.left)

        # Enqueue right child if exists
        if node.right:
            queue.append(node.right)
    # ******* END BFS ******** 

#dicvals: a dictionary with keys including "level", "id" and "Valid"
#level: the level of the current node 
#id: the id of the possible value for current node 0: (0,0); 1: (1,0); 2: (0,1); 3: (1,1) 
#valid: the validity of the current assignment 
def setValueForNode (node, value, dictvals, level, id, carryon, debug=False):
    """
    Recursively set values for nodes in a binary tree representing a boolean function and track valid assignments.
    Parameters:
        node (Node): The current node in the binary tree, to be assigned a value.
        value (int): The value to set for the current node (0 or 1).
        dictvals (dict): A dictionary to store valid assignments for different identifiers.
        level (int): The current level in the binary tree.
        id (int): The identifier for the current assignment (0 to 3).
        carryon (list): A list to track the path of identifiers leading to the current node.
        debug (bool): If True, print debug information.
    Returns:
        None
    """
    if not node:
        return 
    carryon[level] = str(id) 
    if debug:
        node.level = value
        node.display(showlevel=True)
        print("Carry id: {}".format(['0','1','2','3','4','5','6','7']))
        print("Carry on: {}".format(carryon)) 
        
    if not node.left and not node.right: # leaf node 
        # node.exVal = value 
        if debug:
            print("Level {}, id {}: Set value for node {} to {}".format(level, id, node.val, value))
        
        identifier = ''.join(carryon)
        if identifier not in dictvals: 
            dictvals[identifier] = dict()
            dictvals[identifier]['Valid'] = True 
        
        if node.val in dictvals[identifier]: 
            if dictvals[identifier][node.val] != value:
                print("ID {} : Conflict at node {} with value {} and {}".format(identifier, node.val, dictvals[identifier][node.val], value))
                dictvals[identifier]['Valid'] = False
        else:
            print("ID {} : Set value for node {} to {}".format(identifier, node.val, value))
            dictvals[identifier][node.val] = value 

    else:
        if not node.left: # not operator with no left 
            assert node.val.upper() == "NOT", "Uncomplete binary operator"
            setValueForNode(node.right, abs(value - 1), dictvals, level+1, 0, carryon, debug)
        else: # may be nested, may be nuclear operator 
            if node.val == "AND":
                if value == 1:
                    setValueForNode(node.left, 1, dictvals, level+1, 3, carryon, debug)
                    setValueForNode(node.right, 1, dictvals, level+1, 3, carryon, debug)
                else:
                    setValueForNode(node.left, 0, dictvals, level+1, 1, carryon, debug)
                    setValueForNode(node.right, 1, dictvals, level+1, 1, carryon, debug) 
                    print("New possibliity")
                    setValueForNode(node.left, 1, dictvals, level+1, 2, carryon, debug)
                    setValueForNode(node.right, 0, dictvals, level+1, 2, carryon, debug) 
                    print("New possibliity")
                    setValueForNode(node.left, 0, dictvals, level+1, 0, carryon, debug)
                    setValueForNode(node.right, 0, dictvals, level+1, 0, carryon,debug)
            elif node.val == "OR":
                if value == 1:
                    setValueForNode(node.left, 1, dictvals, level+1, 2, carryon, debug)
                    setValueForNode(node.right, 0, dictvals, level+1, 2, carryon ,debug)
                    print("New possibliity")
                    setValueForNode(node.left, 0, dictvals, level+1, 1, carryon, debug)
                    setValueForNode(node.right, 1, dictvals, level+1, 1, carryon, debug) 
                    print("New possibliity")
                    setValueForNode(node.left, 1, dictvals, level+1, 3, carryon, debug)
                    setValueForNode(node.right, 1, dictvals, level+1, 3, carryon, debug)
                else:
                    setValueForNode(node.left, 0, dictvals, level+1, 0, carryon, debug)
                    setValueForNode(node.right, 0, dictvals, level+1, 0, carryon ,debug)

# targetfunction: a node representing the target function 
def propagateFromTarget(table, targetfunction, debug=False):
    
    print("---------Propagate from target---------")
    setLevelTargetFunction(targetfunction, debug)
    if debug:
        print("Target function")
        # targetfunction.display(True)
        # targetfunction.display(showlevel=True)
    # dictvals = dict()
    # setValueForNode(targetfunction, 1, dictvals, 0, 1, ['_', '_', '_','_','_','_','_','_'],  debug)   
    # if debug:
    #     for id, vals in dictvals.items():
    #         print("{}: {}".format(id, vals))
        

    


