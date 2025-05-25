import copy 
from Node import Node
from collections import deque


# get result of a boolean formula given the value of inputs as dictionary  
def getResult(curnode, inputs, debug=False):
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



# function to refine formulas to exclude all the input nodes 
def parseFormula(formula, debug=False):
    # print('\n')
    # print(formula)
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

    if debug:
        print("Equivalent binary boolean network")
        for left, right in biformulas.items():
            print(left, "= ", right)
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

    return toreturn 

# def recReplace (root, cur, biformulas, debug=False):
#     # print(cur)
#     if cur == None:
#         print("None node")
#     if cur.left:
#         # print("Replace 1")
#         if cur.left.val in biformulas:
#             leafterm = parseBinaryOperator(biformulas[cur.left.val].split(' '), debug)
#             cur.left = leafterm 
#             recReplace(cur.left, biformulas, debug)
#     if cur.right:
#         # print("Replace 2")
#         if cur.right.val in biformulas:
#             rightterm = parseBinaryOperator(biformulas[cur.right.val].split(' '), debug)
#             cur.right = rightterm
#             recReplace(cur.right, biformulas, debug)
#     if not cur.left and not cur.right:
#         # print("Replace 3")
#         if cur.val in biformulas:
#             thisterm = parseBinaryOperator(biformulas[cur.val].split(' '), debug) 
#             cur.val = thisterm 
#             recReplace(cur.val, biformulas, debug) 


def expandFormula(biformulas, target, debug=False):
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
        

    


