import networkx as nx
from pyvis.network import Network 
import copy 
import matplotlib.pyplot as plt 

class Node(object):
    def __init__(self, val):
        self.val = val
        self.left = None
        self.right = None
        self.parent = None 
        
    def insert_left(self, child):
        if self.left is None:
            self.left = child
            child.parent = self 
        else:
            child.left = self.left
            self.left = child
            child.parent = self 

    def insert_right(self, child):
        if self.right is None:
            self.right = child
            child.parent = self 
        else:
            child.right = self.right
            self.right = child
            child.parent = self

    def delete_left(self):
        self.val = self.right.val # take right up 
        self.left = self.right.left 
        self.right = self.right.right 
    def delete_right(self): 
        if self.left: # case of binary operator when the left is present 
            self.val = self.left.val # take the left up 
            self.right = self.left.right
            self.left = self.left.left 
    def delete_itself(self):
        self.val = None 
        self = None
        # if self.parent:
        #     if self.parent.left:
        #         self.val = self.parent.left.val 
        #         self.left = self.parent.left.left
        #         self.right = self.parent.left.right
        #     self.parent = self.parent.left
        # else:
        #     self.val = None 

    def display(self):
        lines, *_ = self._display_aux()
        for line in lines:
            print(line)

    def _display_aux(self):
        """Returns list of strings, width, height, and horizontal coordinate of the root."""
        # No child.
        if self.right is None and self.left is None:
            line = '%s' % self.val
            width = len(line)
            height = 1
            middle = width // 2
            return [line], width, height, middle

        # Only left child.
        if self.right is None:
            lines, n, p, x = self.left._display_aux()
            s = '%s' % self.val
            u = len(s)
            first_line = (x + 1) * ' ' + (n - x - 1) * '_' + s
            second_line = x * ' ' + '/' + (n - x - 1 + u) * ' '
            shifted_lines = [line + u * ' ' for line in lines]
            return [first_line, second_line] + shifted_lines, n + u, p + 2, n + u // 2

        # Only right child.
        if self.left is None:
            lines, n, p, x = self.right._display_aux()
            s = '%s' % self.val
            u = len(s)
            first_line = s + x * '_' + (n - x) * ' '
            second_line = (u + x) * ' ' + '\\' + (n - x - 1) * ' '
            shifted_lines = [u * ' ' + line for line in lines]
            return [first_line, second_line] + shifted_lines, n + u, p + 2, u // 2

        # Two children.
        left, n, p, x = self.left._display_aux()
        right, m, q, y = self.right._display_aux()
        s = '%s' % self.val
        u = len(s)
        first_line = (x + 1) * ' ' + (n - x - 1) * '_' + s + y * '_' + (m - y) * ' '
        second_line = x * ' ' + '/' + (n - x - 1 + u + y) * ' ' + '\\' + (m - y - 1) * ' '
        if p < q:
            left += [n * ' '] * (q - p)
        elif q < p:
            right += [m * ' '] * (p - q)
        zipped_lines = zip(left, right)
        lines = [first_line, second_line] + [a + u * ' ' + b for a, b in zipped_lines]
        return lines, n + m + u, max(p, q) + 2, n + u // 2

def deleteNode(cur, node, debug=False):
    if not cur:
        return
    if not cur.val: 
        return 
    if debug:
        print ("Try to delete {} from".format(node))
        cur.display()
        if debug:
            if cur.parent:
                print("With parent is {}".format(cur.parent.val))

    if cur.left and cur.right: # binary operator 
        # print(1)
        if debug: 
            print("Binary operator")
        if cur.left.val.upper() == "NOT" and cur.left.right.val == node:
            # print(2)
            cur.delete_left() 
            if debug:
                print("Just delete some not on the left, now work with")
            deleteNode(cur, node, debug)

        elif cur.right.val.upper() == "NOT" and cur.right.right.val == node:
            # print(3)
            cur.delete_right()
            if debug:
                print("Just delete some not on the right, now work with")
            deleteNode(cur, node, debug)
 
        else:
            if cur.left and cur.left.val == node:
                # print(4)
                if debug:
                    print("Delete node {} on the left".format(cur.left.val))
                cur.delete_left()
                deleteNode(cur, node, debug)
            if cur.right and cur.right.val == node:
                # print(5)
                if debug:
                    print("Delete node {} on the right".format(cur.right.val))
                cur.delete_right()
                deleteNode(cur, node, debug)
            else:
                if debug:
                    print("Now go right")
                deleteNode(cur.right, node, debug)
                if debug:
                    print("Now go left")
                deleteNode(cur.left, node, debug)
    elif cur.right: # unary operator but not to delete 
        # print(10)
        if debug:
            print("Unary operator")
            cur.display()
        if cur.right.val == node:
            cur.val = None
            cur.right = None
            cur = None 
    else:
        # print(11)
        if cur.val == node:
            if debug:
                print("Leaf node {} is to delete".format(cur.val))
            parent = cur.parent 
            # cur.parent.display()
            cur.val == None 
            # cur = None 
            cur.delete_itself()

            deleteNode(parent, node, debug)


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
        
def convertBooleanFormulas2Network(formulas, inputnames, speciesnames, filename, debug=False):
    print("Input are {}".format(inputnames))
    net = nx.DiGraph()

    for spe in speciesnames:
        if spe in inputnames:
            # this node is yellow
            net.add_node(spe, labels = spe, color='#FFFF00') 
        else:
            net.add_node(spe, labels = spe, color='#0000FF')
    
    # now scan formulas to add edges 
    for left, root in formulas.items():
        stack = [root]
        while stack:
            cur = stack.pop()
            if cur.left or cur.right:
                # do left first 
                if cur.left:
                    stack.append(cur.left)
                if cur.right:
                    stack.append(cur.right)
            else: # species at leaf node, create an edge 
                # trace back to root to find number of NOT operator to set edge color
                numnot = 0 
                numand = 0
                leaf = cur.val 
                if leaf not in speciesnames:
                    print("Finding leaf node not in species list {}".format(leaf))
                while cur.parent:
                    if cur.parent.val.upper() == "NOT":
                        numnot += 1
                    if cur.parent.val.upper() == "AND":
                        numand += 1
                    cur = cur.parent 
                if numnot%2 == 0:
                    if numand <= 0: # gray edge 
                        net.add_edge(leaf, left, color="#808080")
                    else:
                        net.add_edge(leaf, left, color="008000") 
                else:
                    net.add_edge(leaf, left, color="#FF0000") # not edges have highest weight 

    nt = Network(directed=True, height="100%", width="100%")
    nt.toggle_physics(False)
    
    arrangeLayers(net, inputnames, formulas, debug)
    nt.from_nx(net)

    nt.show(filename + '.html', notebook=False)

    return net 

def limitGraphAfterNode(net, inputnames, outputname, formulas, filename = 'limitednet', delcycle = False, debug=False):
    # firstly cut all the edges from outputnode, 
    # this may cut a cycle too if the node influences itself 
    try:
        outedges = set(net.out_edges(outputname))
        print(outedges)
    except:
        print("There is no node named {}".format(outputname))             
        # print(edges)

    # this is to test 
    if debug: 
        firstcpt = nx.simple_cycles(net) 
        print("Number of cycles before delete out going edges from output is {}".format(len(list(firstcpt))))

    net.remove_edges_from(outedges) 

    # this is also to test 
    if debug:
        secondcpt = nx.simple_cycles(net) 
        print("Number of cycles after delete out going edges from output is {}".format(len(list(secondcpt))))

    
    # secondly remove all the edges is not related with the output node at all 
    # this also remove cycles
    relatededges = set()
    for inputname in inputnames:  
        allpaths = list(nx.all_simple_edge_paths(net, inputname, outputname))
        thesepaths = set()
        for onepath in allpaths: 
            thesepaths = thesepaths.union(set(onepath))
            # print(set(onepath))
        # print(list(nx.all_simple_edge_paths(net, inputname, outputname)))
        relatededges = relatededges.union(thesepaths)

    alledges = set(net.edges)

    nonrelatededges = alledges.difference(relatededges) 
    print(nonrelatededges)


    print("None related edges to delete:")
    print(nonrelatededges)

    net.remove_edges_from(nonrelatededges)
    
    nt = Network(directed=True, height="100%", width="100%")
    nt.toggle_physics(False)
    
    arrangeLayers(net, inputnames, formulas, debug)
    nt.from_nx(net)

    nt.show(filename + '.html', notebook=False)

    # this is also for testing 
    # find all cycles in this graph 
    if debug:
        thirdcpt = nx.simple_cycles(net)
        print("Number of cycles after delete none related edges to output is {}".format(len(list(thirdcpt))))
    # print(cycles)
    # for cycle in cycles:
    #     print(cycle)

    
    # set of all edges to dele
    todeledges = outedges.union(nonrelatededges) 
    todeledges = nonrelatededges
    # todeledges = {('SOCS1', 'Jak1')}
    # now delete all the edges from the formulas 
    for todeledge in todeledges:
        if debug:
            print("Deleting edge {} from".format(todeledge))
            formulas[todeledge[1]].display()
        deleteNode(formulas[todeledge[1]], todeledge[0], False)
        if debug:
            print("After delete edge {}".format(todeledge))
            formulas[todeledge[1]].display()
            print("\n")
        if formulas[todeledge[1]].val == None:
            if debug:
                print("Node {} with None update function will be removed".format(todeledge[1]))
            del formulas[todeledge[1]] 
    
    # now remove cycles 
    removeCycles(net, outputname, debug)
    
# function to find to make the networks become acyclic 
def removeCycles(net, outputname, debug=False):
    alledges = net.edges(data=True)
    for edge in alledges:
        print(edge)


       


# return a dictionary of binary formulars 
# with keys are species (or temporary node) 
# and values are the corresponding biformulas 
def toBinaryFormulas(formulas, debug=False):
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
    tosort = []
    for term, formula in biformulas.items():
        form = dict()
        if '_XTR_' in term: # this is an extra node 
            count = formula.count('_XTR_') 
            coms = formula.split()
            numcoms = len(coms)
            if count == 0: # there is no extra node in the formulas 
                form['level'] = 0 # highest level 
            else: # there is at least 1 extra node in the formulas 
                if numcoms >= 3: # second highest level 
                    form['level'] = count 
                else:
                    form['level'] = 2
        else:
            form['level'] = 10 # least level of the original node 
        
        form['term'] = term 
        form['formula'] = formula
        tosort.append(form)

    
    # now short the list tosort according to level 
    toreturn = sorted(tosort, key=lambda d: d['level'])
    # print(toreturn)
    if debug:
        for form in toreturn:
            print(form)

    return toreturn 

def convertBiBooleanFormulas2Network(biformulas, inputnames, speciesnames, filename, debug=False):
    print("Input are {}".format(inputnames))
    print("Ordinary nodes are {}".format(speciesnames))
    net = nx.DiGraph() # Networkx 
    # for left, right in biformulas.items():
    for formula in biformulas:
        left = formula['term']
        right = formula['formula']
        if not net.has_node(left):
            if left in inputnames:
                # print("Input node {}".format(left))
                net.add_node(left, labels = left, color='#FFFF00') 
            elif left in speciesnames:
                # print("Original node {}".format(left))
                net.add_node(left, labels = left, color='#0000FF')
            else:
                # print("Added node {}".format(left))
                net.add_node(left, labels = left, color='#808080')
        coms = right.split()
        if len(coms) == 3: # binary 
            if not net.has_node(coms[0]):  
                if coms[0] in inputnames:
                    # print("Input node {}".format(coms[0]))
                    net.add_node(coms[0], labels = coms[0], color='#FFFF00') 
                elif coms[0] in speciesnames:
                    # print("Original node {}".format(coms[0]))
                    net.add_node(coms[0], labels = coms[0], color='#0000FF')
                else:
                    # print("Added node {}".format(coms[0]))
                    net.add_node(coms[0], labels = coms[0], color='#808080')
            if not net.has_node(coms[2]):
                if coms[2] in inputnames:
                    # print("Input node {}".format(coms[2]))
                    net.add_node(coms[2], labels = coms[2], color='#FFFF00') 
                elif coms[2] in speciesnames:
                    # print("Original node {}".format(coms[2]))
                    net.add_node(coms[2], labels = coms[2], color='#0000FF')
                else:
                    # print("Added node {}".format(coms[2]))
                    net.add_node(coms[2], labels = coms[2], color='#808080')
            if coms[1].lower() == "or":
                net.add_edge (coms[0], left, color='#808080', type='arrow', width=3) 
                net.add_edge (coms[2], left, color='#808080', type='arrow', width=3) 
            elif coms[1].lower() == "and":
                net.add_edge (coms[0], left, color='#008000', type='arrow', width=3) 
                net.add_edge (coms[2], left, color='#008000', type='arrow', width=3) 
        elif len(coms) == 2: # not operator (unary)
            if not net.has_node(coms[1]):
                if coms[1] in inputnames:
                    # print("Input node {}".format(coms[1]))
                    net.add_node(coms[1], labels = coms[1], color='#FFFF00') 
                elif coms[1] in speciesnames:
                    # print("Original node {}".format(coms[1]))
                    net.add_node(coms[1], labels = coms[1], color='#0000FF')
                else:
                    # print("Added node {}".format(coms[1]))
                    net.add_node(coms[1], labels = coms[1], color='#808080')
            if coms[0].lower() == "not":
                net.add_edge(coms[1], left, color="#FF0000")
        elif len(coms) == 1: # identical operator 
            if not net.has_node(coms[0]):
                if coms[0] in inputnames:
                    # print("Input node {}".format(coms[0]))
                    net.add_node(coms[0], labels = coms[0], color='#FFFF00') 
                elif coms[0] in speciesnames:
                    # print("Original node {}".format(coms[0]))
                    net.add_node(coms[0], labels = coms[0], color='#0000FF')
                else:
                    # print("Added node {}".format(coms[0]))
                    net.add_node(coms[0], labels = coms[0], color='#808080')
            net.add_edge(coms[0], left, color="#808080")
        else:
            print("Confusing binary operator {}".format(right))
    # nx.draw_circular(net)
    # plt.savefig('plotgraph.png', dpi=300, bbox_inches='tight')
    # plt.show()
    nt = Network(directed=True, height='100%', width='100%')
    nt.toggle_physics(False)
    
    arrangeLayers(net, inputnames, biformulas, debug)
    # nt.show_buttons(filter_=["physics"])
    nt.from_nx(net)
    
    nt.show(filename+ ".html", notebook=False)
    return net 

def arrangeLayers(net, inputnames, formulas, debug=False): 
    print("Arrange layout of the network for visualization")
    nodes = net.nodes()
    layersize = dict()
    for node in nodes:
        if node in inputnames:
            # set layer for this input node is 0
            # print(node)
            net.nodes[node]['layer'] = 0
        else:
        #     # set default layer -1 
            net.nodes[node]['layer'] = -1
    layersize[0] = len(inputnames)
    curs = copy.deepcopy(list(inputnames))
    traversed = set() 

    devoted = set()

    while curs:
        cur = curs.pop()
        # print ("Now working with node {}".format(cur))
        traversed.add(cur)
        # if net.nodes[cur]['layer'] != -1:
        #     traversed.add(cur)

        # get all children nodes of this cur node 
        outedges = list(net.out_edges(cur)) 

        if len(outedges) == 1:
            devoted.add(cur)

        for edge in outedges:
            child = edge[1]
            if child not in traversed and child not in curs:
                curs.append(child)
            # now set layer for child 
            if net.nodes[child]['layer'] == -1: # child doesnt have layer yet 
                # find other parent of child 
                # inedges = net.in_edges(child) # expect to have at most 2 inedges 
                # minlayer = 0
                # maxlayer = 0
                # for inedge in inedges:
                #     par = inedge[0]
                #     if net.nodes[par]['layer'] < minlayer:
                #         minlayer = net.nodes[par]['layer']
                #     if net.nodes[par]['layer'] > maxlayer:
                #         maxlayer = net.nodes[par]['layer']
                
                # only set layer for child node if its parents are set already 
                # if minlayer >= 0:
                # print("Set layer for node {} is {}".format(child, maxlayer + 1))
                net.nodes[child]['layer'] = net.nodes[cur]['layer'] + 1 
                if net.nodes[child]['layer'] not in layersize:
                    layersize[net.nodes[child]['layer']] = 1
                else:
                    layersize[net.nodes[child]['layer']] += 1
                
                # only add child to travers if child's layer is set 
                traversed.add(child)

    print("----Layer information----")
    print(len(layersize))
    print(layersize)

    print("----List of devoted nodes-----")
    print(devoted)

    # split the visualization space for each layer 
    assignednode = dict()
    for node in nodes:
        net.nodes[node]['y'] = (net.nodes[node]['layer'] + 1)*100
        if net.nodes[node]['layer'] not in assignednode:
            assignednode[net.nodes[node]['layer']] = 1
        else:
            assignednode[net.nodes[node]['layer']] += 1
        
        # now set y coordinate 
        y = ((assignednode[net.nodes[node]['layer']]%2)*2 - 1) * 100 * \
            (assignednode[net.nodes[node]['layer']]/2) + 1000 - \
            100 * ((net.nodes[node]['layer']%2)*2 - 1)
        
        net.nodes[node]['x'] = y
        
    if debug:
        nodeswithdata = net.nodes(data=True)
        for node in nodeswithdata:
            print(node)
    
    
            












        
            

