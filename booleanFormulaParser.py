

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
        else: 
            # print("Unary operator")
            self.val = self.right.val 
            self.left = self.right.left 
            self.right = self.right.right

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

def deleteSomeNodes(cur, nodenames, debug=False):
    if not cur:
        return 
    if not cur.val:
        return  
    if debug:
        cur.display()
    if cur.val in nodenames: # leaf node to be delete
        if debug:
            print("Leaf node {} is deleted".format(cur.val))
        parent = cur.parent 
        cur.val = None
        cur = None 
        # cur.delete_myself()
        if parent and debug:
            print("Now process")
            parent.display()
        deleteSomeNodes(parent, nodenames, debug)
    
    if cur:
        if (cur.left and cur.left.val) or (cur.right and cur.right.val): 
            if debug:
                print("There are left or right")
                cur.display()
            if cur.left and cur.left.val in nodenames:
                if debug:
                    print("Delete left node {} of".format(cur.left.val))
                    cur.display()
                cur.delete_left()
                if debug:
                    print("Now process")
                    cur.display()
                deleteSomeNodes(cur, nodenames)
            if cur.right and cur.right.val in nodenames:
                if debug:
                    print("Delete right node {} of".format(cur.right.val))
                    cur.display()
                cur.delete_right()
                if debug:
                    print("Now process")
                    cur.display()
                deleteSomeNodes(cur, nodenames)
            # if cur.left:
            deleteSomeNodes(cur.left, nodenames)
            # if cur.right:
            deleteSomeNodes(cur.right, nodenames)
        else:
            if debug:
                cur.display()
                print("Current node has no left or right")
            if cur.val and cur.val.upper() in ["AND", "OR"]:
                if debug:
                    print("Current node is a binary operator of both None parties")
                cur.val = None


# get result of a boolean formula given the value of inputs as dictionary 

def getResult(cur, inputs, debug=False):
    if not cur:
        return None
    if not cur.left and not cur.right: # leaf node 
        try:
            return inputs[cur.val]
        except:
            print("Cannot find {} in input list, return None", cur.val)
            return None
    else:
        if not cur.left: # not operator with no left 
            if cur.val.upper() != "NOT":
                cur.display()
            assert cur.val.upper() == "NOT", "Uncomplete binary operator"
            # print("Not {}".format(cur.right.val))
            return (not getResult(cur.right, inputs, debug))

        else: # may be nested, may be nuclear operator 
            if not cur.right:
                cur.display()
            assert cur.right, "Uncomplete binary operator"
            if cur.val.upper() == "AND":
                return (getResult(cur.left, inputs, debug) and getResult(cur.right, inputs, debug))
            elif cur.val.upper() == "OR":
                return (getResult(cur.left, inputs, debug) or getResult(cur.right, inputs, debug))
            else:
                assert False, "Non support operator"
            

            # if cur.left.left or cur.left.right: # there is something in the left 
            #     print("Go left")
            #     return getResult(cur.left, inputs, debug)
            
            # if cur.right.left or cur.right.right: # there is something in the right 
            #     print("Go right")
            #     return getResult(cur.right, inputs, debug)
            
            # if cur.val.upper() == "AND":
            #     try:
            #         print("{} and {} is {}".format(cur.left.val, cur.right.val, (inputs[cur.left.val] and inputs[cur.right.val])))
            #         return (inputs[cur.left.val] and inputs[cur.right.val])
            #     except:
            #         cur.display()
            #         print("Cannot find values in input list, return None")
            #         return None
                
            # if cur.val.upper() == "OR":
            #     try:
            #         print("{} or {} is {}".format(cur.left.val, cur.right.val, (inputs[cur.left.val] or inputs[cur.right.val])))
            #         return (inputs[cur.left.val] or inputs[cur.right.val])
            #     except:
            #         cur.display()
            #         print("Cannot find values in input list, return None")
            #         return None
            

    
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
    print('\n')
    print(formula)
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
    # root.display()
    return root 



