class Node(object):
    def __init__(self, val):
        self.val = val
        self.left = None
        self.right = None
        self.parent = None 
        self.level = -1 # this serves only for the target function 
        self.oldVal = None 
        self.exVal = None 
    
        
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

    def display(self, showlevel=False, oldval=False): 
        lines, *_ = self._display_aux(showlevel, oldval)
        for line in lines:
            print(line)

    def _display_aux(self, showlevel=False, oldval=False):
        """Returns list of strings, width, height, and horizontal coordinate of the root."""
        # No child.
        if self.right is None and self.left is None:
            if showlevel:
                line = '%s' % self.val + "_" + str(self.level)
            else:
                line = '%s' % self.val
            if oldval:
                line = '%s' % self.oldVal
            # line = '%s' % self.val + "_" + str(self.level)
            width = len(line)
            height = 1
            middle = width // 2 
            return [line], width, height, middle

        # Only left child.
        if self.right is None:
            lines, n, p, x = self.left._display_aux(showlevel, oldval)
            if showlevel:
                s = '%s' % self.val + "_" + str(self.level)
            else:
                s = '%s' % self.val
            if oldval:
                s = '%s' % self.oldVal
            # s = '%s' % self.val + "_" + str(self.level)
            u = len(s)
            first_line = (x + 1) * ' ' + (n - x - 1) * '_' + s
            second_line = x * ' ' + '/' + (n - x - 1 + u) * ' '
            shifted_lines = [line + u * ' ' for line in lines]
            return [first_line, second_line] + shifted_lines, n + u, p + 2, n + u // 2

        # Only right child.
        if self.left is None:
            lines, n, p, x = self.right._display_aux(showlevel, oldval)
            # s = '%s' % self.val + "_" + str(self.level)
            if showlevel:
                s = '%s' % self.val + "_" + str(self.level)
            else:
                s = '%s' % self.val
            if oldval:
                s = '%s' % self.oldVal
            u = len(s)
            first_line = s + x * '_' + (n - x) * ' '
            second_line = (u + x) * ' ' + '\\' + (n - x - 1) * ' '
            shifted_lines = [u * ' ' + line for line in lines]
            return [first_line, second_line] + shifted_lines, n + u, p + 2, u // 2

        # Two children.
        left, n, p, x = self.left._display_aux(showlevel, oldval)
        right, m, q, y = self.right._display_aux(showlevel, oldval)
        # s = '%s' % self.val + "_" + str(self.level)
        if showlevel:
                s = '%s' % self.val + "_" + str(self.level)
        else:
            s = '%s' % self.val
        if oldval:
            s = '%s' % self.oldVal
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
    
def replaceNode(cur, node, newnode, debug=False):
    if not cur: 
        return
    if not cur.val:
        return
    if debug:
        print ("Try to replace {} with {} from".format(node, newnode))
        cur.display()
        if debug:
            if cur.parent:
                print("With parent is {}".format(cur.parent.val))
    if cur.val == node:
        cur.val = newnode 
    if cur.left:
        replaceNode(cur.left, node, newnode, debug)
    if cur.right:
        replaceNode(cur.right, node, newnode, debug)
    
        

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
