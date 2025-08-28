import networkx as nx
from pyvis.network import Network 
import random
import numpy as np
import copy 
from Node import deleteNode, replaceNode 
# convert formulas to network and visualize it    
def convertBooleanFormulas2Network(formulas, inputnames, speciesnames, filename, debug=False):
    print("Input are {}".format(inputnames))
    net = nx.DiGraph()

    for spe in speciesnames:
        if spe in inputnames:
            # this node is yellow
            net.add_node(spe, labels = spe, color='#FFFF00', size='15') 
        else:
            net.add_node(spe, labels = spe, color='#0000FF', size='15')
    
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

    # nt = Network(directed=True, height="100%", width="100%")
    # nt.toggle_physics(False)
    order = linear_arrangement(net)
    fas, fwas = getFAS(net, order, debug)
    if debug:
        print("Order is {}".format(order))
        print("FAS is {}".format(fas))
    layout(net, inputnames, fas, order, debug)
    # nt.from_nx(net)
    # nt.show(filename + '.html', notebook=False)

    return net, fas 


# def limitGraphAfterNode(net, inputnames, outputname, formulas, filename = 'limitednet', delcycle = False, debug=False):
#     # firstly cut all the edges from outputnode, 
#     # this may cut a cycle too if the node influences itself 
#     try:
#         outedges = set(net.out_edges(outputname))
#         print(outedges)
#     except:
#         print("There is no node named {}".format(outputname))             
#         # print(edges)

#     # this is to test 
#     if debug: 
#         firstcpt = nx.simple_cycles(net) 
#         print("Number of cycles before delete out going edges from output is {}".format(len(list(firstcpt))))

#     net.remove_edges_from(outedges) 

#     # this is also to test 
#     if debug:
#         secondcpt = nx.simple_cycles(net) 
#         print("Number of cycles after delete out going edges from output is {}".format(len(list(secondcpt))))

    
#     # secondly remove all the edges is not related with the output node at all 
#     # this also remove cycles
#     relatededges = set()
#     for inputname in inputnames:  
#         allpaths = list(nx.all_simple_edge_paths(net, inputname, outputname))
#         thesepaths = set()
#         for onepath in allpaths: 
#             thesepaths = thesepaths.union(set(onepath))
#             # print(set(onepath))
#         # print(list(nx.all_simple_edge_paths(net, inputname, outputname)))
#         relatededges = relatededges.union(thesepaths)

#     alledges = set(net.edges)

#     nonrelatededges = alledges.difference(relatededges) 
#     print(nonrelatededges)


#     print("None related edges to delete:")
#     print(nonrelatededges)

#     net.remove_edges_from(nonrelatededges)
    
#     nt = Network(directed=True, height="100%", width="100%")
#     nt.toggle_physics(False)
    
#     manipulateNetwork(net, inputnames, formulas, delcycle, debug)
#     nt.from_nx(net)

#     nt.show(filename + '.html', notebook=False)

#     # this is also for testing 
#     # find all cycles in this graph 
#     if debug:
#         thirdcpt = nx.simple_cycles(net)
#         print("Number of cycles after delete none related edges to output is {}".format(len(list(thirdcpt))))
#     # print(cycles)
#     # for cycle in cycles:
#     #     print(cycle)

    
#     # set of all edges to dele
#     todeledges = outedges.union(nonrelatededges) 
#     todeledges = nonrelatededges
#     # todeledges = {('SOCS1', 'Jak1')}
#     # now delete all the edges from the formulas 
#     for todeledge in todeledges:
#         if debug:
#             print("Deleting edge {} from".format(todeledge))
#             formulas[todeledge[1]].display()
#         deleteNode(formulas[todeledge[1]], todeledge[0], False)
#         if debug:
#             print("After delete edge {}".format(todeledge))
#             formulas[todeledge[1]].display()
#             print("\n")
#         if formulas[todeledge[1]].val == None:
#             if debug:
#                 print("Node {} with None update function will be removed".format(todeledge[1]))
#             del formulas[todeledge[1]] 
    
#     # now remove cycles 
#     # removeCycles(net, outputname, debug)


def signifiedCycles(net, allcycles, debug=False):
    # allarcs = net.edges() 
    for cycle in allcycles:
        for id, node in enumerate(cycle):
            if id < len(cycle): 
                arc = [node, cycle[id +1]]
                # now look for arc in the network 
                net.edges(node, cycle[id+1])


# net is a networkx digraph which has correct color for edges already 
def removeCycles(net, fas, possitiveOnly = True, debug = False):
    # acynet = copy.deepcopy(net) 
    if debug: 
        print("Size of FAS is {}".format(len(fas)))

    if possitiveOnly: # this condition is not in usage anymore, but keep the interface here just in case 
        # firstly, find all the cycles within the graph 
        allCycles = sorted(nx.simple_cycles(net))
        # print(allCycles)

        # secondly, separate possitive and negative cycles 
        for cycle in allCycles: 
            print(cycle) 

        # now we need to find only possitive loops 
        for fa in fas:
            print("Now examining edge {} ".format(fa))
            
            # here, find all cycles containing the feedback arc 

    else:
        net.remove_edges_from(fas) 


    # this is for test only 
    try:
        cycles = nx.find_cycle(net)
        print("Remaining cycles after deleteing feedback arcs: {}".format(cycles))
    except:
        print("No cycle found after deleting feedback arcs")

    # return acynet 

def bfs_arrangement(graph, input_nodes):
    """
    Assigns numbers to nodes in a directed graph using BFS, starting from input nodes.

    Parameters:
    - graph (networkx.DiGraph): A directed graph.
    - input_nodes (list): A list of nodes without incoming edges.

    Returns:
    - list: A list of nodes ordered by their BFS traversal (earlier visited nodes appear first).
    """
    # Initialize the queue with input nodes and a visited set
    queue = deque(input_nodes)
    visited = set(input_nodes)
    node_order = []  # List to store nodes in order of visitation

    while queue:
        node = queue.popleft()  # Get the next node in BFS order
        node_order.append(node)  # Append node to the ordered list

        # Iterate over outgoing neighbors (successors)
        for neighbor in graph.successors(node):
            if neighbor not in visited:
                visited.add(neighbor)
                queue.append(neighbor)

    return node_order

def linear_arrangement(G): # chatGPT code 
    """
    Implements the linear arrangement algorithm based on the provided pseudo-code.

    Parameters:
        G (nx.DiGraph): The input directed graph.

    Returns:
        list: A linear arrangement of the nodes.
    """
    # Copy the graph to avoid modifying the original
    graph = copy.deepcopy(G)

    s1 = []  # List to store the "s1" sequence
    s2 = []  # List to store the "s2" sequence

    while len(graph) > 0:
        # Process sinks
        while True:
            sinks = [node for node in graph.nodes if graph.out_degree(node) == 0]
            if not sinks:
                break
            for sink in sinks:
                s2.insert(0, sink)  # Prepend to s2
                graph.remove_node(sink)

        # Process sources
        while True:
            sources = [node for node in graph.nodes if graph.in_degree(node) == 0]
            if not sources:
                break
            for source in sources:
                s1.append(source)  # Append to s1
                graph.remove_node(source)

        # If the graph still has nodes, find the node with maximum δ(u) = out-degree(u) - in-degree(u)
        if len(graph) > 0:
            max_delta_node = max(graph.nodes, key=lambda u: graph.out_degree(u) - graph.in_degree(u))
            s1.append(max_delta_node)  # Append to s1
            graph.remove_node(max_delta_node)

    # Combine the two sequences
    return s1 + s2

def convertBiBooleanFormulas2Network(biformulas, inputnames, speciesnames, filename, acyclic = False, debug=False, extranodes=None):
    print("Input are {}".format(inputnames))
    print("Ordinary nodes are {}".format(speciesnames))
    print("Extranodes added after removing cycles are {}".format(extranodes))
    net = nx.DiGraph() # Networkx 
    # for left, right in biformulas.items():
    for formula in biformulas:
        left = formula['term']
        right = formula['formula']
        if not net.has_node(left):
            if left in inputnames:
                # print("Input node {}".format(left))
                net.add_node(left, labels = left, color='#FFFF00', size='15') 
            elif left in speciesnames:
                # print("Original node {}".format(left))
                net.add_node(left, labels = left, color='#0000FF',  size='15')   
            else:
                if extranodes:
                    if left in extranodes:
                        # print("PINK")
                        net.add_node(left, labels = left, color='#FFC0CB')
                # print("Added node {}".format(left))
                    else:
                        net.add_node(left, labels = left, color='#808080')
                else:
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
            net.add_edge(coms[0], left, color="#008000")
        else:
            print("Confusing binary operator {}".format(right))


    # nt = Network(directed=True, height='100%', width='100%')
    # nt.toggle_physics(False)
    if extranodes:
        # print("Extranodes is not empty")
        anet, aformular, extra, nodes_positions = manipulateNetwork(net, inputnames.union(set(extranodes)), biformulas, acyclic, True, debug)
    else:
        # print("Extranodes is empty")
        anet, aformular, extra, nodes_positions = manipulateNetwork(net, inputnames, biformulas, acyclic, True, debug)
    # nt.show_buttons(filter_=["physics"])
    # nt.from_nx(anet)
    
    # nt.show(filename+ ".html", notebook=False)
    # showNetwork(anet)
    return anet, nodes_positions 

def showNetwork(net, outname, koinshaps, kiinshap, koshaps, kishaps, filename='network.html', isbound=False):
    """
    Displays the network using pyvis and saves it to an HTML file.

    Parameters:
        net (networkx.DiGraph): The directed graph to display.
        filename (str): The name of the HTML file to save the visualization.
    """
    
    nt = Network(directed=True, height='100%', width='100%')
    nt.toggle_physics(False)
    print(koinshaps)
    if koinshaps:
        shaps = koinshaps
        newlables = dict()
        for node, shap in shaps.items():
            # print(node, shap)
            # net.get_node(node).update({'labels': str(node) + " ko: " + str(shap['pos'] + shap['neg']) +\
            #                            " ki: " + str(kiinshap[outname][node]['pos'] + kiinshap[outname][node]['neg'])})
            # newlables[node] = str(node) + "\n KO: " + str(shap['pos'] + shap['neg'])
            newlables[node] = str(node) + "\n KO: " + str(shap)
            
        if kiinshap:
            for node, shap in kiinshap.items():
                # newlables[node] += " | KI: " + str(kiinshap[node]['pos'] + kiinshap[node]['neg'])
                newlables[node] += " | KI: " + str(kiinshap[node]) 

        if koshaps: 
            # thiskoshaps = koshaps[outname]
            for node, shap in koshaps.items():
                # print(node, shap)
                # net.get_node(node).update({'labels': str(node) + " ko: " + str(shap['pos'] + shap['neg']) +\
                #                            " ki: " + str(kiinshap[outname][node]['pos'] + kiinshap[outname][node]['neg'])})
                if not isbound:
                    # newlables[node] = str(node) + "\n KO: " + str(shap['pos'] + shap['neg'])
                    newlables[node] = str(node) + "\n KO: " + str(shap)
                else:
                    # newlables[node] = str(node) + "\n KO: ( " + str(shap['pos']['lb'] + shap['neg']['ub']) + " to " + str(shap['pos']['ub'] + shap['neg']['lb']) + " )"
                    newlables[node] = str(node) + "\n KO: " + str(shap) 
        if kishaps:
            # print(kishaps)
            # thiskishaps = kishaps[outname]
            for node, shap in kishaps.items():
                # print("DMMMMMM")
                # print(node, shap)
                # net.get_node(node).update({'labels': str(node) + " ko: " + str(shap['pos'] + shap['neg']) +\
                #                            " ki: " + str(kiinshap[outname][node]['pos'] + kiinshap[outname][node]['neg'])})
                if not isbound:
                    # newlables[node] += " | KI: " + str(shap['pos'] + shap['neg'])
                    newlables[node] += " | KI: " + str(shap)
                else:
                    # newlables[node] += " | KI: " + "( " + str(shap['pos']['lb'] + shap['neg']['ub']) + " to " + str(shap['pos']['ub'] + shap['neg']['lb']) + " )"
                    newlables[node] += " | KI: " + str(shap)
        labelednet = nx.relabel_nodes(net, newlables, True)

    else:
        labelednet = net
        
    nt.from_nx(labelednet)
    nt.show(filename, notebook=False) 

def find_index(lst, element):
    """
    Finds the index of the first occurrence of an element in a list.

    Parameters:
        lst (list): The list to search in.
        element: The element to find.

    Returns: 
        int: The index of the element if found, or -1 if not found.
    """
    try:
        return lst.index(element)
    except ValueError:
        return -1  # Return -1 if the element is not in the list


def getFAS(net, order, debug=False):
    FAS = []
    FWAS = []
    edges = net.edges
    for (source, sink) in edges:
        assert find_index(order, source) >= 0, print("Cannot find order of node {}".format(source))
        assert find_index(order, sink) >= 0, print("Cannot find order of node {}".format(sink))
        if find_index(order, source) >= find_index(order, sink):
            FAS.append((source, sink))
        else:
            FWAS.append((source, sink))
    if debug:
        print("APPROXIMATION OF THE FEEDBACK ARC SET IS \n {}".format(FAS))
        print("AND THE REST \n {}".format(FWAS))
    return FAS, FWAS 


# if acyclic is on, remove all the feedback arc and modify also the formulas 
def manipulateNetwork(net, inputnames, formulas, acyclic = False, binary = False, debug=False): 
    # first get the linear arrangement 
    order = linear_arrangement(net) 
    # order = bfs_arrangement(net, inputnames)

    print("ARRANGEMENT IS \n {}".format(order))

    # now find the approximate minimum feedback arc set fas and the 
    fas, fwas = getFAS(net, order, debug)

    # copy a new network and set of formulas and modify this copies, 
    # keep the original network and formulas intact 
    anet = copy.deepcopy(net)
    aformulas = copy.deepcopy(formulas)
    extranodes = []

    if acyclic:
        '''
        # now modify the formulas accordingly 
        # 1st option: this is to simply remove the arc, other option is below 
        for todeledge in fas:
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
        '''
        # 2nd option to manipulate the formulas 
        for todeledge in fas:
            if debug:
                print("Deleting edge {} from".format(todeledge))
                aformulas[todeledge[1]].display()
            # just simply change to name of the original node to a new node with the name is the same with the deleted arc 
            newnodename = str(todeledge[0]) + "_to_" + str(todeledge[1])
            if debug:
                print("Replacing node {} with {}".format(todeledge[0], newnodename))
            replaceNode(aformulas[todeledge[1]], todeledge[0], newnodename, False)
            if debug:
                print("After replace {} with {}".format(todeledge[0], newnodename))
                aformulas[todeledge[1]].display()
                print("\n")
            # now add new node and arc to the networkx graph
            anet.add_node(newnodename, labels = newnodename, color='#FFC0CB', size='25')
            extranodes.append(newnodename)
            newedgecolor = net.get_edge_data(todeledge[0], todeledge[1], "color")
            if not newedgecolor:
                print("Cannot find color of the edge {}".format(todeledge)) 
            else:
                anet.add_edge(newnodename, todeledge[1], color = newedgecolor)
        if debug:       
            print("Extra added nodes in the effort to remove cycles are:")
            print(extranodes)
        # now remove the arc
        removeCycles(anet, fas, possitiveOnly=False, debug=debug) # this function work only with the networkx graph 
    # layout(anet, inputnames, fas, order, False)
    # layout(net, inputnames, fas, order, False)
    node_positions = None
    if acyclic or binary:
        aorder = linear_arrangement(anet)
        # aorder = bfs_arrangement(net, inputnames)
        node_positions = layout_acylic(anet, inputnames.union(set(extranodes)), fas, aorder, False)
    return anet, aformulas, extranodes, node_positions 


def assign_x_coordinates(graph, layer_spacing=100, node_spacing=10):
    """
    Assigns x-coordinates to nodes in a directed graph while minimizing edge lengths.
    Nodes are center-aligned within each layer.

    Parameters:
    - graph (networkx.DiGraph): A directed graph with nodes having a 'layer' attribute.
    - layer_spacing (int): Minimum vertical distance between layers.
    - node_spacing (int): Additional spacing between adjacent nodes.

    Returns:
    - dict: A dictionary mapping each node to its (x, y) position.
    """
    # Group nodes by layer
    layer_dict = {}
    for node, data in graph.nodes(data=True):
        layer = data.get("layer", 0)
        if layer not in layer_dict:
            layer_dict[layer] = []
        layer_dict[layer].append(node)

    # Sort layers by increasing order
    sorted_layers = sorted(layer_dict.keys())

    # Store node positions
    positions = {}

    # Track maximum width to center align nodes
    max_layer_width = 0

    # First pass: Determine x-coordinates and find max width
    layer_widths = {}
    for layer in sorted_layers:
        nodes_in_layer = layer_dict[layer]

        # Sort nodes based on predecessors' average x-coordinates
        nodes_in_layer.sort(key=lambda node: np.mean(
            [positions[pred][0] for pred in graph.predecessors(node) if pred in positions]
        ) if any(pred in positions for pred in graph.predecessors(node)) else 0)

        # Compute x-coordinates for this layer
        x_positions = []
        x_position = 0
        for node in nodes_in_layer:
            x_positions.append(x_position)
            x_position += layer_spacing + node_spacing

        # Store width and update max width
        layer_widths[layer] = x_positions[-1] if x_positions else 0
        max_layer_width = max(max_layer_width, layer_widths[layer])

    # Second pass: Assign positions with center alignment
    for layer in sorted_layers:
        nodes_in_layer = layer_dict[layer]
        layer_width = layer_widths[layer]

        # Center align nodes within max width
        x_offset = (max_layer_width - layer_width) / 2  # Centering offset

        for i, node in enumerate(nodes_in_layer):
            positions[node] = (x_offset + (layer_spacing + node_spacing) * i, layer * layer_spacing)

    return positions
 
def layout_acylic(net, inputnames, fas, order, debug=False):
    print("Now arrange layout of the original network for visualization")
    # print("Starting from input set: {}".format(inputnames))
    nodes = net.nodes()

    node_layers = {}

    # Dictionary to store the number of nodes in each layer
    layersize = {}

    # List to store devoted nodes (nodes with exactly one outgoing edge)
    devoted = []

    # Perform a topological sort to process nodes in dependency order
    order = list(nx.topological_sort(net))
    # print("\nOrder is {}".format(order))
    for node in order:
        # Assign layer as max(layer of predecessors) + 1
        if net.in_degree(node) == 0:
            node_layers[node] = 0  # Input nodes (no dependencies) go to layer 0
        else:
            node_layers[node] = max(node_layers[pred] for pred in net.predecessors(node)) + 1

        # Count nodes per layer
        layer = node_layers[node]
        layersize[layer] = layersize.get(layer, 0) + 1

        # Check if node has exactly one outgoing edge
        if net.out_degree(node) == 1:
            devoted.append(node)
    
    for node in nodes:
        net.nodes[node]['layer'] = node_layers[node]
    
    print("----Layer information----") 
    print(node_layers)
    
    # split the visualization space for each layer 
    assignednode = dict()

    # save the positions of nodes 
    nodes_positions = {node: (0, 0) for node in nodes}

    for node in nodes:
        net.nodes[node]['y'] = (net.nodes[node]['layer'] + 1)*100 + random.randint(-30, 30)
        nodes_positions[node] = (0, net.nodes[node]['y'])  # Initialize x to 0, y to layer position
        if net.nodes[node]['layer'] not in assignednode:
            assignednode[net.nodes[node]['layer']] = 1
        else:
            assignednode[net.nodes[node]['layer']] += 1
        
        # now set y coordinate 
        # y = ((assignednode[net.nodes[node]['layer']]%2)*2 - 1) * 250 * \
        #     (assignednode[net.nodes[node]['layer']]/2) + 1000 - \
        #     250 * ((net.nodes[node]['layer']%2)*2 - 1) + random.randint(-30, 30) 
        
        # net.nodes[node]['x'] = y
    
    # now assign the x coordinate of the node, iterate layer by layer 
    x_positions = assign_x_coordinates(net)
    # print("\n X positions are {}".format(x_positions))

    for node in nodes:
        net.nodes[node]['x'] = x_positions[node][0] + random.randint(-20, 20)
        nodes_positions[node] = (net.nodes[node]['x'], net.nodes[node]['y'])  # Update x position
        # net.nodes[node]['y'] = x_positions[node][1]


    if debug:
        nodeswithdata = net.nodes(data=True)
        for node in nodeswithdata:
            print(node)
            
    return nodes_positions 

    

def count_paths_to_target(graph, target):
    """
    Computes the number of paths from each node to the target node in a Directed Acyclic Graph (DAG).

    Parameters:
    - graph (networkx.DiGraph): A directed acyclic graph.
    - target (node): The target node.

    Returns:
    - dict: A dictionary where keys are nodes and values are the number of paths from that node to the target.
    """
    # Ensure the graph is a Directed Acyclic Graph (DAG)
    if not nx.is_directed_acyclic_graph(graph):
        raise ValueError("The input graph must be a Directed Acyclic Graph (DAG).")

    # Get topological ordering of nodes (ensures valid dependency processing)
    topo_order = list(nx.topological_sort(graph))

    # Dictionary to store the number of paths to the target
    path_count = {node: 0 for node in graph.nodes}
    
    # The target node has exactly one path to itself
    path_count[target] = 1  

    # Process nodes in reverse topological order (from target to source)
    for node in reversed(topo_order):
        if node == target:
            continue  # The target node is already set to 1
        # Sum the paths from successors
        path_count[node] = sum(path_count[succ] for succ in graph.successors(node))

    return path_count

def layout(net, inputnames, fas, order, debug=False):
    print("Now arrange layout of the original network for visualization")
    # print("Starting from input set: {}".format(inputnames))
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

    while curs: # curs is a list 
        cur = curs.pop() 
        # print ("Now working with node {}".format(cur))
        traversed.add(cur)
        # if net.nodes[cur]['layer'] != -1:
        #     traversed.add(cur)

        # get all children nodes of this cur node 
        outedges = list(net.out_edges(cur)) 

        if len(outedges) == 1:
            devoted.add(cur)
        
        tem = [] # stores the indexes of sink nodes from current node 
        for edge in outedges:
            if edge in fas:
                print("Encounter feedback arc {}".format(edge))
                continue
            else:
                tem.append(find_index(order, edge[1])) # get the sink vertex 
        tem = sorted(tem)
        # print(tem)
        for sinkid in tem:
            sink = order[sinkid]
            if sink not in traversed and sink not in curs:
                curs.append(sink)
            if net.nodes[sink]['layer'] == -1:
                net.nodes[sink]['layer'] = net.nodes[cur]['layer'] + 1 
                if net.nodes[sink]['layer'] not in layersize:
                    layersize[net.nodes[sink]['layer']] = 1
                else:
                    layersize[net.nodes[sink]['layer']] += 1
                traversed.add(sink)

    # after all, the remaining nodes without layer are the one only be the sink in feedback arcs 
    # adding an extra layer for this node 
    nextlayer = len(layersize)  
    layersize[nextlayer] = 0      
    for node in list(net.nodes):
        if net.nodes[node]['layer'] == -1:
            print("Couldn't assign layer to node {}".format(node)) 
            net.nodes[node]['layer'] = nextlayer
            layersize[nextlayer] += 1
    
    # # layerindex = dict() 
    # # for node in list(net.nodes):
    # #     if net.nodes[node]['layer'] not in layerindex:
    # #         layerindex[net.nodes[node]['layer']] = [node]

    # #     else: 
    # #         layerindex[net.nodes[node]['layer']].append(node)
    
    # # for i in sorted(layerindex.keys()):
    # #     print("Layer {} includes:{}".format(i,layerindex[i]))


    print("----Layer information----")
    print(len(layersize))
    print(layersize)

    print("----List of devoted nodes-----")
    print(devoted)

    # split the visualization space for each layer 
    assignednode = dict()
    for node in nodes:
        net.nodes[node]['y'] = (net.nodes[node]['layer'] + 1)*150 + random.randint(-30, 30)
        if net.nodes[node]['layer'] not in assignednode:
            assignednode[net.nodes[node]['layer']] = 1
        else:
            assignednode[net.nodes[node]['layer']] += 1
        
        # now set y coordinate 
        # y = ((assignednode[net.nodes[node]['layer']]%2)*2 - 1) * 250 * \
        #     (assignednode[net.nodes[node]['layer']]/2) + 1000 - \
        #     250 * ((net.nodes[node]['layer']%2)*2 - 1) + random.randint(-30, 30) 
        
        # net.nodes[node]['x'] = y
    
    # now assign the x coordinate of the node, iterate layer by layer 
    x_positions = assign_x_coordinates(net)

    for node in nodes:
        net.nodes[node]['x'] = x_positions[node]


    if debug:
        nodeswithdata = net.nodes(data=True)
        for node in nodeswithdata:
            print(node)
