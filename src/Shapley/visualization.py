import networkx as nx
from pyvis.network import Network 

def showNetwork(net, koinshaps, kiinshap, koshaps, kishaps, filename='network.html'):
    """
    Displays the network using pyvis and saves it to an HTML file.

    Parameters:
        net (networkx.DiGraph): The directed graph to display.
        filename (str): The name of the HTML file to save the visualization.
    """
    
    nt = Network(directed=True, height='100%', width='100%')
    nt.toggle_physics(False)
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
                newlables[node] = str(node) + "\n KO: " + str(shap) 
        if kishaps:
            # print(kishaps)
            # thiskishaps = kishaps[outname]
            for node, shap in kishaps.items():
                # print("DMMMMMM")
                # print(node, shap)
                # net.get_node(node).update({'labels': str(node) + " ko: " + str(shap['pos'] + shap['neg']) +\
                #                            " ki: " + str(kiinshap[outname][node]['pos'] + kiinshap[outname][node]['neg'])})
                newlables[node] += " | KI: " + str(shap)
        labelednet = nx.relabel_nodes(net, newlables, True)

    else:
        labelednet = net
        
    nt.from_nx(labelednet)
    nt.show(filename, notebook=False) 