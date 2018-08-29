## We are going to use Networkx to manipulate and process networks##
import networkx as nx


### 
# This function is used to format the pairwise comparison file to an edge list.
# The edge list is used to make a networkx graph object.
# The input is a path to the file (produced by the process_spectra_similarity() function.
###
def generate_edgeList(path):
    fileHandle = open(path, "rU")
    Edge_list = []   # The edge list variable
    
    for line in fileHandle:
        
        if line[0].isdigit():
        
        # Note: info in line is in the format ID, ID, Score  
        #Score is the edge attribute, ID's are the start and end nodes (edges)
            data = line.split()
            edge_attribute = {}    # dictionary for storing edge attributes (scores)
            edge_attribute["score"] = data[2]
            tups = (data[0], data[1], edge_attribute)    # This contains the edges and edge_attribute dictionary
            Edge_list.append(tups)
        # This condition is used to filter low scoring edges
         #   if float(data[2]) >= threshold:
                
        
    return Edge_list
	
	
###
# This function filters the top k edges of each node.
# The ranking of each edge is determined by the score
# Note: Edges are removed ONLY if the lowest scoring edge of one node is not in the top k of the neighbouring node
# This means that some nodes may keep more than k edges
# Input is a networkx graph object, and the argument for k (which must be an integer)
###
def filter_topk (graph, k):
    
    # The graph dictionary stores nodes as keys and list of tuples as values
    # The each tuple consists of the incident node, and the score of the corresponding incident edge
    graph_dictionary = {}
    
    kept_edges = []  # List for storing edges that are kept
    discarded_edges = []   # List for storing edges that are to be removed
    unfiltered_edges = []   # List for storing all edges in a graph
    
    # Creating graph dictionary
    for node in graph.nodes():
        
        # Note: the incident node will always be the second element in the edge tuple
        incident_edges = graph.edges(node)
        node_score = []
        for edge in incident_edges:
            node_score.append((edge[1], graph.get_edge_data(edge[0], edge[1])["score"]))
        
        graph_dictionary.setdefault(node,node_score)
    
    
    # Filtering top k ..
    ## Each node (key) in the graph dictionary has its list of tuples sorted according to
    ## the score of the edge (DESCENDING order).
    ## Edges not part of the top k edges of both node and incident node are removed
    
    for node in graph_dictionary.keys():
        graph_dictionary[node].sort(key= lambda x:x[1], reverse=True)
        
        for edge in graph_dictionary[node][:k]:
           # other_node_edges = ggd[edge[0]]
            graph_dictionary[edge[0]].sort(key= lambda x:x[1], reverse=True)
            
            for other_edge in graph_dictionary[edge[0]][:k]:
                if edge[1] == other_edge[1]:
                    if {edge[0],node} not in kept_edges:
                        kept_edges.append({node,edge[0]})
                        
    
    # Removing the discarded edges from the graph...
    ## First make a list of edges in the graph (make sure the edges are sets NOT tuples)
    for edge in list(graph.edges()):
        unfiltered_edges.append(set(edge))
    
    ## Second, generate list of discarded (or unkept) edges
    for edge in unfiltered_edges:
        if edge not in kept_edges:
            discarded_edges.append(tuple(edge))
           
            
    
    #Remove the discarded edges from the graph
    graph.remove_edges_from(discarded_edges)
	
	
###
# This function ensures that connecteing components are not greater than a certain size threshold
# It keeps removing the lowest scoring edge in a connecting component until it size has broken down.
# Inputs are a networkx graph object and the size threshold (which is an integer)
###

def manage_component_size(graph, size):
    
    # The boolean that determines whether or not the components are within the size threshold
    fin = False
    

    while fin == False:
        to_be_processed = [] # List for storing subgraphs (created via nodes of connecting components)
        fin = True
        
        for component in nx.connected_components(graph):
            
            # This condition checks if any component is above the size threshold
            if len(component) > size:
                fin = False
            
            subG = graph.subgraph(list(component))
            to_be_processed.append(subG)
       
    
        for i,g in enumerate(to_be_processed):
                
        # Condition for removing lowest scoring edge from component    
            if len(g.nodes()) > size:
                ordered_edge_list = sorted(nx.get_edge_attributes(g, "score").items(),
                                           key= lambda edge:edge[1], reverse=True) 
                
               # node = ordered_edge_list[0][0][0]
                edge_to_remove = ordered_edge_list.pop()[0]
               

                graph.remove_edge(edge_to_remove[0],edge_to_remove[1])