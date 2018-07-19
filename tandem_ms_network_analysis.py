# Make sure you have the scpectrum_alignment module in your search path
from spectrum_alignment import score_alignment 
## We are going to use Networkx to manipulate and process networks##
import networkx as nx



###
# This function parses mgf files into two dictionaries
# One dictionary for the metadata, the other for the fragment ion m/z and intesity
# It takes the path to the mgf file as an input
# The limit argument is set to True if you want to process a specific number of spectra
###

def parseSpectra(path, NUM_OF_SPECTRA=10, limit=False):
    print "Running..."
    fileHandle = open(path, "rU")    # Store file handle object in variable
    
    ID = ""    # Defined ID variable (which will be our key)
    fragments = []    # Defined list (for storing tuples of m/z and intensity values)
    field_descriptions = []    # Defined list (for storing tuples of field descriptions, e.g (pepmass, 270))
    metadata = {}     # Defined dictionary (for storing metadata a.k.a field descriptions)
    spectra = {}      # Defined dictionary (for storing fragment ion info for one spectra)
    
    counter = 0    

    # Loop for handling each line in the mgf file
   
    for line in fileHandle:
        
        if limit == True:
            if counter > NUM_OF_SPECTRA:
                print "done"
                print "%d spectra parsed" % (NUM_OF_SPECTRA)
                break
            
        if line.startswith("BEGIN"):
            counter+=1
            
            #For each spectra, the list of tuples are re-initialized
            field_descriptions = [] 
            fragments =[]
            
        elif line.startswith("FEATURE"):
            feature_id = line.split("=")
            ID = feature_id[1][:-1]
            
            
        ## Creating metadata dictionary ##
        elif "=" in line:          # All the field descriptions have an "=" sign
            data = line.split("=")
            
            # The slicing notation is to remove the trailing "/n"
            field_tups = (data[0], data[1][:-1])    
            field_descriptions.append(field_tups)
            metadata[ID] = field_descriptions
            
            
        
        
        ## Creating mass spectra dictionary  ##
        elif not(line.isspace()) and "END" not in line:
            ion_data = line.split()
            ion_tups = (float(ion_data[0]), float(ion_data[1]))
            fragments.append(ion_tups)
            spectra[ID] = fragments
            
    if limit == False:
        print "%d spectra parsed" % (counter)"
    
    print "Finished running"
    return metadata, spectra
	

###	
# Function that scores the similarity between two spectra using the imported "score_alignment()" function
# It also writes the results to a file
# It takes the outputs of parseSpectra(), and a path for the result file as arguments
###

def process_spectra_similarity(metadata, spectral_data, path):
    fileout = open(path,"w")
    keys = metadata.keys()
    
    header = "The format is in ID, ID, similarity score" + "\n"
    fileout.write(header + "\n")
    
    # This lists holds tuples of indices
    # This list is used to prevent repeating comparisons
    reverse_order = []       
        
    for i,x in enumerate(keys):
        for j,y in enumerate(keys):
            
        # This condition prevents self comparison and repeating comparisons (happening in reverse order)
            if i != j and (i, j) not in reverse_order:             
                
                reverse_order.append((j,i))
        
        # The value is a list of tuples
        # Index 0 is the first tuple in the list (containing parent mass info)
        # Index 1 of the tuple is the parent mass value
                parent_mass1 = float(metadata[x][0][1])
                parent_mass2 = float(metadata[y][0][1]) 
                
                score, reported_alignments = score_alignment(spectral_data[x], 
                                                             spectral_data[y], parent_mass1, parent_mass2, 0.3)
                
                
                line = x + " " + y + " " + str(score) + "\n"
                fileout.write(line)
                
	
	###
	# This function is used to arrange the data in the file (produced from the process_spectra_similarity() function)  into an edge list 
	# There is an optioinal argument to specify the threshold for filtering edges based on their scores
	###
def generate_edgeList(path, threshold=0):
    fileHandle = open(path, "rU")
    Edge_list = []   # The edge list variable
    
    for line in fileHandle:
        
    # Note: info in line is in the format ID, ID, Score  
    # Score is the edge attribute, ID's are the start and end nodes (edges)
        data = line.split()
        edge_attribute = {}    # dictionary for storing edge attributes (scores)
        edge_attribute["score"] = data[2]
        
    # This condition is used to filter low scoring edges
        if float(data[2]) >= threshold:
            tups = (data[0], data[1], edge_attribute)    # This contains the edges and edge_attribute dictionary
            Edge_list.append(tups)
        
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
    ## Edges is not part of the top k edges of incident nodes are removed
    
    for node in graph_dictionary.keys():
        
        for edge in graph_dictionary[node]:
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