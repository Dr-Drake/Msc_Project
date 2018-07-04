# Make sure you have the scpectrum_alignment module in your search path
from spectrum_alignment import score_alignment 
## We are going to use Networkx to manipulate and process networks##
import networkx as nx


     # variable dictating the number of spectra processed in the mgf file

###
# This function parses mgf files into two dictionaries
# One dictionary for the metadata, the other for the fragment ion m/z and intesity
# It takes the path to the mgf file as an input
# It parses 10 spectra on default
###

def parseSpectra(path, NUM_OF_SPECTRA=10 ):
    fileHandle = open(path, "rU")    # Store file handle object in variable
    
    ID = ""    # Defined ID variable (which will be our key)
    fragments = []    # Defined list (for storing tuples of m/z and intensity values)
    field_descriptions = []    # Defined list (for storing tuples of field descriptions, e.g (pepmass, 270))
    metadata = {}     # Defined dictionary (for storing metadata a.k.a field descriptions)
    spectra = {}      # Defined dictionary (for storing fragment ion info for one spectra)
    
    counter = 0    

    # Loop for handling each line in the mgf file
   
    for line in fileHandle:
        
        if counter > NUM_OF_SPECTRA:
            print "done"
            break
            
        elif line.startswith("BEGIN"):
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
            
        
    return metadata, spectra
	

###	
# Function for scores the similarity between two spectra using the imported "score_alignment()" function
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
# This function pairs  the nodes in connected components with their edge attributes (scores) #
# It is a function designed for to help filter top n nodes 
###
def get_component_edge_attributes(graph):    # put graph as an argument
    component_dict = {}
    
    # This generator yields a set of nodes that make up a connected component
    for component in nx.connected_components(graph):
        scores = []
        
        for edge in graph.edges():
            
        # If there is an intersection between the set of nodes (component) and the set of edges
        # then the score of that edge is paired with component nodes (in a dictionary)
            if set(edge) & component != set([]) :
                
                scores.append(float(nx.get_edge_attributes(graph, "score")[edge]))
                component_dict.setdefault(tuple(component),scores)
            
    return component_dict 
	
	
####
# This function removes the lowest scoring edge from the of connected components present in the network (a.k.a graph)
# It returns an edge_list that has been top n filtered
###

def filter_top_n(graph):
    
    # Get the component : edge attribute pair of the graph
    component_edge_attributes = get_component_edge_attributes(graph)
    
    # Convert the graph to an edge list
    # Note that this edge list is now visualised as an edge data view rather than a list
    el = nx.to_edgelist(graph)
    
    # We convert the view to a list
    list_el = list(el)
    
    # List of the components (components are tuples of nodes)
    component_list = component_edge_attributes.keys()
    lowest_scores = []     # Stores the scores of the lowest scoring edge from each connected component
    
    for scores in component_edge_attributes.values():
        
        # Condition to filter edges in components with more than one edge
        if len(scores) > 1:
            scores.sort()
            scores.reverse()
            lowest_scores.append(scores.pop())
    
    for component in component_list:
        
        for i,edge in enumerate(list_el):
            edge_set = set([edge[0], edge[1]])     # This set contains the edge in the edge list (el)
            edge_score = float(edge[2].values()[0])  
            
            # Edge has to be present in component and has to have a score present in
            # lowest scores list, before it can be removed
            # It prevents edges with identical score from being removed accidentally
            if  edge_set & set(component) != set([]) and edge_score in lowest_scores:
                
                del list_el[i]  # We can only delete from a list, rather than an EdgeDataView
        
    
    return list_el
	
	