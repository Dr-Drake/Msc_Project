from collections import namedtuple

# To caluclate the similarity (cosine) scores we make use of these two modules
# Make sure you add them to your search path
import scoring_functions as sf
import spectrum_alignment as sa

## We are going to use Networkx to manipulate and process networks##
import networkx as nx

# Variable (attribute) for storing spectral library paths
LIBRARIES = []

###
# This function adds paths to a libraries to the LIBRARIES list.
# It can take multiple paths as an input
###

def set_libraries(*library_paths):
    del LIBRARIES[:]
    for path in library_paths:
        LIBRARIES.append(path)
		

###
# This function parses mgf files into two dictionaries
# One dictionary for the metadata, the other for the fragment ion m/z and intesity
# It takes the path to the mgf file as an input
# The limit argument is set to True if you want to process a specific number of spectra
###

def parseSpectra(path, NUM_OF_SPECTRA=10, limit=False):
    print "Parsing..."
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
            ID = ""
            
        elif line.startswith("FEATURE") or line.startswith("SPECTRUMID"):
            feature_id = line.split("=")
            ID = feature_id[1][:-1]
            
            if len(field_descriptions) > 0:
                metadata[ID] = field_descriptions
            
        ## Creating metadata dictionary ##
        elif "=" in line:
            if not(line.startswith("FEATURE")) or not(line.startswith("SPECTRUMID")):
            
                data = line.split("=")
            
            # The slicing notation is to remove the trailing "/n"
                field_tups = (data[0], data[1][:-1])    
                field_descriptions.append(field_tups)
                
                if ID != "":
                    metadata[ID] = field_descriptions
            
            
        
        
        ## Creating mass spectra dictionary  ##
        elif line[0].isdigit():
            ion_data = line.split()
            ion_tups = (float(ion_data[0]), float(ion_data[1]))
            fragments.append(ion_tups)
            spectra[ID] = fragments
            
    
    print "%d spectra parsed" % (counter)
    
    print "Finished parsing"
    return metadata, spectra
	

###	
# This function performs a pairwise calculation of similarity scores. 
# It takes metadata and spectral data libraries generated from the parseSpectra() function.
# The tolerance for matching peaks, when calculating the similarity score, has default value of 0.3.
# The minimum number of matched peaks between spectra, when calculating the similarity score, has a default
# value of 6.
# The parent mass difference for indicating a possible identical match (or consensus) has default value of 0.2 Da.
# The output is a report of the ID's and their cosine score. The output path should be specified.
###

def process_spectra_similarity(metadata, spectral_data, path, tolerance=0.3, min_match=6, pm_range=0.2):
    
    # It's converted to an integer to accommodate the range() function
    k = pm_range * 1000 
    
    fileout = open(path,"w")
    keys = metadata.keys()
    
    # Named tuple for storing peak information necessary for calculating spectra similarity
    spec = namedtuple("spec",["n_peaks", "normalised_peaks", "parent_mz"])
    
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
                
                # Condition to prevent potential comparison of the same compound
            # Here, it is believed that compounds with parent massess within 0.2 Da of each other are
            # the same compound
                if parent_mass1 * 1000 not in range(int((parent_mass2 * 1000)) - int(k), 
                                                    int((parent_mass2 * 1000)) + int(k) + 1): 
                    
                    # Normalize spectra
                    spec1_n = sa.sqrt_normalize_spectrum(sa.convert_to_peaks(spectral_data[x]))
                    spec2_n = sa.sqrt_normalize_spectrum(sa.convert_to_peaks(spectral_data[y]))
                
                    # Formatting tuple for spec1 and spec2
                    nt_spec1 = (len(spectral_data[x]), spec1_n, parent_mass1)
                    nt_spec2 = (len(spectral_data[y]), spec2_n, parent_mass2)
                
                        # Calculating cosine score
                    score, used_matches = sf.fast_cosine_shift(spec(*nt_spec1), spec(*nt_spec2), 
                                                               tolerance, min_match)
                
                #score, reported_alignments = score_alignment(spectral_data[x], 
                                                            # spectral_data[y], parent_mass1, parent_mass2, 0.3)
                
                line = x + " " + y + " " + str(score) + "\n"
                fileout.write(line)
                
	
###
# This function returns the name, and the pubmed ID, of a compound.
# It is intended for thr find_library_match function, returning these values for spectra found to match
# with the query spectra dataset.
# It takes a list of metadata as input
###

def identify_match_compound(metadata_list):
    
    compound_name = ""
    pubmed_id = ""
    
    for feature in metadata_list:
        if feature[0] == "NAME":
            compound_name = feature[1]
        
        if feature[0] == "PUBMED": 
            pubmed_id = feature[1]
    
    return compound_name, pubmed_id
	


###
# This function finds spectral library hits. The library searched are the ones present in the LIBRARIES list.
# It takes query metadata and spectral data libraries generated from the parseSpectra() function.
# The tolerance for matching peaks, when calculating the similarity score, has default value of 0.3
# The minimum number of matched peaks between spectra, when calculating the similarity score, has a default
# value of 6.
# The parent mass difference for indicating a possible identical match (or consensus) has default value of
# 0.2 Da
# The output is a report of the matches found. The output path should be specified.
# If no matches are found, then the report will only have the header.
###


def find_library_match(query_metadata, query_spec, path, tolerance=0.3, min_match=6, pm_range=0.2, cosine_threshold=0.95):
    
    print "Looking for hits..."
    
    # It's converted to an integer to accommodate the range() function
    k = pm_range * 1000 
    
    # Named tuple for storing peak information necessary for calculating spectra similarity
    spec = namedtuple("spec",["n_peaks", "normalised_peaks", "parent_mz"])
    
    fileout = open(path,"w")
    header = "ID" + "\t" + "Identification status" + "\t" + "Compound name" + "\t" + "Match ID" + "\t" + "Pubmed ID" + "\t" + "Match score" + "\n"
    fileout.write(header)
    
    # Each spectral library is parsed and each peak info is iterated over until a match is found
    for library in LIBRARIES:
            l_mdata, l_sdata = parseSpectra(library)
           # print l_mdata.keys()
            
            for x in query_metadata.keys():
                q_parent_mass = float(query_metadata[x][0][1])
                
                for l_keys in l_mdata.keys():
                    #print l_keys
                    l_parent_mass = float(l_mdata[l_keys][0][1])
                    
                    
                    try:    
                    # Similarity scores are only calculated if the parent masses are within a certain
                    # range of each other
                
                         if l_parent_mass * 1000 in range(int((q_parent_mass * 1000)) - int(k),
                                                     int((q_parent_mass * 1000)) + int(k) + 1):
                            
                            # Normalize spectra
                            l_spec = sa.sqrt_normalize_spectrum(sa.convert_to_peaks(l_sdata[l_keys]))
                            q_spec = sa.sqrt_normalize_spectrum(sa.convert_to_peaks(query_spec[x]))
                
                            # Formatting tuple for spec1 and spec2
                            t_Lspec = (len(l_sdata[l_keys]), l_spec, l_parent_mass)
                            t_Qspec = (len(query_spec[x]), q_spec, q_parent_mass)
                
                            # Calculating cosine score
                            score, used_matches = sf.fast_cosine_shift(spec(*t_Lspec), spec(*t_Qspec), 
                                                               tolerance, min_match)
                        
                            # It's only a match if the cosine score above a 0.95 threshold
                            # A report of the match is given in a file
                            if score >= cosine_threshold:
                                compound_name, pubmed_id = identify_match_compound(l_mdata[l_keys])
                            
                                line = x + "\t" + "identified" + "\t" + compound_name + "\t" + l_keys + "\t" + pubmed_id + "\t" + str(score) + "\n"
                                fileout.write(line)
                    except KeyError:
                        print "%s does not have peak data" % (l_keys)
            
    print "Finished looking for hits"



	
### 
# This function is used to format the pairwise comparison file to an edge list.
# The edge list is used to make a networkx graph object.
# The input is a path to the file (produced by the process_spectra_similarity() function.
# There is an optioinal argument to specify the threshold for filtering edges based on their scores. The default is 0.5
###
def generate_edgeList(path, threshold=0.5):
    fileHandle = open(path, "rU")
    Edge_list = []   # The edge list variable
    
    for line in fileHandle:
        
        if line[0].isdigit():
        
        # Note: info in line is in the format ID, ID, Score  
        #Score is the edge attribute, ID's are the start and end nodes (edges)
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
               