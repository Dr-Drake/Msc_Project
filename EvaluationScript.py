# Neccesary packages
import networkx as nx
from ms2network import ms_profile as mp
from ms2network import graph_tools as gt


#Input files
mzml_query1 = "Beer_multibeers_Pooled_T10_POS.mzML"
mzml_query2 = "Beer_multibeers_2_T10_POS.mzML"
mzml_query3 = "Beer_multibeers_3_T10_POS.mzML"

# Libraries
GNPS = "C:\\Users\\2060576e\\Documents\\ALL_GNPS.mgf"
HMDB ="C:\\Users\\2060576e\\Documents\\HMDB.mgf"
MASSBANK = "C:\\Users\\2060576e\\Documents\\MASSBANK.mgf"

# Output files
net_cluster = "computed_network.txt"
Lib_annot = "Library_annotation.txt"
Annot = "additional_annotations.txt"
network = "computed_network.graphml"


# Create dataset
dataset = mp.MSdata()

# Load dataset
dataset.add_file(mzml_query1, "mzml", "extract3")
dataset.parse(limit=True, NUM_OF_SPECTRA=500)
dataset.add_file(mzml_query2, "mzml", "extract1")
dataset.parse(limit=True, NUM_OF_SPECTRA=500)
dataset.add_file(mzml_query3, "mzml", "extract2")
dataset.parse(limit=True, NUM_OF_SPECTRA=500)

# Number of specta in dataset
print "Number of MS/MS spectra in the dataset is",
print len(dataset)

# Build consensus
dataset.build_consensus(pm_range=0.5,score_threshold=0.7)


# Some clustering statistics
print "The number of clusters (or consensuses/unique spectra) is",
print len(dataset.consensus_data)
print "The number of clusters with more than one spectra is",
print len(dataset.clustered)
print "The number of clusters with only one spectra is",
print len(dataset.unclustered)


# Load libraries
dataset.add_library(HMDB, "mgf", "HMDB")
dataset.parse_library()
dataset.add_library(GNPS, "mgf", "GNPS")
dataset.parse_library()
dataset.add_library(MASSBANK, "mgf", "MASSBANK")
dataset.parse_library()

# Total number of spectra in the libraries
print "The number of MS/MS spectra parsed from the libraries is",
print len(dataset.library_spec_data)

# Generate annotation files
dataset.find_library_match(Lib_annot,score_threshold=0.7,pm_range=0.5)
dataset.spectrum_annotation_file(Annot)


# Compute molecular network
dataset.compute_molecular_network(net_cluster,edge_filer=0.5)


## Simplifiy molecular network ##
# Convert to networkx graph object
e_list = gt.generate_edgeList(net_cluster)
Graph = nx.from_edgelist(e_list)

# Process the graph object (Filter edges)
gt.filter_topk(Graph, 10)
gt.manage_component_size(Graph, 100)

# Convert graph object to graphml file
nx.write_graphml(Graph, network)




