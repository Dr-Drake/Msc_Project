from spec import Spectrum
import consensus

# To caluclate the similarity (cosine) scores we make use this modules
import scoring_functions as sf

# To parse files in th e mzml format, we use the pymzml package 
import pymzml

# Matplot lib and numpy packages are needed for ploting spectrum
import matplotlib.pyplot as plt
import numpy as np


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



##################################################################################################################

# The name of the script will be called ms_profile
class MSdata(object):
    def __init__(self): #group, limit=False, NUM_OF_SPECTRA=10):
        self.files = []
        self.file_formats = ["mzml","mgf"]
        self.spectral_data = []
        self.index = 0
        self.parsed = []
        self.consensus_data = []
        self.clustered = []
        self.unclustered = []
        self.libraries = []
        self.parsed_libraries = []
        self.library_spec_data = []
    
    
    def add_library(self, library_paths, file_format, name):
        if (library_paths,file_format,name) in self.libraries:
            raise Exception("You have already added a library with the exact same specifications")
        else:
            self.libraries.append((library_paths,file_format,name))
    
    
    def load_consensus(self, consensus):
        self.consensus_data = consensus
        
    def assess_top_peaks(self, spec1, spec2, tolerance=0.3):
        spec1_top = spec1.top_peaks
        spec2_top = spec2.top_peaks
        
        #answer = None
        
        for peak1 in spec1_top:
            for peak2 in spec2_top:
                if abs(peak1[0]-peak2[0]) <= tolerance:
                    answer = True
                    break
                else:
                    answer = False
                    
        return answer
    
    
    def reset_library(self):
        self.libraries = []
        self.parsed_libraries = []
        self.library_spec_data = []
    
    def add_file(self, path, file_format, group):
        if (path, file_format, group) in self.files:
            raise Exception("You have already added a file with the exact same specifications")
        else:
            self.files.append((path, file_format, group))
    
    
        
        
    def read_mzml(self, file, limit, NUM_OF_SPECTRA, tolerance, library=False):
        metadata = []     # List of tuples for storing metadata a.k.a field descriptions)
        peaks = []      # List of tuples for storing m/z and intensity values)
        spec = None  # Used as Spectrum object
        
        counter = 0
        
        
        # Find a way to keep obo_version updated
        msrun = pymzml.run.Reader(file[0],obo_version = "4.0.1",
                                  extraAccessions=[('MS:1000016', ['value', 'unitName'])])
    
        for spectrum in msrun:
            ID = ""
                
            if spectrum["ms level"] == 2:
                counter+=1
                    
                for tags in spectrum.keys():
                    if tags.startswith("MS") == False:
                        feat_tups = (tags, spectrum[tags])
                        metadata.append(feat_tups)
                            
                    
                if limit == True:
                    if counter > NUM_OF_SPECTRA:
                        #print "%d spectra parsed" % (counter)
                        break
                    
                    
                   
                ID = str(spectrum["id"])
                spec = Spectrum(ID)
                spec.set_group(file[2])
                spec.set_metadata(metadata)
                spec.set_peaks(spectrum.peaks)
                spec.filter_precursor(tolerance)
                spec.normalize_peaks()
                
                
                if library == True:
                    self.library_spec_data.append(spec)
                else:
                    self.spectral_data.append(spec)
        
    
        
    
    
    
    
    def read_mgf(self, file, limit, NUM_OF_SPECTRA, tolerance, library=False):
        metadata = []     # List of tuples for storing metadata a.k.a field descriptions)
        peaks = []      # List of tuples for storing m/z and intensity values)
        spec = None  # Used as Spectrum object
        
        
        fileHandle = open(file[0], "rU")    # Store file handle object in variable
        counter = 0
            
            
        for line in fileHandle:
            is_ID = False
                
            
                
            if limit == True:
                if counter > NUM_OF_SPECTRA:
                    #print "%d spectra parsed" % (counter)
                    break
                
            if line.startswith("BEGIN"):
                counter+=1
            
                #For each spectra, the list of tuples are re-initialized
                metadata = [] 
                peaks =[]
                ID = ""
                    
            
            elif line.startswith("FEATURE") or line.startswith("SPECTRUMID"):
                feature_id = line.split("=")
                ID = feature_id[1][:-1]
                spec = Spectrum(ID)
                is_ID = True
                
                
            
                
            ## Creating metadata List ##
            elif "=" in line and is_ID == False:
                    
                data = line.split("=")
            
                # The slicing notation is to remove the trailing "/n"
                field_tups = (data[0], data[1][:-1])    
                metadata.append(field_tups)
                
            ## Creating peak List ##
            elif line[0].isdigit():
                ion_data = line.split()
                ion_tups = (float(ion_data[0]), float(ion_data[1]))
                peaks.append(ion_tups)
                    

            elif line.startswith("END"):
                spec.set_group(file[2])
                spec.set_metadata(metadata)
                spec.set_peaks(peaks)
                spec.filter_precursor(tolerance)
                spec.normalize_peaks()
                    
                if library == True:
                    self.library_spec_data.append(spec)
                else:
                    self.spectral_data.append(spec)
    
    
    
    
    def parse(self, limit=False, NUM_OF_SPECTRA=10, tolerance=17):
        
       # self.spectral_data = []
        #self.parsed = []
        if len(self.files) == 0:
            raise LookupError("No files have been added")
            
        
        for file in self.files:
            file_format = file[1]
            
            
            if file_format not in self.file_formats:
                raise ValueError("Valid formats are mzml or mgf")
                
            
            if file[0] not in self.parsed:
                self.parsed.append(file[0])
        
                if file_format == "mzml":
                    return self.read_mzml(file, limit, NUM_OF_SPECTRA, tolerance)
            
           
        
                elif file_format == "mgf":
                    return self.read_mgf(file, limit, NUM_OF_SPECTRA, tolerance)
            
                self.parsed.append(file[0])
                
                
    
    def __iter__(self):
        return self
    
    
    def next(self):
        if self.index == len(self.spectral_data):
            self.index = 0
            raise StopIteration
            
        val = self.spectral_data[self.index]
        self.index+=1
        return val
    
    def get_spectrum(self, **params):
        found = False
        for spec in self.spectral_data:
            if params["ID"] == spec.ID and params["group"] == spec.experiment_group:
                found = True
                return spec
                break
        
        if found == False:
            raise NameError("Spectrum not found")
    
    def spectrum_from_library(self, ** params):
        found = False
        for lspec in self.library_spec_data:
            if params["ID"] == lspec.ID and params["group"] == lspec.experiment_group:
                found = True
                return lspec
                break
        
        if found == False:
            raise NameError("Spectrum not found")
    
    def get_consensus(self, **params):
        found = False
        for con in self.consensus_data:
            if params["ID"] == con.ID and params["group"] == con.experiment_group:
                found = True
                return con
                break
        
        if found == False:
            raise NameError("Consensus not found")
    
    
    def __len__(self):
        return len(self.spectral_data) 
    
    
    def spectrum_annotation_file(self, path):
        fileout = open(path,"w")
        header = "ID\tGroup\tPreffered Label\tPrecursor mass\n"
        fileout.write(header)
        
        for spec in self.spectral_data:
            line = "%s.%s\t%s\t%s\t%s" % (spec.ID,spec.experiment_group, spec.experiment_group,
                                          spec.ID, spec.parent_mz)
            fileout.write(line + "\n")
        fileout.close() 
    
    
    def comparison_plot(self, spec1, spec2, path, tolerance=0.3, min_match=6, shift=False):
        
        if shift == False:
            score, used_matches = sf.fast_cosine(spec1, spec2, tolerance, min_match)
        
        elif shift == True:
            score, used_matches = sf.fast_cosine_shift(spec1, spec2, tolerance, min_match)
        
        plt.figure(figsize=(30,10))
        
        spec1_matched_positions = []
        spec2_matched_positions = []
        
        for matches in used_matches:
            spec1_matched_positions.append(matches[0])
            spec2_matched_positions.append(matches[1])
            
        
        
        spec1_sorted_intensity = sorted(spec1.normalised_peaks, key= lambda x:x[1], reverse=True)
        spec2_sorted_intensity = sorted(spec2.normalised_peaks, key= lambda x:x[1], reverse=True)
        
        # Calculating factors for scaling y axis (maximum value is 100)
        spec1_max_intensity = spec1_sorted_intensity[0][1]
        spec2_max_intensity = spec2_sorted_intensity[0][1]
        
        spec1_factor = spec1_max_intensity/100
        spec2_factor = spec2_max_intensity/100
        

        for i,pk in enumerate(spec1.normalised_peaks):
            if i in spec1_matched_positions:
                plt.plot((pk[0],pk[0]), (0,pk[1]/spec1_factor), 'b-')
            else:
                plt.plot((pk[0],pk[0]), (0,pk[1]/spec1_factor), 'r-')

        for j,pk2 in enumerate(spec2.normalised_peaks):
            if j in spec2_matched_positions:
                plt.plot((pk2[0], pk2[0]), (0,-pk2[1]/spec2_factor), 'b-')
            else:
                plt.plot((pk2[0], pk2[0]), (0,-pk2[1]/spec2_factor), 'g-')
                
        sorted_peaks = sorted(spec1.normalised_peaks, key= lambda x:x[0])
        
        xlower = sorted_peaks[0][0]-50
        xupper = sorted_peaks.pop()[0] +50
        title = "%s (%s) vs %s (%s)" % (str(spec1.parent_mz), str(spec1.ID),
                                        str(spec2.parent_mz), str(spec2.ID))

        plt.axis([xlower,xupper,-100,100])
        #plt.xticks(np.arange(50,300, step=5))
        plt.xticks(np.arange(int(xlower),xupper, step=50))
        plt.xlabel("mass-to-charge ratio (m/z)", fontsize=20)
        plt.ylabel("relative intensity", fontsize=20)
        plt.title(title, fontsize=20)
        plt.axhline(0,color="black", lw=2)
        plt.savefig(path)
        plt.close()



    def compute_molecular_network(self, path, tolerance=0.3, min_match=6, pm_range=1.0, edge_filer=0.7):
        
        fileout = open(path,"w")
    
        header = "ID1"+" "+"ID2"+" "+"score"+"\n"
        fileout.write(header + "\n")
    
    
        # This lists holds tuples of indices
        # This list is used to prevent repeating comparisons
        reverse_order = [] 
        
        for con1 in self.consensus_data:
            for con2 in self.consensus_data:
        
            # This condition prevents self comparison and repeating comparisons (happening in reverse order)
                if con1.ID != con2.ID and (con1.ID, con2.ID) not in reverse_order:             
                    
                    reverse_order.append((con2.ID,con1.ID))
                    
                    
                    
                    if abs(con1.parent_mz - con2.parent_mz) <= pm_range:
                        
                        # Calculating cosine score
                        score, used_matches = sf.fast_cosine_shift(con1, con2, tolerance, min_match)
                        
                        if score >= edge_filer:
                            line = "%s.%s %s.%s %s" % (str(con1.ID),str(con1.experiment_group),
                                                       str(con2.ID),str(con2.experiment_group),
                                                       str(score))
                            fileout.write(line +"\n") 
                        
    
    def initialise_consensus(self):
        self.consensus_data = []
        self.clustered = []
        self.unclustered = []
        con = None
        
        for spec in self.spectral_data:
            con = consensus.Consensus(spec)
            self.consensus_data.append(con)
    
    def separate_clusters(self):
        for con in self.consensus_data:
            if len(con) >= 2:
                self.clustered.append(con)
            else:
                self.unclustered.append(con)
                
    def condense_duplicate_ID(self):
        for j,con1 in enumerate(self.consensus_data):
            matched = []
            for k,con2 in enumerate(self.consensus_data[1:]):
                #if con1 != con2:
                if con1.ID == con2.ID and con1.experiment_group == con2.experiment_group:
                    new_con = consensus.merge_consensus(con1,con2)
                    for i,c in enumerate(self.consensus_data):
                        if c==con2:
                            self.consensus_data[i] = new_con
                            #break
                                
                                
                    if con1 not in matched:
                        matched.append(con1)
                elif con1 == con2:
                    self.consensus_data.remove(self.consensus_data[k])
            
            if len(matched)!= 0:
                for cons in matched:
                    if cons in self.consensus_data:
                        self.consensus_data.remove(cons)
    
    def condense_consensus(self, pm_range):
        for con1 in self.consensus_data:
            matched = []
            for con2 in self.consensus_data:
                if con1 != con2 and abs(con1.parent_mz - con2.parent_mz) <= pm_range:
                    new_con = consensus.merge_consensus(con1,con2)
                    self.consensus_data.remove(con2)
                    self.consensus_data.append(new_con)
                    if con1 not in matched:
                                matched.append(con1)
            
            for cons in matched:
                self.consensus_data.remove(cons)
    
    def build_consensus(self, tolerance=0.3, min_match=6, pm_range=0.5, score_threshold=0.9):
        self.initialise_consensus()
        
        new_con = None
        
        reverse_order = []
        reference = []
        r =0
        
       
            
        
        
        for j,con1 in enumerate(self.consensus_data):
            matched = []
            for con2 in self.consensus_data[1:]:
                    
                
            
            # This condition prevents self comparison and repeating comparisons (happening in reverse order)
                if con1.ID != con2.ID and (con1.ID, con2.ID) not in reverse_order:
                            #comparisons+=1
                    
                    reverse_order.append((con2.ID,con1.ID))
                    
                    if abs(con1.parent_mz - con2.parent_mz) <= pm_range:
                        
                            # Calculating cosine score
                        score, used_matches = sf.fast_cosine(con1, con2, tolerance, min_match)
                        if score >= score_threshold:
                            new_con = consensus.merge_consensus(con1,con2)
                            for i,c in enumerate(self.consensus_data):
                                if c==con2:
                                    self.consensus_data[i] = new_con
                                    #break
                                
                                
                            if con1 not in matched:
                                #print self.consensus_data[j]
                                matched.append(con1)
            #print len(matched)                   
            if len(matched)!= 0:
                for cons in matched:
                    self.consensus_data.remove(cons)
                            
                            
                                            
    
            #self.condense_consensus(pm_range)
        self.condense_duplicate_ID()

            #r+=1
            
        self.separate_clusters()
                        
     
    
    def parse_library(self, limit=False, NUM_OF_SPECTRA=10, tolerance=17):
        #self.spectral_data = []
        #self.parsed = []
        if len(self.libraries) == 0:
            raise LookupError("No libraries have been added")
            
        #self.spectral_data = []
        for file in self.libraries:
            file_format = file[1]
            
            
            if file_format not in self.file_formats:
                raise ValueError("Valid formats are mzml or mgf")
                
            
            if file[0] not in self.parsed_libraries:
                self.parsed_libraries.append(file[0])
        
                if file_format == "mzml":
                    return self.read_mzml(file, limit, NUM_OF_SPECTRA, tolerance, library=True)
                
                elif file_format == "mgf":
                    return self.read_mgf(file, limit, NUM_OF_SPECTRA, tolerance, library=True)
    
    
    
    

    def find_library_match(self,path,tolerance=0.3, min_match=6, pm_range=0.5, score_threshold=0.90):
        #self.parse_library()
        if len(self.consensus_data) == 0:
            raise LookupError("Consensus has not yet been built")
        
        elif len(self.library_spec_data) == 0:
            raise LookupError("Library files have not been parsed")
        
        print "Looking for hits..."
        
    
        fileout = open(path,"w")
        header = "ID" + "\t" + "Compound name" + "\t" + "Compound mass" + "\t"+"Library" + "\t" + "Match ID" + "\t"+ "Identification status" + "\t" + "Match score" + "\n"
        
        fileout.write(header)
    
    
        # Each spectral library is parsed and each peak info is iterated over until a match is found
        for con in self.consensus_data:
            for lib_spec in self.library_spec_data:
                if abs(con.parent_mz - lib_spec.parent_mz) <= pm_range:
            
                
                    # Calculating cosine score
                    score, used_matches = sf.fast_cosine(con, lib_spec, tolerance, min_match)
                        
                    
                    # It's only a match if the cosine score above a 0.90 threshold
                    # A report of the match is given in a file
                            
                    if score >= score_threshold:
                        compound_name = lib_spec.get_compound_name()
                            
                        line = "%s.%s\t%s\t%s\t%s\t%s\tIdentified\t%s" % (con.ID, con.experiment_group,
                                                                          compound_name,lib_spec.parent_mz,
                                                                          lib_spec.experiment_group,
                                                                          lib_spec.ID,str(score))
                        
                        fileout.write(line +"\n")
                    
                    #except KeyError:
                        #print "%s does not have peak data" % (l_keys)
            
        print "Finished looking for hits"                   


##################################################################################################################

def store_library_annotation(path):
    filehandle = open(path,"rU")
    header = filehandle.readline()
    d = dict()
    seen = []
    #extra
    
    for line in filehandle:
        data = line.split("\t")
        
        
        if data[0] in seen:
            extras = []
            if type(d[data[0]]) != "list":
                extras.append(d[data[0]])
                d[data[0]] = extras
            else:
                extras.extend(d[data0])
                d[data[0]] = extras
        else:
            seen.append(data[0])
            d.setdefault(data[0], data[1])
                
        
    return d
                