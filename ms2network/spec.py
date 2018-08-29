# To help normalize the peaks
import math
from collections import namedtuple
Peak = namedtuple('Peak',['mz','intensity'])

# Matplot lib and numpy packages are needed for ploting spectrum
import matplotlib.pyplot as plt
import numpy as np



class Spectrum(object):
    def __init__(self, ID):
        self.ID = ID
        self.experiment_group = None
        self.metadata = None
        self.parent_mz = None
        self.peaks = None
        self.n_peaks = None
        self.normalised_peaks = None
        self.top_peaks = None
        
    
    def set_group(self, experiment_group):
        self.experiment_group = experiment_group
        
    def set_metadata(self, metadata):
        self.metadata = metadata
        
        for feature in metadata:
            if feature[0] == "selected ion m/z":
                self.parent_mz = float(feature[1])
                
            elif feature[0] == "PEPMASS":
                self.parent_mz = float(feature[1])
                
        #self.parent_mz = parent_mz
                
        #self.parent_mz = parent_mz
    
    def convert_to_peaks(self, peak_tuples):
        #using the splat we can handle both size 2 lists and tuples
        return [Peak(*p) for p in peak_tuples]

    def sqrt_normalize_spectrum(self, spectrum):
        output_spectrum = []
        intermediate_output_spectrum = []
        acc_norm = 0.0
        for s in spectrum:
            sqrt_intensity = math.sqrt(s.intensity)
            intermediate_output_spectrum.append(Peak(s.mz,sqrt_intensity))
            acc_norm += s.intensity
        normed_value = math.sqrt(acc_norm)
        for s in intermediate_output_spectrum:
            output_spectrum.append(Peak(s.mz,s.intensity/normed_value))
        return output_spectrum
    
    
    
    
    
    def filter_precursor(self, tolerance):
        for peak in self.peaks:
            if abs(peak[1] - self.parent_mz) <= tolerance:
                self.peaks.remove(peak)
    
    
    def set_peaks(self, peaks):
        self.peaks = peaks
        self.n_peaks = len(peaks)
        
    def normalize_peaks(self):
        self.normalised_peaks = self.sqrt_normalize_spectrum(self.convert_to_peaks(self.peaks))
        sorted_intensity = sorted(self.normalised_peaks, key= lambda x:x[1], reverse=True)
        self.top_peaks = sorted_intensity[:5]
        
    def get_compound_name(self):
        for feature in self.metadata:
            if feature[0] == "NAME":
                name = feature[1]
        
        return name

    def plot_spectrum(self, path):
        
        sorted_intensity = sorted(self.normalised_peaks, key= lambda x:x[1], reverse=True)
        
        # Calculating factor for scaling y axis (maximum value is 100)
        max_intensity = sorted_intensity[0][1]
        scaling_factor = max_intensity/100
        
        plt.figure(figsize=(30,10))

        for pk in self.normalised_peaks:
            plt.plot((pk[0],pk[0]), (0,pk[1]/scaling_factor), color="black")
        
        sorted_peaks = sorted(self.normalised_peaks, key= lambda x:x[0])
        
        xlower = sorted_peaks[0][0]-20
        xupper = sorted_peaks.pop()[0] +30
        
        
        plt.axis([xlower,xupper,0,100])
        plt.xticks(np.arange(int(xlower),xupper, step=50))
        plt.xlabel("mass-to-charge ratio (m/z)", fontsize=20)
        plt.ylabel("relative intensity", fontsize=20)
        plt.title(str(self.parent_mz), fontsize=20)
        plt.savefig(path)
        plt.close()
        
        
    
    
    def total_peak_intensity(self):
        intensity = 0
        for peak in self.normalised_peaks:
            intensity+=peak[1]
               
        return intensity