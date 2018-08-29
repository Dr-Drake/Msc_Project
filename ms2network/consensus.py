from spec import Spectrum







class Consensus(object):
    def __init__(self, spec):
        self.ID = spec.ID
        self.parent_mz = spec.parent_mz
        #self.pm_range = 0.5
        #self.score_threshold = 0.95
        self.n_peaks = spec.n_peaks
        self.normalised_peaks = spec.normalised_peaks
        self.experiment_group = spec.experiment_group
        self.spectral_data = [spec]
        self.total_peak_intensity = spec.total_peak_intensity()
        self.top_peaks = spec.top_peaks
        self.index = 0
        
    def __iter__(self):
        return self
    
    
    def __next__(self):
        if self.index == len(self.spectral_data):
            self.index = 0
            raise StopIteration
            
        val = self.spectral_data[self.index]
        self.index+=1
        return val
    
    
    
    #def set_pm_range(self, value):
    #    self.pm_range = value

    
        
    def add(self, new_spec):
        self.spectral_data.append(new_spec) 
        
    def __contains__(self,spec):
        for s in self.spectral_data:
            if spec.ID == s.ID and spec.experiment_group == s.experiment_group:
                val = True
                break
            else:
                val = False
        
        return val
    
    #def within_range(self, spec1):
        #eligible = False
       # for spec2 in self.spectral_data:
        #    if abs(spec1.parent_mz - spec2.parent_mz) <= self.pm_range:
           #     eligible = True
        #    else:
        #        eligible = False
        #        break
                
       # return eligible
    
    def update(self):
        #highest_intensity = self.spectral_data[0].highest_peak_intensity()
        #consensus_ID = self.spectral_data[0].ID
        
        for spec in self.spectral_data:
            if spec.total_peak_intensity() > self.total_peak_intensity:
                self.total_peak_intensity = spec.total_peak_intensity()
                self.ID = spec.ID
                self.parent_mz = spec.parent_mz
                self.normalised_peaks = spec.normalised_peaks
                self.n_peaks = spec.n_peaks
                self.experiment_group = spec.experiment_group
                self.top_peaks = spec.top_peaks
                
                
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
                
                
    def __len__(self):
        return len(self.spectral_data)
		
				
#################################################################################################

def merge_consensus(con1, con2):
    spectral_data = []
    sdata_set = set()
    new_con = None
    
    spectral_data.extend(con1.spectral_data)
    spectral_data.extend(con2.spectral_data)
    
    sdata_set.update(spectral_data)
    
    spectral_data = list(sdata_set)
    new_con = Consensus(spectral_data.pop())
    
    for s in spectral_data:
        new_con.add(s)
    
    new_con.update()
    
    return new_con