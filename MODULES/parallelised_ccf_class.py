import numpy as np
import multiprocessing
	
class ccf_:

	def __init__(self):
		#self.data_arrays = data_arrays
		pass
    
	def shift_flux_fn(self, data_arrays):
	
		def getnearpos(array,value):
			array = np.asarray(array)
			idx = (np.abs(array-value)).argmin()
			return idx 
	
		
		# unpack arrays
		flux_mask, wlen_obs, wlen_mask = data_arrays
		
		# perform shift 
		flux_shifted = [flux_mask[getnearpos(wlen_obs, w)] \
						for w in wlen_mask]
		
		# return shifted array
		flux_shifted = np.asarray(flux_shifted)
		return flux_shifted