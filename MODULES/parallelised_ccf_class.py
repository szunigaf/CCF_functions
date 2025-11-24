import numpy as np

class ccf_:
    """
    CCF utility class for spectral flux shifting operations.
    
    Note: This class provides methods for CCF calculations that are designed
    to be used with external multiprocessing (via pathos.multiprocessing).
    The actual parallelization is handled in calc_rv_new_chapman.py.
    """

    def __init__(self):
        """Initialize the CCF utility class."""
        pass
    
    def shift_flux_fn(self, data_arrays):
        """
        Shift flux data to match wavelength mask.
        
        This method is designed to be called by external multiprocessing
        (via pathos.multiprocessing.ProcessingPool).
        
        Parameters:
        -----------
        data_arrays : tuple
            Tuple containing (flux_mask, wlen_obs, wlen_mask)
            
        Returns:
        --------
        numpy.ndarray
            Shifted flux array matched to wavelength mask
        """
        
        def getnearpos(array, value):
            """Find index of nearest value in array."""
            array = np.asarray(array)
            idx = (np.abs(array - value)).argmin()
            return idx 

        # unpack arrays
        flux_mask, wlen_obs, wlen_mask = data_arrays
        
        # perform shift 
        flux_shifted = [flux_mask[getnearpos(wlen_obs, w)] 
                        for w in wlen_mask]
        
        # return shifted array
        flux_shifted = np.asarray(flux_shifted)
        return flux_shifted