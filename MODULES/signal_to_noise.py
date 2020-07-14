import sys
sys.path.append("./MODULES/")

import matplotlib as mpl

#mpl.use('Agg')
import math
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pf
								  
from matplotlib import rc
import matplotlib as mpl
import glob


mpl.rcParams['legend.numpoints'] = 1

rc('text', usetex=True)
rc('font',**{'family':'serif','serif':['Times'], 'size':'12','weight':'bold'})

def poisson(k, lamb):
    return (lamb**k/factorial(k)) * np.exp(-lamb)

def getnearpos(array,value):
	array = np.asarray(array)
	idx = (np.abs(array-value)).argmin()
	return idx     

def ReadSpec_UVES_fits(fits_file):

	"""Function to read in fits file of high
	resolution spectra and then return wlen 
	and flux arrays"""
	
	# open fits file
	image=pf.open(fits_file)
	
	# retrieve flux and header
	flux=image[0].data
	header=image[0].header
	
	# create wavelength array
	wlen=-header['CRPIX1']*header['CDELT1'] + \
		  header['CRVAL1'] + \
		  np.arange(1,len(flux)+1)*header['CDELT1']
		  
	image.close()
	
	return wlen, flux

direc = "/Users/pelliott/Desktop/STANDALONE_CCF_ROUTINE/UVES_FITS/"
fits_file = "WXCol_2007-07-31T09:33:18.898_MERGED_SCIENCE_REDL.fits"



def SigtoNoise_calc(wlen, flux):
	print "\nWlen range: %i - %i\n" % (min(wlen), max(wlen))
	print "Wlen Step: %.6f" % (wlen[1]-wlen[0])
	
	wlen_mid = np.mean([max(wlen), min(wlen)])
	
	wlen_mid_min, wlen_mid_max = wlen_mid - 200., wlen_mid + 200.
	
	idx_min = getnearpos(wlen, wlen_mid_min)
	idx_max = getnearpos(wlen, wlen_mid_max)
	
	# make a new arrays of just the sub section around the 
	# mid point of the wlen
	wlen, flux = wlen[idx_min:idx_max], flux[idx_min:idx_max]
	
	
	wlen_step = wlen[1]-wlen[0]
	angstrom_chunk = 1.
	box_size = int(angstrom_chunk / wlen_step)
	
	print "Will use pixel chunking of %i, for %i A" % (box_size, angstrom_chunk)
	
	w=0
	
	wlen_chunked, snr_array=[], []
	while w < (len(wlen) - box_size):
		
		if np.median(flux[w:w+box_size]) > 0:
			snr = (np.median(flux[w:w+box_size]) / np.std(flux[w:w+box_size]))
		
			if math.isnan(snr) != True:
				snr_array.append(snr)
				wlen_chunked.append(np.median(wlen[w:w+box_size]))
		
		w=w+box_size
	
# 	lin_bins = np.arange(0, 100, 2)
# 	plt.hist(snr_array, bins = lin_bins)
# 	plt.show()
	print "Median: %.1f, Mean: %.1f" % (np.median(snr_array), np.mean(snr_array))
	
	return np.median(snr_array)
	



wlen, flux = ReadSpec_UVES_fits(direc+fits_file)
snr = SigtoNoise_calc(wlen, flux)


direc = "/Users/pelliott/Desktop/STANDALONE_CCF_ROUTINE/UVES_FITS/"
fits_list = glob.glob(direc+"*fits")


master_output = "/Users/pelliott/Desktop/STANDALONE_CCF_ROUTINE/OUTPUTS/signal_to_noise.txt"

for fits_file in fits_list:
	print "%i / %i" % (fits_list.index(fits_file)+1, len(fits_list))
	wlen, flux = ReadSpec_UVES_fits(fits_file)
	snr = SigtoNoise_calc(wlen, flux)
	
	output_string = "%s	%.1f\n" % (fits_file.replace(direc, ""), snr)
	
	text_file = open(master_output, "a")
	text_file.write(output_string)
	text_file.close()
	
	
	





