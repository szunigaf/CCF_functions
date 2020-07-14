import sys
sys.path.append("./MODULES/")

import matplotlib as mpl

#mpl.use('Agg')

import numpy as np
import matplotlib.pyplot as plt
import glob
import scipy as sp
import astropy.io.fits as pf
from scipy import interpolate
from scipy.stats import chisquare
from scipy.optimize import curve_fit
import time
import parallelised_ccf_class
import pathos.multiprocessing as mp
from scipy.stats import linregress
								  
from matplotlib import rc
import matplotlib as mpl
mpl.rcParams['legend.numpoints'] = 1

rc('text', usetex=True)
rc('font',**{'family':'serif','serif':['Times'], 'size':'12','weight':'bold'})


# global parameters
# *******************************************

# speed of light in km/s
c = 299792.458

# define velocity range of potential RVs
#vel_step = 10.
vel_step = 1.
vel_range = np.arange(-250, 251, vel_step)
#vel_range = np.arange(-50, 51, vel_step)

# vsini step - that of the Gray profiles
#vsini_step = 10.
vsini_step = 1

	
# resolution of wavelength grid in Angstroms
# better resolution significantly slows down
# the CCF calculation
wlen_step = 0.02
# best value as a compromise of resolution
# and time is 0.02


# Number of CPUs to use
NCPU = 8

# *******************************************




def add_subplot_axes(ax,rect,axisbg='w'):
	fig = plt.gcf()
	box = ax.get_position()
	width = box.width
	height = box.height
	inax_position  = ax.transAxes.transform(rect[0:2])	
	transFigure = fig.transFigure.inverted()
	infig_position = transFigure.transform(inax_position)	
	x = infig_position[0]
	y = infig_position[1]
	width *= rect[2]
	height *= rect[3]  # <= Typo was here
	subax = fig.add_axes([x,y,width,height])#,axisbg=axisbg)
	x_labelsize = subax.get_xticklabels()[0].get_size()
	y_labelsize = subax.get_yticklabels()[0].get_size()
	x_labelsize *= rect[2]**0.5
	y_labelsize *= rect[3]**0.5
	subax.xaxis.set_tick_params(labelsize=x_labelsize)
	subax.yaxis.set_tick_params(labelsize=y_labelsize)
	return subax



def apply_sigma_cut(input_array, sigma):
	# get 3 sigma values
	upper_val = np.median(input_array) + (sigma * np.std(input_array))
	lower_val = np.median(input_array) - (sigma * np.std(input_array))
	
	# apply clipping
	output_array = np.clip(input_array, lower_val, upper_val)	
	
	return output_array
	

def getnearpos(array,value):
	t = time.time()
	array = np.asarray(array)
	idx = (np.abs(array-value)).argmin()
	return idx 	


def get_mask(SpT):
	mask_direc = "./masks/"
	# find mask and return the wlen + profile arrays
	if glob.glob(mask_direc+SpT+".mas") != []:
	
		data = np.loadtxt(mask_direc+SpT+".mas", unpack=True)
		wlen1, wlen2, prof, wlen_av = data
		
	else:
		print "\nMask not found for Spt: %s\n"	% SpT
	
	wlen_arr, flux_arr = [], []
	for w in np.arange(0, len(wlen1)-1, 1):
		
		wlen_arr.append(wlen1[w])
		wlen_arr.append(wlen2[w])
		wlen_arr.append(wlen2[w]+0.001)
		wlen_arr.append(wlen1[w+1]-0.001)
		
		flux_arr.append(1)
		flux_arr.append(1)
		flux_arr.append(0)
		flux_arr.append(0)

	
	return wlen_arr, flux_arr


def get_mask_early_type(SpT):
	mask_direc = "/Users/pelliott/ASTRO/Data/CORREL/cors1d_py/MASKS/"

	# find mask and return the wlen + profile arrays
	if glob.glob(mask_direc+SpT+".mas") != []:
	
		data = np.loadtxt(mask_direc+SpT+".mas", unpack=True)
		wlen, x, x2, flux = data
		
	return wlen, flux





def ReadSpec_UVES_fits(fits_file):

	"""Function to read in fits file of high
	resolution spectra and then return wlen 
	and flux arrays"""
	
	# open fits file
	image=pf.open(fits_file)
	
	# retrieve flux and header
	flux=image[0].data
	header=image[0].header
	
	# get object name
	star_name = header['HIERARCH ESO OBS TARG NAME'].replace(" ", "")
	
	# get coords
	ra_val = header['RA']
	dec_val = header['DEC']
	
	# create wavelength array
	wlen=-header['CRPIX1']*header['CDELT1'] + \
		  header['CRVAL1'] + \
		  np.arange(1,len(flux)+1)*header['CDELT1']
	
	# file type to distinguish between 	  
	file_type = header['HIERARCH ESO PRO CATG']
	
	# date of observations
	mjd_obs = header['MJD-OBS']
		  
	# get barycentric correction	  
	bary_correc = header['HIERARCH ESO QC VRAD BARYCOR']	  
		  
	image.close()
	
	# define star_name again with file type and obs date
	star_name = "%s_%s" % (star_name.replace(" ", ""), \
	file_type)
	
	return star_name, ra_val, dec_val, mjd_obs, wlen, flux, bary_correc
	

def combine_spec(fits1, fits2):

	"""Function to read in fits file of high
	resolution spectra and then return wlen 
	and flux arrays"""
	
	# open fits file
	image1=pf.open(fits1)
	image2=pf.open(fits2)
	
	# retrieve flux and header
	flux1=image1[0].data
	header1=image1[0].header
	
	flux2=image2[0].data
	header2=image2[0].header
		
	
	# get object name
	star_name = header1['HIERARCH ESO OBS TARG NAME'].replace(" ", "")
	
	# get coords
	ra_val = header1['RA']
	dec_val = header1['DEC']
	
	# create wavelength array
	wlen1=-header1['CRPIX1']*header1['CDELT1'] + \
		  header1['CRVAL1'] + \
		  np.arange(1,len(flux1)+1)*header1['CDELT1']
		  
	wlen2=-header2['CRPIX1']*header2['CDELT1'] + \
		  header2['CRVAL1'] + \
		  np.arange(1,len(flux2)+1)*header2['CDELT1']
	
	# concatenate arrays
	wlen = np.concatenate((wlen1, wlen2), axis=0)
	flux = np.concatenate((flux1, flux2), axis=0)
	
	
	# file type to distinguish between 	  
	file_type = header1['HIERARCH ESO PRO CATG']
	
	# date of observations
	mjd_obs = header1['MJD-OBS']
		  
	# get barycentric correction	  
	bary_correc = header1['HIERARCH ESO QC VRAD BARYCOR']	  
		  
	image1.close()
	image2.close()
	
	# define star_name again with file type and obs date
	star_name = "%s_%s" % (star_name.replace(" ", ""), \
	file_type)
	
	return star_name, ra_val, dec_val, mjd_obs, wlen, flux, bary_correc	


def ReadSpec_FEROS_fits(fits_file):

	"""Function to read in fits file of high
	resolution spectra and then return wlen 
	and flux arrays"""
	
	# open fits file
	image=pf.open(fits_file)
	
	# retrieve flux and header
	flux=image[0].data
	header=image[0].header
	
	# get object name
	star_name = header['OBJECT']
	
	# get coords
	ra_val = header['RA']
	dec_val = header['DEC']
	
	# create wavelength array
	wlen=-header['CRPIX1']*header['CDELT1'] + \
		  header['CRVAL1'] + \
		  np.arange(1,len(flux)+1)*header['CDELT1']
		  
	print "Min. - Max. wlen: %i - %i \n" % (wlen[0], wlen[-1])
	
	# file type to distinguish between 	  
	file_type = header['HIERARCH ESO PRO CATG']
	
	# date of observations
	mjd_obs = header['MJD-OBS']
		  
	# barycentric correction already applied for FEROS
	bary_correc = 0
		  
	image.close()
	
	# define star_name again with file type and obs date
	star_name = "%s_%s" % (star_name.replace(" ", ""), \
	file_type)
	
	return star_name, ra_val, dec_val, mjd_obs, wlen, flux, bary_correc

def ReadSpec_UVES_PH3_fits(fits_file):

	"""Function to read in fits file of high
	resolution spectra and then return wlen 
	and flux arrays"""
	
	# open fits file
	image=pf.open(fits_file)
	
	# retrieve flux and header
	#data = image[1].data
	wlen, flux, err = image[1].data['WAVE'][0],image[1].data['FLUX'][0],image[1].data['ERR'][0]
	
	header=image[0].header
	
	# get object name
	star_name = header['OBJECT']
	
	# get coords
	ra_val = header['RA']
	dec_val = header['DEC']
	
	print "Min. - Max. wlen: %i - %i \n" % (wlen[0], wlen[-1])
	
	# file type to distinguish between 	  
	file_type = header['PRODCATG']
	
	# date of observations
	mjd_obs = header['MJD-OBS']
		  
	# barycentric correction already applied for FEROS
	bary_correc = header['HIERARCH ESO QC VRAD BARYCOR']
	#bary_correc = 0
		  
	image.close()
	
	# define star_name again with file type and obs date
	star_name = "%s_%s" % (star_name.replace(" ", ""), \
	file_type)
	
	print star_name, ra_val, dec_val, mjd_obs, wlen, flux, bary_correc	
	
	return star_name, ra_val, dec_val, mjd_obs, wlen, flux, bary_correc	


	

def ReadSpec_HARPS_fits(fits_file):

	"""Function to read in fits file of high
	resolution spectra and then return wlen 
	and flux arrays"""
	
	# open fits file
	image=pf.open(fits_file)
	
	# retrieve flux and header
	data = image[1].data
	wlen, flux, err = data[0]
	
	header=image[0].header
	
	# get object name
	star_name = header['OBJECT']
	
	# get coords
	ra_val = header['RA']
	dec_val = header['DEC']
	
	print "Min. - Max. wlen: %i - %i \n" % (wlen[0], wlen[-1])
	
	# file type to distinguish between 	  
	file_type = header['PRODCATG']
	
	# date of observations
	mjd_obs = header['MJD-OBS']
		  
	# barycentric correction already applied for FEROS
	bary_correc = 0
		  
	image.close()
	
	# define star_name again with file type and obs date
	star_name = "%s_%s" % (star_name.replace(" ", ""), \
	file_type)
	
	print star_name, ra_val, dec_val, mjd_obs, wlen, flux, bary_correc	
	
	return star_name, ra_val, dec_val, mjd_obs, wlen, flux, bary_correc	




	
def ReadSpec_txt(txt_file):

	"""Read in spectrum from tab separated txt file in Angstroms"""
	
	data = np.loadtxt(txt_file, unpack=True)
#	wlen, flux, err = data
	wlen, flux = data
	err = 0.05*flux
	
	star_name = txt_file.split('/')[-1].replace(".txt", "")
	
	return star_name, wlen, flux, 0
	


#def make_ccf(wlen_obs, flux_obs, wlen_mask, flux_mask, \
#			bary_correc, wlen_unit='A', plot_flag=False):
def make_ccf(wlen_obs, flux_obs, wlen_mask, flux_mask, \
			bary_correc, wlen_unit='A', plot_flag=False):
	
	# wlen_unit is important should either be 'A' for Angstrom
	# or 'N' for nanometers.  Mask should always be in Angstroms
	
	# for mask
	wlen_mask_km = np.asarray(wlen_mask) * 1e-13
	
	# for observed spectrum
	if wlen_unit == 'A':
		wlen_obs_km = np.asarray(wlen_obs) * 1e-13
	if wlen_unit == 'N':
		wlen_obs_km = np.asarray(wlen_obs) * 1e-12
	
	# first get the range of wlen values that match
	# between the mask and the observed spectrum
	
	abs_min = max([min(wlen_obs), min(wlen_mask)])
	abs_max = min([max(wlen_obs), max(wlen_mask)])
	
	
	# create reduced arrays using min and max values
	# for both observed and mask spectra
	wlen_obs_new, flux_obs_new = [], []
	for w_obs, f_obs in zip(wlen_obs, flux_obs):
		if w_obs > abs_min and w_obs < abs_max:
			wlen_obs_new.append(w_obs)
			flux_obs_new.append(f_obs)
	
	
	wlen_mask_new, flux_mask_new = [], []
	for w_mask, f_mask in zip(wlen_mask, flux_mask):
		if w_mask > abs_min and w_mask < abs_max:
			wlen_mask_new.append(w_mask)
			flux_mask_new.append(f_mask)

	# define abs min and max again to avoid 
	# interpolation	problem
	abs_min = max([wlen_mask_new[0], wlen_obs_new[0]])
	abs_max = min([wlen_mask_new[-1], wlen_obs_new[-1]])

	# create a linear wlen array using min and max values
	wlen_lin = np.arange(abs_min, abs_max, wlen_step)	
	
	# define central wavelength and wavelength range
	central_wlen = wlen_lin[int(len(wlen_lin) / 2.)]
	wlen_range = abs_max - abs_min

	
	# print minimum and maximum wavelength
	print "\nWill compute CCF using data \
	\nbetween %.1f and %.1f A\n" % (abs_min, abs_max)
	
	# perform interpolation for obs spectrum
	f = interpolate.interp1d(wlen_obs_new, flux_obs_new, 'linear')
	flux_obs_lin = f(wlen_lin)
	
	# 'flatten' spectrum
	f = np.poly1d(np.polyfit(wlen_lin, flux_obs_lin, 6))
	poly_flux = f(wlen_lin)
	
	# divide by poly_flux to 'flatten' spectrum
	flux_obs_lin = flux_obs_lin / poly_flux
	
	# set the continuum to 1
	median_flux = np.median(flux_obs_lin)
	flux_shift_val = 1 - median_flux
	flux_obs_lin = np.asarray(flux_obs_lin) + flux_shift_val
	
	
	# plot spectrum if plot_flag is True
	if plot_flag == True:
		plt.plot(wlen_lin, flux_obs_lin, alpha=0.5)
		plt.show()
	
	
	# perform interpolation for mask	
	f = interpolate.interp1d(wlen_mask_new, flux_mask_new, 'linear')
	flux_mask_lin = f(wlen_lin)
	
	# make sure values are zero or 1 as the interpolation
	# can create values between
	flux_mask_lin_new=[]
	for f in flux_mask_lin:
		if f != 0 and f != 1:
			flux_mask_lin_new.append(0)
		else:
			flux_mask_lin_new.append(f)
	
	# move 'new' array to the original mask flux array
	flux_mask_lin = flux_mask_lin_new
	
	
	
	# calculate signal to noise of observed spectrum
	# using area of the spectrum where the mask is set to zero
	non_abs_spec, non_abs_wlen=[], []
	for f_obs, f_mask, w in zip(flux_obs_lin, flux_mask_lin, wlen_lin):
		if f_mask == 0:
			non_abs_spec.append(f_obs)
			non_abs_wlen.append(w)
	
	# define bin size
	bin_size = int(len(non_abs_spec) / 50.)

	
	# get central wavelength
	cent_wlen = wlen_lin[int(len(wlen_lin)/2)]
	
	# create array using binning
	t = 0
	snr_arr, wlen_binned_arr = [], []
	while t < len(non_abs_spec):
		# SNR value
		if np.std(non_abs_spec[t:t+bin_size]) != 0:
			snr_val = np.mean(non_abs_spec[t:t+bin_size]) / np.std(non_abs_spec[t:t+bin_size])
		else:
			snr_val = 0
		wlen_bin_av = np.mean(non_abs_wlen[t:t+bin_size])
		
			
		# append to array
		snr_arr.append(snr_val)
		wlen_binned_arr.append(wlen_bin_av)
		
		# add bin_size on to t
		t = t + bin_size
	
	idx_cent_wlen = getnearpos(wlen_binned_arr, cent_wlen)
	
	# take value at central wavelength
	snr_median = np.median(snr_arr)
	
	# print SNR
	print "\nS/N ratio estimated at: %.1f..." % snr_median
	
	# convert the wavelength to km
	wlen_lin_km = np.asarray(wlen_lin) * 1e-13	
	
	
	# now for each velocity in that range
	all_xcorr, all_xcorr_c = [], []
	print "\nNow will compute CCF for each velocity, this may take some time...\n"
	
			  	
	print "Will go from %i to %i km/s in velocity..." % (vel_range[0], vel_range[-1])

		  	
	for vel in vel_range:
		# make sure it's floats
		vel = np.float(vel)
		
		
		# print out something every 20 km/s
		if vel % 20 == 0:
			print "Velocity: %i km/s" % vel
			
		# make sure the arrays are numpy arrays
		flux_mask_lin = np.asarray(flux_mask_lin, dtype=float)
		wlen_lin = np.asarray(wlen_lin, dtype=float)
		flux_obs_lin = np.asarray(flux_obs_lin, dtype=float)
		

# ******************************************************************************************
# ******************************************************************************************
# ******************************************************************************************
#

		#Python way to make CCF		
		t0=time.time()	
		
		#make shifted wavelength array		
		wlen_mask_km_shifted = \
		np.asarray([np.float(w)-((vel * w)/c) for w in wlen_lin_km])


		# split up array for parallelisation
		N = int(len(wlen_mask_km_shifted)/NCPU)
        
        # make chunks of arrays for each input array
		wlen_mask_chunk = [wlen_mask_km_shifted[i:i + N] \
		for i in range(0, len(wlen_mask_km_shifted), N)]
        
		wlen_obs_chunk = [wlen_lin_km[i:i + N] \
		for i in range(0, len(wlen_lin_km), N)]
        
		flux_mask_chunk = [flux_mask_lin[i:i + N] \
		for i in range(0, len(flux_mask_lin), N)]
		
		# make instance of ccf object
		obj = parallelised_ccf_class.ccf_()
			        
        
		# make argument inputs consisting of the chunked arrays
		Toparallel = []
		for i in range(len(wlen_mask_chunk)):			
			Toparallel.append((flux_mask_chunk[i], \
							   wlen_obs_chunk[i], \
							   wlen_mask_chunk[i]))
		
	
		# use map to perform parallelization
		flux_shifted_=[]
		p = mp.ProcessingPool(len(Toparallel))
		flux_shifted_ = p.map(obj.shift_flux_fn, Toparallel)
		
		# remove last element (not sure why an array [1.]
		# is being added to the output array)
		flux_shifted_ = flux_shifted_[:-1]
		
		# clear processes 
		p.clear()
		
		
		# make it into a numpy array and flatten
		flux_shifted_ = np.asarray(flux_shifted_)
		flux_shifted_ = flux_shifted_.flatten()
						
		t2 = time.time() 
		
		#do cross correlation
		xcorr_python = \
		[(f_obs*f_mask) for f_obs, f_mask in zip(flux_obs_lin, flux_shifted_)]

		# sum array
		xcorr_sum_python = sum(xcorr_python)
                #print xcorr_sum_python
#
#
# ******************************************************************************************
# ******************************************************************************************
# ******************************************************************************************


		# append sum of xcorr to master array
		all_xcorr.append(xcorr_sum_python)
		
	
	
	# normalise by the median of the CCF i.e. where
	# the match is 'perfect' i.e. at 1
        #print norm_xcorr
	norm_xcorr = all_xcorr / np.nanmedian(all_xcorr)
	norm_xcorr = (norm_xcorr - np.nanmedian(norm_xcorr)) * -1
	
	# apply barycentric correction
	vel_range_correc = vel_range + bary_correc
	
	
	# plot		
	if plot_flag:
		plt.plot(vel_range_correc, norm_xcorr)
		plt.show()
	
		
	# return arrays
	return vel_range_correc, norm_xcorr, snr_median, central_wlen, wlen_range




def get_gray_profiles(min_vel, max_vel):
	
	"""Return an array of all the Gray profiles found in the 
	defined 'gray_direc' directory"""

	# directory of gray profiles (profiles already convolved with gaussians)
	gray_direc = "./vsini_templates/"
	
	all_vsini_vels, all_vsini_profiles, all_vsini_values=[], [], []
	
	for gray_file in glob.glob(gray_direc+"qqcc*.dat"):
		# read in data as string
		data = np.loadtxt(gray_file, unpack=True, dtype='S20')
		vsini_vel, vsini_profile = data
		
		vsini_profile_value = int(gray_file[-7:-4])
		
		# convert weird string with 'D' to float type
		# also do (value - 1) * -1 so that both the CCF profile
		# and the v sin i profiles have continua set to zero
		# and positive increases 
		
		vsini_profile = [(float(v.replace("D", "e"))-1) * -1  for v in vsini_profile]
		vsini_vel = [float(v) for v in vsini_vel]
		
		# only keep values between min and max velocities
		vsini_vel_new, vsini_profile_new = [], []
		for vel, prof in zip(vsini_vel, vsini_profile):
			if vel >= min_vel and vel <= max_vel:
				vsini_vel_new.append(vel)
				vsini_profile_new.append(prof)

		# append to output_arrays
		all_vsini_vels.append(vsini_vel_new)
		all_vsini_profiles.append(vsini_profile_new)
		all_vsini_values.append(vsini_profile_value)
	
	
	# return the output arrays
	return all_vsini_vels, all_vsini_profiles, all_vsini_values
		


def get_vsini(vel, ccf):
	
	"""A function that fits the inputted
	ccf profile with a set of Gray profiles
	in order to determine the vsini(i) value"""
	
	# 'stretch' CCF to match with GRAY profiles
	stretch_factor = (1 - max(ccf)) / max(ccf)
	ccf_new = ccf * stretch_factor
	
	# use ccf_new as ccf
	ccf = ccf_new
	
	
	# get the maximum of the ccf profile
	# which should correspond to the radial velocity
	idx = list(ccf).index(max(ccf))
	rv = vel[idx]
	
	# shift velocity
	vel_shifted = vel - rv
	
	# get minimum and maximum velocities
	min_vel, max_vel = int(vel_shifted[0]), int(vel_shifted[-1])
	
	# use min and max values in the retrieval of the
	# Gray profiles
	vsini_vels, vsini_profiles, \
	vsini_values = get_gray_profiles(min_vel, max_vel)
	
	residuals_sum_all, stretch_factors=[], []
	for vsini_prof, vsini_vel, vsini_val in zip(vsini_profiles, vsini_vels, vsini_values):
				
		# make ccf on the same grid as the vsini profile
		# using linear interpolation

		lin_vel = np.arange(min_vel, max_vel+1, vsini_step)
		f = interpolate.interp1d(vel_shifted, ccf, 'linear')
		ccf_lin = f(lin_vel)
		
		
		# get the minimum and maximum velocities to consider in the fit
		zero_point_idx = list(vsini_prof).index(max(vsini_prof))

		
		# for max. wavelength		
		r = zero_point_idx
		new_max_vel=[]
		while r <= len(vsini_prof)-1:
			if vsini_prof[r] < 0.01:
				new_max_vel.append(r)
			r=r+1
			
		if new_max_vel != []:
			new_max_vel_val = min(new_max_vel)
		else:
			new_max_vel_val = len(vsini_prof)-1
		
		
		
		# for minimum wavelength
		r = zero_point_idx		
		new_min_vel=[]
		while r >= 1:
			if vsini_prof[r] < 0.01:
				new_min_vel.append(r)
			r=r-1
		
		if new_min_vel != []:
			new_min_vel_val = max(new_min_vel)
		else:
			new_min_vel_val = 0
		
		# now use min_vel_val and max_vel_val to constrain
		# how much of the profile is fitted
		new_min_vel_val, new_max_vel_val = int(new_min_vel_val), \
											int(new_max_vel_val)
		
		
	
		# make array of residuals squared for each 'stretch' factor
		stretch_fac = np.arange(0.1, 10, 0.1)
		
		# now perform fit
		residuals_sum=[]
		for x in stretch_fac:
			vsini_prof_stretch = np.asarray(vsini_prof) * x
			
			prof_reduced = vsini_prof_stretch[new_min_vel_val:new_max_vel_val]
			ccf_reduced = ccf_lin[new_min_vel_val:new_max_vel_val]
			
			residuals = [(p-c)**2 for p, c in zip(prof_reduced, ccf_reduced)]
			residuals_sum.append(sum(residuals))
		
		
		# get minimum of residuals for each stretch factor
		idx_min = list(residuals_sum).index(min(residuals_sum))
		
		# get 'best' stretch factor
		best_stretch = stretch_fac[idx_min]
			
		# now append to master arrays
		residuals_sum_all.append(residuals_sum[idx_min])
		stretch_factors.append(best_stretch)
	
	
	# print warning message of the lengths don't match
	if len(ccf_lin) == len(vsini_prof):
		print "\nNow computing residuals from rotational profile fitting...\n"
	else:
		print "\nWarning: CCF not on the same grid as vsini profile\n"	
	
	# Get best value from the minimum of the residuals with stretch factor
	idx_vsini = list(residuals_sum_all).index(min(residuals_sum_all))
	min_resid_val = min(residuals_sum_all)	
	
	stretch_fact = stretch_factors[idx_vsini]
	calc_vsini_val = vsini_values[idx_vsini]
# 	
# 	print "\nBest fit vsini value: %i km/s\n" % calc_vsini_val
	
	#return calc_vsini_val, min_resid_val
	
	return vel_shifted, ccf, vsini_vels[idx_vsini], \
	np.asarray(vsini_profiles[idx_vsini]) * stretch_factors[idx_vsini], \
	vsini_values, residuals_sum_all, calc_vsini_val, min_resid_val, stretch_fact
	



def fit_gaussian(obj_name, mask_type, data_x, data_y, mjd_obs):
	
	# function to fit Gaussian to the xcorr
	# data in order to compute RV and sigma 
	
	# define Gaussian function
	def gaussian(x, amp, cen, wid):
		return amp * np.exp(-(x-cen)**2 /wid)
		
	
	# first try and remove any residual slope in the CCF
	f = np.poly1d(np.polyfit(data_x, data_y, 1))
	lin_slope = f(data_x)
	
	# remove slope
	data_y = (data_y - lin_slope)
	
	# try and better set continuum to zero
	data_y = data_y - np.median(data_y)
	
	# initial guess parameters
	init_vals = [max(data_y), data_x[list(data_y).index(max(data_y))], 100] 
	
	# make linear array of x_vals
	x_lin = np.arange(min(data_x), max(data_x), 0.1)
	
	# do curve fit and get error
	# make sure the iteration isn't stopped
	# if there is a run time error
	try:
		best_vals, covar = curve_fit(gaussian, data_x, data_y, p0=init_vals)
	except RuntimeError:
		print("Error - curve_fit failed")
		# dummy values
		best_vals, covar = [0.5, 50, 50], [50]
	
	# get parameters out of 'best_vals' array
	amplitude, centroid, width = best_vals
	sigma = np.sqrt(width)
	
	# define error on amplitude
	perr = np.sqrt(np.diag(covar))[0] * amplitude
	
	
	# get depth of CCF
	idx_depth = getnearpos(x_lin, centroid)
	depth_val = gaussian(x_lin, amplitude, centroid, width)[idx_depth]
	
	# perform span analysis:
	# from continuum (0) to depth of CCF
	depth_array = np.arange(0, amplitude, amplitude / 50.)
	
	
	v1, v2 = centroid - (2 * sigma), centroid + (2 * sigma)
	
	# make two different arrays, one for each side of the CCF
	vel_span1 = np.arange(min([v1, v2]), centroid, 1)
	vel_span2 = np.arange(centroid, max([v1, v2]), 1)
	
	ccf_span1 = [data_y[getnearpos(v, data_x)] for v in vel_span1]
	ccf_span2 = [data_y[getnearpos(v, data_x)] for v in vel_span2]

	# make 'fine' arrays using linear interpolation
	vel_fine1 = np.arange(vel_span1[0], vel_span1[-2], 0.01)
	f1 = interpolate.interp1d(vel_span1, ccf_span1, 'linear')
	ccf_fine1 = f1(vel_fine1)

	vel_fine2 = np.arange(vel_span2[0], vel_span2[-2], 0.01)
	f2 = interpolate.interp1d(vel_span2, ccf_span2, 'linear')
	ccf_fine2 = f2(vel_fine2)

	
	bisector_array=[]
	
	for c_ in depth_array:
		# for each side get the corresponding velocity
		idx1 = getnearpos(ccf_fine1, c_)
		vel1 = vel_fine1[idx1]
		
		idx2 = getnearpos(ccf_fine2, c_)
		vel2 = vel_fine2[idx2]
		
		# get the mean value
		mean_vel = np.mean([vel1, vel2])
		
		#append to the bisector array
		bisector_array.append(mean_vel)
	
	
	# define upper and lower sectors
	l1, l2 = 0.10 * amplitude, 0.40 * amplitude
	u1, u2 = 0.55 * amplitude, 0.9 * amplitude
	
	lower_section, upper_section = [], []
	for c_, v_ in zip(depth_array, bisector_array):
		# lower boundary
		if c_ >= l1 and c_ <= l2:
			lower_section.append(v_)
			
		# upper boundary
		if c_ >= u1 and c_ <= u2:
			upper_section.append(v_)
	
	# get means (vt and vb) - vel. top and vel. bottom
	vt, vb = np.mean(upper_section), np.mean(lower_section)
	
	# for bisector slope
	b1, b2 = 0.25 * amplitude, 0.8 * amplitude
	
	bisector_slope_array, bisector_slope_depth = [], []
	for c_, v_ in zip(depth_array, bisector_array):	
		if c_ >= b1 and c_ <= b2:
			bisector_slope_array.append(v_)
			bisector_slope_depth.append(c_)
	
	if bisector_slope_array == []:
		bisector_slope_array = [999,999]
		bisector_slope_depth = [999,999]

	# do linear regression 
	inverse_m, intercept, rval, \
	pval, stderr = linregress(bisector_slope_depth, bisector_slope_array)
	
	
	# for 'curvature' measure (as defined in T. H. Dall et al. 2006)
	sec_array = [[0.2, 0.3],[0.4, 0.55],[0.75, 1]]
	
	curvature_array=[]
	for pair in sec_array:
		c_array=[]
		p1, p2 = pair[0], pair[1]
		for c_, v_ in zip(depth_array, bisector_array):
			if c_ >= (p1 * amplitude) and c_ <= (p2 * amplitude):
				c_array.append(v_)
		m = np.mean(c_array)
		curvature_array.append(m)
	
	z1, z2, z3 = curvature_array
	
	# calculate curvature measurement ((z3 - z2) - (z2 - z1))
	curvature_val = (z3 - z2) - (z2 - z1)
	

	# Anderson-Darling test
	
	a1, a2 = int(centroid - (1 * sigma)), int(centroid + (1 * sigma))
	vel_span_ad = np.arange(min([a1, a2]), max([a1, a2]), 1)
	ccf_span_ad = [data_y[getnearpos(v, data_x)] for v in vel_span_ad]
	
# 	plt.plot(vel_span_ad, ccf_span_ad)
# 	plt.show()

	ad_stat, crit_vals, sig_level = sp.stats.anderson(ccf_span_ad, dist='norm')
	
	# get significance level of AD stat
	# using crit_vals and sig_level
	
	sig_levels=[]
	for val, sig in zip(crit_vals, sig_level):
		if ad_stat / val < 1:
			sig_levels.append((sig / 100.))
	
	if sig_levels == []:
		sig_level = "$<$ 0.01"
		sig_level_num = 0
	else:
		sig_level = "%s" % (max(sig_levels))
		sig_level_num = max(sig_levels)
			
	

	
	# vsini part
	vel_shifted, ccf, vsini_vel, \
	vsini_profile, vsini_values, residuals_sum_all, \
	vsini_calc_val, min_chi, stretch_fact = get_vsini(data_x, data_y)
	
	
	
	# now make a plot of the results
	fig=plt.figure(num=None, figsize=(11, 11), dpi=80, facecolor='w', \
	edgecolor='k', tight_layout=True)
	
	# first panel
	ax=fig.add_subplot(2,2,1)
	plt.plot(data_x, data_y, '-', color='black')
	# plot Gaussian fit
	plt.plot(x_lin, gaussian(x_lin, amplitude, centroid, width), '-', color='blue', \
	label='RV: %.1f\nDepth: %.3f\nWidth: %.1f\nAD stat: %.3f (%s)\nMJD: %.2f' % (centroid, depth_val, \
																				np.sqrt(width), ad_stat, \
																				sig_level, mjd_obs))
	plt.xlabel("Velocity (km/s)",fontsize=12)
	plt.xlim(min(x_lin), max(x_lin))
	plt.ylabel("CCF",fontsize=12)
	plt.ylim(-0.01, max(depth_array)+0.01)
	plt.legend(frameon=False, loc='best', fontsize=12)
	plt.gca().invert_yaxis()
	
	# second panel
	ax=fig.add_subplot(2,2,2)
	
	# plot left and right hand side of CCF	
	plt.plot(vel_span1, ccf_span1, color='red')
	plt.plot(vel_span2, ccf_span2, color='blue')

	# plot bisector
	plt.errorbar(bisector_array, depth_array, fmt='o', mfc='black', mec='black', alpha=0.5, ms=2)
	plt.xlabel("Velocity (km/s)",fontsize=12)
	plt.ylabel("CCF",fontsize=12)
	plt.xlim(min(vel_span1), max(vel_span2))
	plt.ylim(-0.01, max(depth_array)+0.01)
	plt.gca().invert_yaxis()
	
	subpos = [0.6,0.6,0.5,0.5]
	ax3=fig.add_subplot(2,2,3)
	
	plt.xlim(min(vel_span1), max(vel_span2))
	plt.plot(data_x, ccf, color='black', lw=1)
	plt.plot(np.asarray(vsini_vel)+centroid, vsini_profile, '-', color='red', lw=1.5, \
	label='Best fit: %i km/s\nStretch: %.1f' % (vsini_calc_val, stretch_fact))	
	plt.ylim(-0.8, 1)
	plt.gca().invert_yaxis()
	plt.xlabel("Velocity (km/s)",fontsize=12)
	plt.ylabel("Relative depth",fontsize=12)
	plt.legend(frameon=False, loc='lower left', fontsize=12)
	
	subax1 = add_subplot_axes(ax3,subpos)
	plt.plot(vsini_values, residuals_sum_all, color='blue', lw=1.5)
	
	idx_res = getnearpos(vsini_values, (2 * vsini_calc_val))
	res_lim = residuals_sum_all[idx_res]
	
	plt.xlabel("vsin(i) value",fontsize=12)
	plt.ylabel("Profile residual",fontsize=12)
	plt.xlim(0, (2 * vsini_calc_val))
	plt.ylim(0, res_lim)
	
	
	
	
	ax=fig.add_subplot(2,2,4)
		
	# plot boundaries for upper and lower sectors
	# for bisector
	plt.plot([min(vel_span1), max(vel_span2)], [l1, l1], '--', color='grey', lw=0.8)
	plt.plot([min(vel_span1), max(vel_span2)], [l2, l2], '--', color='grey', lw=0.8)
	plt.plot([min(vel_span1), max(vel_span2)], [u1, u1], '--', color='grey', lw=0.8)
	plt.plot([min(vel_span1), max(vel_span2)], [u2, u2], '--', color='grey', lw=0.8)
	
	# fill areas for clarity on plot
	
	plt.fill_between([min(vel_span1), max(vel_span2)], l1, l2, \
	facecolor='SkyBlue', interpolate=True, alpha=0.2)
	
	plt.fill_between([min(vel_span1), max(vel_span2)], u1, u2, \
	facecolor='SkyBlue', interpolate=True, alpha=0.2)	
	
	# plot text of bisector span, slope and curvature
	spacing = max(depth_array) * 0.0625
	
	plt.text(min(bisector_array)+0.1, max(depth_array) * 0.7, \
	r"$\bar{v}_t - \bar{v}_b$ = %.3f" % (vt-vb), fontsize=12)
	
	plt.text(min(bisector_array)+0.1, max(depth_array) * 0.7 + spacing, \
	r"$b_b$ = %.3f" % (inverse_m))

	plt.text(min(bisector_array)+0.1, max(depth_array) * 0.7 + (2* spacing), \
	r"$c_b$ = %.3f" % curvature_val)
		
	
	fit_array = np.asarray([b1, b2])
	plt.plot((fit_array * inverse_m)+intercept, fit_array, '-', color='red', lw=1.5)
	plt.errorbar(bisector_array, depth_array, fmt='o', mfc='black', mec='black', alpha=0.5, ms=4)
	plt.ylim(-0.01, max(depth_array)+0.01)
	plt.ylabel("CCF",fontsize=12)
	plt.gca().invert_yaxis()
	
	plt.xlim(min(bisector_array)-0.2, max(bisector_array)+0.5)
	#plt.legend(frameon=False, loc='best', fontsize=12)
	plt.xlabel("Velocity (km/s)",fontsize=12)
	

	# plot output
	# plot_output = "./OUTPUTS/PDFS/ccf_%s_%s.pdf" % (obj_name, mask_type)
	plot_output = "./OUTPUTS/PDFS/ccf_%s_%s_%s.pdf" % (obj_name, mask_type, str(mjd_obs))
	plt.savefig(plot_output)
	print "Made a plot at: %s..." % plot_output
	
	plt.close()
	
	# print RV and sigma
	print "\n**************\nRV: %.1f" % centroid
	print "Sigma: %.2f\n**************" % sigma
	
	# define radial velocity and sigma
	rv = centroid

	
	return rv, sigma, depth_val, vsini_calc_val, min_chi, \
		   (vt-vb), inverse_m, stderr, curvature_val, sig_level_num
	

# now define one function that calls all the other functions

def master_make_ccf_vsini(file_location, flag='fits', mask_type='K0', instrument='UVES'):

	print "\nStarting CCF and Vsin(i) procedure..."

	# first read in file and get spectrum
	if flag == 'fits':
		if instrument == 'UVES_PH3':
			obj_name, ra_deg, dec_deg, mjd_obs, \
			wlen_obs, flux_obs, bary_correc = ReadSpec_UVES_PH3_fits(file_location)
		if instrument == 'UVES':
			obj_name, ra_deg, dec_deg, mjd_obs, \
			wlen_obs, flux_obs, bary_correc = ReadSpec_UVES_fits(file_location)
		if instrument == 'UVES_combine':
			obj_name, ra_deg, dec_deg, mjd_obs, \
			wlen_obs, flux_obs, bary_correc = combine_spec(file_location[0], file_location[1])
		if instrument == 'FEROS':
			obj_name, ra_deg, dec_deg, mjd_obs, \
			wlen_obs, flux_obs, bary_correc = ReadSpec_FEROS_fits(file_location)			
		if instrument == 'HARPS':	
			obj_name, ra_deg, dec_deg, mjd_obs, \
			wlen_obs, flux_obs, bary_correc = ReadSpec_HARPS_fits(file_location)			
			
						
	if flag == 'txt':
		obj_name, wlen_obs, flux_obs, bary_correc = ReadSpec_txt(file_location)	
		
	print "File: %s" % file_location

	print "Got info for object: %s..." % obj_name
	print "Barycentric correction is %.3f km/s" % bary_correc

	# now get mask (default value is 'M4')
	wlen_mask, flux_mask = get_mask(mask_type)
	
	print "Got Mask of type %s..." % mask_type
	
	# now make ccf profile
	print "Will now perform CCF fit for velocities..."
	velocity, ccf_profile, snr_median, central_wlen, wlen_range = \
	make_ccf(wlen_obs, flux_obs, wlen_mask, \
	         flux_mask, bary_correc, wlen_unit='A', plot_flag=False)
	
	print "\nNow going to fit CCF with a Gaussian..."
	
	# fit Gaussian to CCF
	rv, width, depth, vsini_calc_val, min_chi, \
	BIS, b_b, stderr, c_b, ad_sig = fit_gaussian(obj_name, mask_type, velocity, ccf_profile, mjd_obs)


	print "\nFinally I will fit the CCF with vsini rotational profiles..."

	
	
	print "************************************"
	print "Summary of results:"
#	print "\nRV: %.1f, sigma: %.1f, vsini: %i km/s" % (rv, width, vsini_calc_val)
	print "************************************\n"
	
	return obj_name, ra_deg, dec_deg, mjd_obs, rv, width, depth, vsini_calc_val, \
	min_chi, snr_median, central_wlen, wlen_range, BIS, b_b, stderr, c_b, ad_sig


