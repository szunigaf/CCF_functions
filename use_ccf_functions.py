import sys
sys.path.append("./MODULES")

from calc_rv_new_chapman import *
import glob


# cp calc_rv_new.py calc_rv_new.pyx
# setup.py (contains calc_rv_new module)
# then:
# python setup.py build_ext --inplace

# define data directory
#data_direc = './example_XSHOOTER_PH3/'
#data_direc = './example_UVES_combine/'
#data_direc = './example_UVES_PH3/'
data_direc = './INPUT/'

# define the mask
#mask_name='R50G2'
#mask_name = 'M4'
mask_name = 'K0'


# make a list of all the file names
data_files = glob.glob(data_direc+"*fits")
#data_files = glob.glob(data_direc+"*ascii")

#data_files =["r.FEROS.2012-11-26T06:48:05.136.1081.fits", ]

# make a list of all the file names without the directory
data_files_names = [f.replace(data_direc, "") for f in data_files]
#data_files_names = ['Sz91.ascii']
print (data_files_names)
# if running again
# *******************************


# define master output for properties
master_output = "./OUTPUTS/master_output.txt"

existing_data = np.loadtxt(master_output, dtype={'names':('file','obj_name','ra_deg','dec_deg',\
                                                          'mjd_obs','rv','width','depth','vsini_calc_val',\
                                                          'min_chi','snr_median','mask','central_wlen',\
                                                          'wlen_range','BIS','b_b', 'b_b_err','c_b','ad_sig'), \
                                                          'formats':('S100', 'S50', 'f8', 'f8', 'f8', 'f8', \
                                                          'f8', 'f8', 'f8', 'f8', 'f8', 'S2', 'f8', \
                                                          'f8', 'f8', 'f8', 'f8', 'f8', 'f8')}, \
                                                          skiprows=1, unpack=True)

# get the list of files in the table
#commented out because there are no result since the first run of the code
#existing_files = existing_data[0]

# make a new list of all the files that aren't in common
# between data_files_names and existing_files

#files_without_ccf = []
#for file_name in data_files_names:
#       if file_name not in existing_files:
#               files_without_ccf.append(data_direc+file_name)

#First iteration where there is nothing in the master output
files_without_ccf = data_files

# print (something that tells the user how many)
# files are still left to process
print ("\n\n%i / %i still not in output table\n\n" % \
(len(files_without_ccf), len(data_files_names)))

#instrument='UVES_combine'
#instrument='XSHOO_PH3'
#instrument='UVES'
#instrument='UVES_PH3'
#instrument='FEROS'
instrument='HARPS'


i=0
while i < len(files_without_ccf):
#for file in files_without_ccf:
#        print (file)
        if instrument != 'UVES_combine':
                obj_name, ra_deg, dec_deg, \
                mjd_obs, rv, width, depth, \
                vsini_calc_val, min_chi, \
                snr_median, central_wlen, \
                wlen_range, BIS, b_b, stderr, \
                c_b, ad_sig = master_make_ccf_vsini(files_without_ccf[i], flag='fits', mask_type=mask_name, instrument=instrument)
                output_string = "%s     %s      %.5f    %.5f    %.3f    %.2f    %.2f    %.3f    %i      %.3f    %.2f    %s      %.1f    %.1f    %.3f    %.3f    %.3f    %.3f    %.2f\n" % (files_without_ccf[i].replace(data_direc, ""), obj_name, ra_deg, dec_deg, \
                mjd_obs, rv, width, depth, vsini_calc_val, min_chi, snr_median, mask_name, central_wlen, wlen_range, BIS, b_b, stderr, c_b, ad_sig)
                text_file = open(master_output, "a")
                text_file.write(output_string)
                text_file.close()
                i+=1        
        else:
                obj_name, ra_deg, dec_deg, \
                mjd_obs, rv, width, depth, \
                vsini_calc_val, min_chi, \
                snr_median, central_wlen, \
                wlen_range, BIS, b_b, stderr, \
                c_b, ad_sig = master_make_ccf_vsini([files_without_ccf[i], files_without_ccf[i+1]], flag='fits', mask_type=mask_name, instrument=instrument)
                output_string = "%s     %s      %.5f    %.5f    %.3f    %.2f    %.2f    %.3f    %i      %.3f    %.2f    %s      %.1f    %.1f    %.3f    %.3f    %.3f    %.3f    %.2f\n" % (files_without_ccf[i].replace(data_direc, "")+'combined', obj_name, ra_deg, dec_deg, \
                mjd_obs, rv, width, depth, vsini_calc_val, min_chi, snr_median, mask_name, central_wlen, wlen_range, BIS, b_b, stderr, c_b, ad_sig)
                text_file = open(master_output, "a")
                text_file.write(output_string)
                text_file.close()
                i+=2        
