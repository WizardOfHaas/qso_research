import requests
import matplotlib.pyplot as plt
import numpy as np
import sys
import json

############################################Functional Definitions
def smooth(norm_flux, box_pts):       
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(norm_flux, box, mode='same')
    return y_smooth

def group_consecutives(vals, step=1):	
    result = [[vals[0]]]
    group = 0

    for i in range(1, len(vals)):
    	if vals[i +-1] == vals[i] - 1:
    		result[group].append(vals[i])
    	else:
    		group += 1
    		result.append([vals[i]])
    
    return result
##################################################################

#Set some constants... for math
wavelength_CIV_emit = 1550.7700
wavelength_NV_emit = 1242.8040

#load config file
config_file = sys.argv[1]
f = open(config_file, 'r')
config = f.readlines()

for spec in config:	
	#Set search params
	spec = spec.split(',')

	plateid = spec[0].zfill(4)
	mjd = spec[1].zfill(5)
	fiber = spec[2].zfill(4)
	z = float(spec[3])
	
	wavelength_CIV_shifted = (z + 1) * (wavelength_CIV_emit) #Shift CIV
	wavelength_NV_shifted = (z + 1) * (wavelength_NV_emit) #Shift NV
	
	#Send HTTP request to skyserve to get csv spectra
	#print "http://api.sdss3.org/spectrum?mjd=" + mjd + "&fiber=" + fiber + "&plateid=" + plateid + "&type=json"
	r = requests.get("http://api.sdss3.org/spectrum?mjd=" + mjd + "&fiber=" + fiber + "&plate=" + plateid + "&format=json")
	
	data = json.loads(r.content)
	wavelength = np.asarray(data['wavelengths'])
	flux = np.asarray(data['flux'])
	fit = np.asarray(data['best_fit'])
	error = np.asarray(data['wavelength_dispersion'])
	sky_flux = np.asarray(data['sky_flux'])

	"""
	lines = r.content.split('\r\n') #Split by line
	lines.pop(0) #Remove header
	
	data = [] #Place to store spectra
	
	wavelength = []
	flux = []
	fit = []
	error = []
	
	#Dirty csv parse, ok since the data is in a known format
	for line in lines:
		fields = line.split(',')
		if len(fields) == 4:
			fields = map(float, fields) #Make it usable
	
			#Filter for the region of interest
			if fields[0] >= wavelength_NV_shifted and fields[0] <= wavelength_CIV_shifted:
				data.append(fields) #Push full object
				wavelength.append(fields[0]) #Push wavelengths
				flux.append(fields[1]) #Push flux
				fit.append(fields[2]) #Push bestfit
				error.append(fields[3]) #Push error
	
	#Make everything a damned array...
	wavelength = np.asarray(wavelength)
	flux = np.asarray(flux)
	fit = np.asarray(fit)
	error = np.asarray(error)
	"""

	wavelength = wavelength / (z + 1) #Shift to rest-frame

	#Normalize it!
	flux = flux / fit
	error = error / fit
	erorr = error ** (1. / 3)

	#Convert! (bad pizels)
	flux_final = flux
	error_final = error
	for i in range (1, len (flux) - 1):
		if error[i] / flux[i] > 1:
			#flux_final[i] = (flux[i - 1] + flux[i + 1]) / 2
			#error_final[i] = (error[i - 1] + error[i + 1]) / 2

			flux_final[i] = flux[i - 1]
			error_final[i] = error[i - 1]

	"""
	for i in range (1, len (flux) - 1):
		if abs(flux[i + 1] - flux[i]) > 0.5 and error[i + 1] > 0.25:
			flux_final[i + 1] = flux[i]
			error_final[i + 1] = error[i]

		if error[i] > 0.5:
			flux_final[i] = flux[i - 1]
			error_final[i] = error[i - 1]
	"""

	flux = flux_final
	error = error_final

	#Smoothe it!
	flux = smooth(flux, 3)

	data = (wavelength, flux, error)
	data =(np.transpose(data))

	np.savetxt('/home/sean/qso_data/data/dr9_flux/sdss_norm/spec-' + plateid + "-" + mjd + "-" + fiber + '.norm.DR9', data)

	plt.plot(wavelength, flux)
	plt.plot(wavelength, error, color = "black")
	plt.axhline(y = 1)
	plt.title(plateid + "-" + mjd + "-" + fiber + ",z=" + `z`)
	plt.savefig("plots/" + plateid + "-" + mjd + "-" + fiber + ".pdf")
	plt.clf()

	
	print `plateid` + "," + `mjd` + "," + `fiber` + "," + `z`