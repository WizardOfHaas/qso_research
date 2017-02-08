import requests
import matplotlib.pyplot as plt
import numpy as np
import sys

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

	plateid = int(spec[0])
	mjd = int(spec[1])
	fiber = int(spec[2])
	z = float(spec[3])
	
	wavelength_CIV_shifted = (z + 1) * (wavelength_CIV_emit) #Shift CIV
	wavelength_NV_shifted = (z + 1) * (wavelength_NV_emit) #Shift NV
	
	#Send HTTP request to skyserve to get csv spectra
	r = requests.get("http://dr9.mirror.sdss3.org/csvSpectrum?mjd=" + `mjd` + "&amp;fiber=" + `fiber` + "&amp;plateid=" + `plateid` + "&amp;reduction2d=v5_4_45")
	
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
	
	#Some speed shit
	beta = []
	avr_CIV_doublet = 1549. #Average wavelength of CIV doublet
	z_absC = (np.asarray(wavelength)/avr_CIV_doublet)-1.
	RC = (1.+z)/(1.+z_absC)
	betaC = ((RC**2.)-1.)/((RC**2.)+1.)
	betaa = -betaC*(300000.)
	for ll in betaa:
		betas = round (ll,4)
		beta.append (betas)
	beta = np.asarray(beta)
	
	#Normalize it!
	flux = flux / fit
	
	#Smoothe it!
	flux = smooth(flux, 3)
	fit = smooth(fit, 3)
	
	#Calculate BI
	bracs = 1 - (flux / 0.9) #Calculate inside brackets
	troughs = np.where(bracs > 0) #Where is there absorption?
	troughs = group_consecutives(troughs[0]) #Group it into descreet instances
	
	bal_cutoff = 1000
	bis = []
	bals = []
	
	for trough in troughs:
		betas = beta[trough]
		width = abs(betas[0] - betas[-1])
	
		if width > bal_cutoff:
			bals.append(trough)
			bi = 0
	
			for i in trough:
				if len(beta) > i + 1:
					bi += bracs[i] * abs(beta[i] - beta[i + 1])
	
			bis.append(bi)
	
	total_bi = sum(bis)
	
	#Plot curves
	plt.plot(beta, flux)
	for trough in bals:
		plt.plot(beta[trough], flux[trough], color = "red")
	
	#plt.plot(wavelength, fit)
	#plt.plot(wavelength, depths)
	#plt.plot(wavelength, error)
	
	#Plot markers
	#plt.axvline(x = wavelength_CIV_shifted)
	#plt.axvline(x = wavelength_NV_shifted)
	plt.axhline(y = 1)
	
	plt.title(`plateid` + "-" + `mjd` + "-" + `fiber` + ",z=" + `z` + ",BI=" + `round(total_bi, 2)`)

	#Show it!
	#plt.show()
	plt.savefig("plots/" + `plateid` + "-" + `mjd` + "-" + `fiber` + ".png")
	plt.clf()
	
	print `plateid` + "," + `mjd` + "," + `fiber` + "," + `z` + "," + `total_bi`