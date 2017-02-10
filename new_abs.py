import requests
import matplotlib.pyplot as plt
import numpy as np
import sys
from pylab import *
from sympy import sympify
from scipy import *
from astropy import *

############################################Functional Definitions
def smooth(norm_flux, box_pts):       
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(norm_flux, box, mode='same')
    return y_smooth

def group_consecutives(vals, step=1):	
    result = [[vals[0]]]
    group = 0

    for i in range(1, len(vals)):
    	if (vals[i - 1] == vals[i] - 1) or (vals[i - 2] == vals[i] - 2) or (vals[i - 3] == vals[i] - 3):
    		result[group].append(vals[i])
    	else:
    		group += 1
    		result.append([vals[i]])
    
    return result
##################################################################

#Set some constants... for math
wavelength_CIV_emit = 1550.7700
wavelength_NV_emit = 1242.8040

#Velocity cutoffs
maxvel= 0.#original
minvel= -60000#original

#load config file
config_file = sys.argv[1]
f = open(config_file, 'r')
config = f.readlines()

for spec in config:
	#Set search params
	spec = spec.split(',')

	specname = spec[0]
	z = float(spec[1])
	
	wavelength_CIV_shifted = (z + 1) * (wavelength_CIV_emit) #Shift CIV
	wavelength_NV_shifted = (z + 1) * (wavelength_NV_emit) #Shift NV
	
	#Load spectral file
	norm_spec = loadtxt('/home/sean/qso_data/data/dr9_flux/norm/' + specname + 'norm.DR9')

	#Make everything a damned array...
	wavelength = np.asarray(norm_spec[:,0])
	flux = np.asarray(norm_spec[:,1])
	error = np.asarray(norm_spec[:,2])
	
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
	
	#Smoothe it!
	flux = smooth(flux, 3)
	
	#Calculate BI
	bracs = 1 - (flux / 0.9) #Calculate inside brackets
	troughs = np.where(bracs > 0) #Where is there absorption?
	troughs = group_consecutives(troughs[0]) #Group it into descreet instances
	#troughs = [x for x in troughs if len(x) > 3]

	bal_cutoff = 1000
	bis = []
	ews = []
	bals = []
	measured_bals = []
	vlims = []
	ion_lines = []
		
	#Draw line for possible SiIV...
	#Calculate EW & BI...
	for trough in troughs:
		betas = np.sort(beta[trough])
		width = abs(betas[0] - betas[-1])			

		if (betas[0] > minvel) and (betas[-1] < maxvel):			
			bi = 0
			ew = 0
			bi_width = 0
			bal = []

			for i in reversed(trough):
				if len(beta) > i + 1:
					dv = abs(beta[i] - beta[i + 1])
					bi_width += dv
					if bi_width > bal_cutoff:
						bi += bracs[i] * dv
						bal.append(i)

			if len(bal) > 0:
				for i in bal:
					ew += (1 - flux[i]) * abs(beta[i] - beta[i + 1])

			#Find location(?) of possible peaks
			SiIV_emitted = 1397.#this line and the four lines below, I'm just trying to find the corresponding wavelenth value of CIV absorption feature
			z_absSiIV = (wavelength[i]/SiIV_emitted)-1#
			obs_wavelength_C=(z_absSiIV+1)*(1549.4825)#
			obs_wavelength_C_index =np.min (where (wavelength>obs_wavelength_C))#
                    
			obs_wavelength_C_final=beta[obs_wavelength_C_index]+600
			ion_lines.append(obs_wavelength_C_final)

			if len(bal) > 0:
				bals.append(trough) #Add to list of BALs
				vlims.append([betas[-1], betas[0]]) #Add to list of vmin/vmax
				measured_bals.append(bal) #Add to list of measured BALs
				bis.append(bi) #Append BI
				ews.append(ew) #Append EW

	#Get total BI	
	total_bi = sum(bis)
	total_ew = sum(ews)
	
	#Get vmin/vmax
	vlims = np.asarray(vlims)
	vlims_flat = vlims.flatten()

	if len(vlims_flat) > 0:
		vmax = np.max(vlims_flat)
		vmin = np.min(vlims_flat)
	else:
		vmax = 0
		vmin = 0

	#Scale the error
	final_norm_error = error / 3**(0.5)

	#Plot curves
	fig, ax = plt.subplots()
	plt.plot(beta, flux, color = "green") #Plot flux

	for i, trough in enumerate(bals): #Plot absorption features/ew lines
		plt.plot(beta[trough], flux[trough], color = "red")

		#plt.plot(
		#		[
		#			beta[trough[0]],
		#			beta[trough[0]] + ews[i]
		#		],[
		#			1.2,
		#			1.2
		#		],
		#		color = "black")

		ax.axvspan(beta[trough[0]], beta[trough[0]] + ews[i], 0, 0.3)

		#plt.axvline(ion_lines[i], color = "black")
		plt.axvline(beta[trough][-1], color = "red")	

	for bal in measured_bals:
		plt.plot(
			[
				beta[bal[0]],
				beta[bal[-1]]
			],[
				0.95,
				0.95
			],
			color = "black")

	plt.plot(beta, final_norm_error, 'k-') #Plot error

	#Plot markers	
	plt.axhline(y = 1) #Plot guideline
	plt.axhline(y = 0.9, ls = "dashed") #BI line

	plt.axvline(vmax, color = "black", ymax = 0.1) #Plot vmax
	plt.axvline(vmin, color = "black", ymax = 0.1) #Plot vmin
	ax.axvspan(-75000, vmin, alpha=0.5, color='grey')
	ax.axvspan(vmax, 0, alpha=0.5, color='grey')
	
	plt.text(-70000, 2.5, 
		"vmin = " + `round(vmin, 2)` + "\n" +
		"vmax = " + `round(vmax,2)` + "\n" +
		"BI = " + `round(total_bi, 2)` + "\n" +
		"EW = " + `round(total_ew, 2)`)
	
	#Set plot range
	plt.xlim(-75000, 0)
	plt.ylim(0, 3)

	plt.title(`specname` + ",z=" + `z`)

	#Show it!
	plt.savefig("plots/" + specname + ".pdf")
	plt.clf()
	plt.close(fig)
	
	print specname + "," + `z` + "," + `total_bi` + "," + `vmin` + "," + `vmax` + "," + `total_ew`
	#print ews
	#print bis
	#exit()