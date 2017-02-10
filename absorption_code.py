import sys

from pylab import*
from numpy import*
from sympy import sympify
from scipy import*
from astropy import*
from matplotlib.pyplot import*
from scipy.optimize import curve_fit
from matplotlib.backends.backend_pdf import PdfPages

pp=PdfPages('absorption_dr10_dr10.pdf')
fig=figure()

#redshifts = loadtxt('confs/DR9_norm_zem.lst')
#redshifts_action = redshifts[:]

#spectra = loadtxt('confs/DR9_norm_specnames.lst',dtype='str')
#spectra_action = spectra[:]

#######USE AN ACTUAL CONFIG FILE!!!!!
config_file = sys.argv[1] #Set cfg file path (csv with spec_name,z,snr...)
config = loadtxt(config_file,
    delimiter = ",",
    dtype={
        'names': ('spec', 'z', 'snr'),
        'formats': ('|S35', np.float, np.float)
    })

spectra_action = [row[0] for row in config]
redshifts_action = [row[1] for row in config]
snr_action = [row[2] for row in config]

brac_all=[]
deltav_all=[]
thousand=[]

wavelength_CIV_emit1=1550.7700
wavelength_CIV_emit2=1548.1950
count=0
vmins=[]
vmins_all=[]
vmaxs_all=[]
vmaxs=[]
BI_all=[]
loop1=[]
BI=0
BI_adding=[]
thousand_final=[]
mnim_i=[]
mnim_j=[]
all_final_depth=[]
all_all_final_depth=[]
BI_individual=[]
BI_individual_appended=[]
BI_all_individual=[]
BI_for_individual=[]

EW_individual=[]
EW_ind=[]
EW_all_individual=[]
vlast=[]







for i,j,snr in zip (spectra_action, redshifts_action, snr_action): #Iterate over lists of sdss_ids and redshifts
    snr = round(snr, 5)

    count=count+1
    print count
    
    
    BI_individual_appended=[]
    EW_individual=[]
    beta=[]
    
    
    
    

    


    
    
        
    
        
    
    normalized_dr9 =  loadtxt('/home/sean/qso_data/data/dr9_flux/norm/'+i+'norm.DR9') #Load in normalized spectrum
    thousand=0
    
    wavelength = normalized_dr9[:,0] ###################wavelength
    norm_flux = normalized_dr9[:,1] ####################normalized flux
    norm_error = normalized_dr9[:,2] ###################normalized error
    BI_mid=[]
    #################################################Defining smoothing function
    
    def smooth(norm_flux, box_pts):
        
        box = np.ones(box_pts)/box_pts
        y_smooth = convolve(norm_flux, box, mode='same')
        return y_smooth
        

    
    
    z=j #Rename redshift
    z=round (z,5) #Round the redshift
    sm_flux=smooth(norm_flux,3) #Smooth the spectrum (3 point boxcar)
#############################################New (not using in the code)
    avr_SiIV_doublet = 1397.
    z_absS = (wavelength/avr_SiIV_doublet)-1
    obs_wavelength_C=(z_absS+1)*(1549.4825)
    count2=0
    
    

    
    
    RS=(1+z)/(1+z_absS)
    betaS=((RS**2)-1)/((RS**2)+1)
    beta1=-betaS*(300000)
    
    wavelength_CIV_abs1=(z+1)*(wavelength_CIV_emit1)
    wavelength_CIV_abs2=(z+1)*(wavelength_CIV_emit2)
 ########################################################New  (not using in the code) 
                               

    #################################What paola had in her code
    avr_CIV_doublet = 1549. #Average wavelength of CIV doublet
    z_absC = (wavelength/avr_CIV_doublet)-1.
    RC=(1.+z)/(1.+z_absC)
    betaC=((RC**2.)-1.)/((RC**2.)+1.)
    betaa = -betaC*(300000.)
    for ll in betaa:
        betas=round (ll,4)
        beta.append (betas)
    beta=array(beta)
    
    ################################What paola had in her code

    
    #maxvel and minvel are the limits, in velocity, where we want to find absorption and assume that they are SiIV
    maxvel= -25000.#original
    minvel=-60000#original

    fst = 0

    if beta.any():
        try:
            fst = np.max(where (beta <= maxvel))#index value of the starting point (on the very left) -- index value of minvel
        except:
            #fst = np.max(where (beta == maxvel))
            fst = 0

    try:
        lst = np.min(where(beta >=minvel))#index value of the ending point (on the very right) -- index value of maxvel
    except:
        lst = where (beta == np.min(beta))
    
    deltav=0#change in velocity
    part=0
    bb=-1
    
    countBI=1000 #this is the lower limit of how we are categorizing BALs, so if the width of an absorption feature is above 600, we are considering it to be a BAL
    
    jjj=arange (lst,fst)
    jjj=array(jjj)
    jjj = jjj[::-1]
    
    figure(count) #everything indented below is for the graph 
    
    vmins=[]
    vmaxs=[]
    all_final_depth=[]
    
    for jjjs in jjj:
        
        
        
        
            
        
        C=0
        
        
        #beta[jjj] is the range of velocity (velocity=x-axis) where we want to find BAL
        #              [1 - f(v)/0.9]
        brac=(1. - (norm_flux[jjjs]/0.9))# brac is 1 when norm_flux is negative (or when it is an absorption feature)
         
        if brac >0: #if brac >0, we have an absorption feature
            deltav=beta[jjjs+1]-beta[jjjs]  # TEST
            part=part+deltav#initially, part=0, for every positive brac, part is the last value of part+deltav, so when that adds up to 600, we have got our absorption
            bb=0 #forget about bb
            brac_all.append(brac)
            deltav_all.append(deltav)

            
            D=1# this plays the role of C but for equivalent width
            EW_individuals = (brac*D)*(deltav)
            EW_individuals = round (EW_individuals,4)
            EW_ind.append (EW_individuals)

            
                        
            if part >=countBI and bb==0:#code will go through this only when hte absorption is wider than 600km/s.
                index_depth_final=[]
                flux_depth=[]
                
                
                    
                
                
            
                
                C=1
                
                if beta[jjjs] < -3000 and beta[jjjs] >-30000:
                    BI=(brac*C)*(deltav)
                    BI=round (BI,4)
                    BI_mid.append(BI)

                    

                # the following two lines are for individual BI
                BI_for_individuals = (brac*C)*(deltav)
                BI_for_individuals=round (BI_for_individuals,4)
                BI_for_individual.append (BI_for_individuals)

                
                
                
                
                    
                    
                
                plot((beta[jjjs+1],beta[jjjs]),(1.5,1.5),'k-')# i changed this line
                if count2==0:  # i have made count2 initially zero
                    thousand = thousand+1
                    if thousand == 1:
                        thousand_final.append(thousand)
                        mnim_i.append(i)
                        mnim_j.append(j)
                

                    #plot((beta[jjjs+1],beta[jjjs]),(1.5,1.5),'k-')#you could just beta[jjjs] but veclocity (and wavelength) are not defined at all points, so to have a line, you need to do this, nothing to do with integral
                #plot(beta,norm_flux)
                    
                    
                ####################################################################################################code counts from the left
                
                    
                    SiIV_emitted = 1397.#this line and the four lines below, I'm just trying to find the corresponding wavelenth value of CIV absorption feature
                    z_absSiIV = (wavelength[jjjs]/SiIV_emitted)-1#
                    obs_wavelength_C=(z_absSiIV+1)*(1549.4825)#
                    obs_wavelength_C_index =np.min (where (wavelength>obs_wavelength_C))#
                    
                    obs_wavelength_C_final=beta[obs_wavelength_C_index]+600
                    plot((obs_wavelength_C_final, obs_wavelength_C_final),(-1,2),'k-')#I have added this line (11/08/2016)

                    
                    vvmins=np.min(where(beta > beta[jjjs]+600))#the logic of this line is: the value of speed the code is at right now is obviously 600km/s to the right of the vmin of the absorption feature (this is the index value)
                    
                    plot ((beta[vvmins], beta[vvmins]) , (-1,2),'r-')#NEW LINE: PLOTTING RED LINE WHERE THE ABSORPTION STARTS, NOT WHERE 600KM/S STARTS
                    vvvmins=beta[vvmins]
                    vvvmins=round (vvvmins,4)
                    vmins.append(vvvmins)
                if (1. - (norm_flux[jjjs-1]/0.9)) <0: # logic of this line is: if the next value of brac is negative (meaning the absorption feature has ended) then the current point of jjjs is our Vmax.
                    vvmaxs =beta[jjjs]
                    vvmaxs=round (vvmaxs,4)
                    vmaxs.append(vvmaxs)
                    
                    BI_individual = sum (BI_for_individual)
                    BI_individual=round(BI_individual,4)
                    BI_individual_appended.append (BI_individual)# this array contains one single BI value of each absortopn feature in a single spectrum
                    BI_for_individual=[]
                    
                    EW_ind_sum = sum (EW_ind)
                    EW_ind_sum=round(EW_ind_sum,4)
                    EW_individual.append (EW_ind_sum)
                    EW_ind=[]

                    
                    temp_index_vmin = where (beta ==beta[vvmins])
                    temp_index_vmax = where (beta ==vvmaxs )
                    

                    depth_array= beta [temp_index_vmin[0]:temp_index_vmax[0]]#original
                    
                    if len (depth_array) <1:
                        depth_array= beta [temp_index_vmax[0]:temp_index_vmin[0]]
                        
                    
                    for depth in depth_array:
                        index_depth = where (depth == beta)
                        index_depth_final.append(index_depth)
                        
                    for flux_depths in index_depth_final:
                        flux_depth1 = sm_flux[flux_depths]
                        flux_depth.append (flux_depth1)
                    final_depth = 1 - (np.min (flux_depth))
                    final_depth=round (final_depth,4)
                    all_final_depth.append (final_depth)
                    


                    
                if beta [jjjs] == beta[lst]: #the logic of this line is: if it happens that the absorption feature that the code is currently at continues even after our range (range of where we want to find the absorption feature), then just take the last point of our range as the Vmax
                    vvmaxs=(beta[jjjs])
                    vvmaxs=round (vvmaxs,4)
                    vmaxs.append(vvmaxs)

                    BI_individual = sum (BI_for_individual)
                    BI_individual=round(BI_individual,4)
                    BI_individual_appended.append (BI_individual)# this array contains one single BI value of each absortopn feature in a single spectrum
                    BI_for_individual=[]
                    
                    EW_ind_sum = sum (EW_ind)
                    EW_ind_sum=round(EW_ind_sum,4)
                    EW_individual.append (EW_ind_sum)
                    EW_ind=[]

                    
                    temp_index_vmin = where (beta==beta[vvmins ])
                    temp_index_vmax = where (beta ==vvmaxs )
                    

                    depth_array= beta [temp_index_vmin[0]:temp_index_vmax[0]]#original
                    
                    if len (depth_array) <1:
                        depth_array= beta [temp_index_vmax[0]:temp_index_vmin[0]]
                        
                    
                    for depth in depth_array:
                        index_depth = where (depth == beta)
                        index_depth_final.append(index_depth)
                        
                    for flux_depths in index_depth_final:
                        flux_depth1 = sm_flux[flux_depths]
                        flux_depth.append (flux_depth1)
                    final_depth = 1 - (np.min (flux_depth))
                    final_depth=round (final_depth,4)
                    all_final_depth.append (final_depth)

                
                   
                    
                        
                count2=1 #(original) logic of this line: now the code will go through the if satement on line 205 only when the code goes over a new absorption feature
                

            
        
        else: #if the brac value is not more than zero (so if we don't have absorption feature)
            part=0 # this is so b/c we do not want to keep counting the width of the absorption feature if it is not wider than 600km/s
            count2=0# this is so b/c if the ocdee encounters an other absorption feature which is wider than 600km/s, the code is going to go thorugh the if statement on line 205
            EW_ind=[]
            

        
    
        
        if jjjs == lst:
            BI_adding= sum(BI_mid)
            BI_adding=round(BI_adding,4)
            
            BI_all.append(BI_adding)#BI_all considers all absorption features where brac >0.


            
            BI_all_individual.append (BI_individual_appended)# contains multiple arrays, each array contains BI values for each absorption feature of a single spectrum
            EW_all_individual.append (EW_individual)
        
        



            
    
    yes=('name: ' + `i` + '\n' + 'BI (-3000 > v > -30,000): ' + `BI_adding` + '\n' +  'vmins: ' + `vmins` + '\n' + 'vmaxs: '+`vmaxs`+ '\n' + 'BI_individual: '+`BI_individual_appended`+ '\n' + 'EW_individual: '+`EW_individual`+ '\n' + 'Depth: '+`all_final_depth`+'\n'+'\n')
    vlast.append (yes)





            
    all_all_final_depth.append (all_final_depth)
    
    final_norm_error = norm_error /3**(0.5)
    xlim(np.min(beta),0)# this is just seting how wide the graph should be (so we are setting the domain)
    title('Normalized flux vs velocity')
    xlabel('Velocity (km/s)')
    ylabel('Normalized flux')
    plot((np.min(beta),np.max(beta)),(1,1))
    plot(beta, smooth(norm_flux,3), 'g-')
    plot (beta, final_norm_error,'k-')
    ylim(0,3)
    text(-40000, 2, `i`+',     z='+`z`+' snr='+`snr`, rotation=0, fontsize=9)
    
    pp.savefig()
    
    close(count)
    vmins_all.append(vmins)
    vmaxs_all.append(vmaxs)
    
    
BI_all= array(BI_all)

vmins = array(vmins)
vmaxs = array(vmaxs)
pp.close()
vmaxs_final=[]
vmins_final=[]

for loop in range (0, len (vmaxs_all)):
    
    
    vmaxs_final.append (`vmaxs_all[loop]`+ ',' )

for loop2 in range (0, len(vmins_all)):
    vmins_final.append (`vmins_all[loop2]`+ ',' )
                    

vmaxs_final = array(vmaxs_final)
vmins_final = array(vmins_final)    
savetxt('BI_EW_depth_dr10_dr10.txt',vlast,fmt='%s')




                    

    
        

       




    
