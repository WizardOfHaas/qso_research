from pylab import*
from numpy import*
from sympy import sympify
from scipy import*
from astropy import*
from matplotlib.pyplot import*
from scipy.optimize import curve_fit
from matplotlib.backends.backend_pdf import PdfPages

normalized_dr9 =  loadtxt('data/dr9_flux/raw/spec-4743-55645-0118.dr9') #Load in normalized spectrum
thousand=0

wavelength_CIV_emit1=1550.7700
wavelength_CIV_emit2=1548.1950

wavelength = normalized_dr9[:,0] ###################  wavelength
norm_flux = normalized_dr9[:,1] ##################### normalized flux
norm_error = normalized_dr9[:,2] ###################normalized error
BI_mid=[]
#################################################Defining smoothing function
    
def smooth(norm_flux, box_pts):
       
    box = np.ones(box_pts)/box_pts
    y_smooth = convolve(norm_flux, box, mode='same')
    return y_smooth
        

    
    
z=1.9 #Rename redshift to something useful
z=round (z,5) #Round the redshift
sm_flux=smooth(norm_flux,3) #Smooth the spectrum
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
beta=[]
avr_CIV_doublet = 1549. #Average wavelength of CIV doublet
z_absC = (wavelength/avr_CIV_doublet)-1.
RC=(1.+z)/(1.+z_absC)
betaC=((RC**2.)-1.)/((RC**2.)+1.)
betaa = -betaC*(300000.)
for ll in betaa:
    betas=round (ll,4)
    beta.append (betas)
beta=array(beta)

maxvel= 0.#original
minvel=-60000#original

fst=np.max(where (beta <= maxvel))#index value of the starting point (on the very left) -- index value of minvel
try:
    lst=np.min(where(beta >=minvel))#index value of the ending point (on the very right) -- index value of maxvel
except:
    lst = where (beta == np.min(beta))

jjj=arange (lst,fst)
jjj=array(jjj)
jjj=jjj[::-1]

vs = []
for j in jjj:
	v = beta[j + 1] - beta[j]
	vs.append(v)

plot((beta[jjj[10] + 1],beta[jjj[10]]))
show()