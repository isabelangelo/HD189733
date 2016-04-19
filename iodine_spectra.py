import numpy as np
import pyfits as pf
from astropy import table
import matplotlib.pyplot as plt
plt.ion()

#open spectra and store the data and header in arrays
iwavfile=pf.open('keck_iwav.fits')[0].data[0]
flat_o=pf.open('ij55.254.fits')[0].data[0]
flat_t=pf.open('ij97.2716.fits')[0].data[0]
star_o=pf.open('ij90.734.fits')[0].data[0]
star_t=pf.open('ij90.744.fits')[0].data[0]

iwavfilec=pf.open('keck_iwav.fits')[0].data[0,500:800]
flat_oc=pf.open('ij55.254.fits')[0].data[0,500:800]
flat_tc=pf.open('ij97.2716.fits')[0].data[0,500:800]
star_oc=pf.open('ij90.734.fits')[0].data[0,500:800]
star_tc=pf.open('ij90.744.fits')[0].data[0,500:800]


#fit a parabola to the entire thing and divide by it
def normalize(x):
	a,b,c=np.polyfit(iwavfile,x,2)
	quadratic_fit=[]
	for i in iwavfile:
		quadratic_fit.append(a*i**2+b*i+c)
	flatspec=x/quadratic_fit
	return flatspec/np.median(flatspec)

def normalize_c(x):
	a,b,c=np.polyfit(iwavfilec,x,2)
	quadratic_fit=[]
	for i in iwavfilec:
		quadratic_fit.append(a*i**2+b*i+c)
	flatspec=x/quadratic_fit
	return flatspec/np.mean(flatspec)


#overplots normalized full spectra
def plot_raw_spectra():
	plt.suptitle('Full Spectra')
	plt.subplot(221)
	plt.title('Quartz Lamp With Iodine')
	plt.xlabel('Wavelength (Angstroms)')
	plt.ylabel('Flux')
	plt.plot(iwavfile,normalize(flat_o))

	plt.subplot(222)
	plt.title('Quartz Lamp Without Iodine')
	plt.xlabel('Wavelength (Angstroms)')
	plt.ylabel('Flux')
	plt.plot(iwavfile,normalize(flat_t))

	plt.subplot(223)
	plt.plot(iwavfile,normalize(star_o))
	plt.title('B Star With Iodine')
	plt.xlabel('Wavelength (Angstroms)')
	plt.ylabel('Flux')

	plt.subplot(224)
	plt.plot(iwavfile,normalize(star_t))
	plt.title('B Star Without Iodine')
	plt.xlabel('Wavelength (Angstroms)')
	plt.ylabel('Flux')

#overplots clipped spectra
def plot_clipped_spectra():
	plt.suptitle('Clipped Spectra')
	plt.subplot(221)
	plt.title('Quartz Lamp With Iodine')
	plt.xlabel('Wavelength (Angstroms)')
	plt.ylabel('Flux')
	plt.plot(iwavfilec,flat_oc)

	plt.subplot(222)
	plt.title('Quartz Lamp Without Iodine')
	plt.xlabel('Wavelength (Angstroms)')
	plt.ylabel('Flux')
	plt.plot(iwavfilec,flat_tc)

	plt.subplot(223)
	plt.title('Star With Iodine')
	plt.xlabel('Wavelength (Angstroms)')
	plt.ylabel('Flux')
	plt.plot(iwavfilec,star_oc)

	plt.subplot(224)
	plt.title('Star Without Iodine')
	plt.xlabel('Wavelength (Angstroms)')
	plt.ylabel('Flux')
	plt.plot(iwavfilec,star_tc)

#overplots normalized clipped spectra
def overplot_spectra():
	plt.subplot(121)
# 	x1=6560
# 	y1=0.96
# 	x2=6566
# 	y2=1.04
# 	x=np.arange(0,len(iwavfilec),1)
# 	y=((y2-y1)/(x2-x1))*(x-x1)+y1
#	ynorm=y/np.median(y)
	plt.title('Quartz Lamp Spectra')
	plt.xlabel('Wavelength (Angstroms)')
	plt.ylabel('Flux')
	plt.plot(iwavfilec,normalize_c(flat_oc),label='iodine')
	plt.plot(iwavfilec,normalize_c(flat_tc),label='no iodine')
	plt.xlim(6559,6567)
	plt.ylim(0.98,1.017)
	plt.legend()
	
	plt.subplot(122)
# 	x1=6559
# 	y1=1.02
# 	x2=6567
# 	y2=1.09
# 	x=np.arange(0,len(iwavfilec),1)
# 	y=((y2-y1)/(x2-x1))*(x-x1)+y1
	plt.title('B Star Spectra')
	plt.xlabel('Wavelength (Angstroms)')
	plt.ylabel('Flux')
	plt.plot(iwavfilec,normalize_c(star_oc),label='iodine')
	plt.plot(iwavfilec,normalize_c(star_tc),label='no iodine')
	plt.xlim(6559,6567)
	plt.ylim(0.2,1.4)
	plt.legend()
	
#want to use test method for diagnostics because we are working with spectra
def test():
	plt.suptitle('Full Spectra')
	plt.subplot(121)
	plt.title('Quartz Lamp With Iodine')
	plt.xlabel('Wavelength (Angstroms)')
	plt.ylabel('Flux')
	plt.plot(iwavfile,normalize(flat_o))
	plt.plot(iwavfile,normalize(flat_t))
	plt.xlim(6559,6567)

	plt.subplot(122)
	plt.plot(iwavfile,normalize(star_o))
	plt.title('B Star With Iodine')
	plt.xlabel('Wavelength (Angstroms)')
	plt.ylabel('Flux')
	plt.plot(iwavfile,normalize(star_t))
	plt.xlim(6559,6567)

def final_plot():
#everything except star spectrum from original overplot, star spectrum from test
	plt.suptitle('Iodine Spectra')
	plt.subplot(121)
	plt.title('Quartz Lamp')
	plt.plot(iwavfilec,normalize_c(flat_oc),label='iodine')
	plt.plot(iwavfilec,normalize_c(flat_tc),label='no iodine')
	plt.xlim(6559,6567)
	plt.ylim(0.98,1.017)
	plt.xlabel('Wavelength (Angstroms)')
	plt.ylabel('Normalized Flux')
	plt.legend(bbox_to_anchor=(1.05,1.05))
	
	plt.subplot(122)
	plt.title('B Star')
	plt.plot(iwavfile[500:850],normalize(star_o)[500:850],label='iodine')
	plt.plot(iwavfilec,normalize_c(star_tc),label='no iodine')
	plt.xlim(6559,6567)
	plt.xlabel('Wavelength (Angstroms)')
	plt.ylabel('Normalized Flux')
	plt.legend(bbox_to_anchor=(1.05,1.05))
	
	
	
	
	
	