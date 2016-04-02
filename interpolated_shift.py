import numpy as np
import scipy as sp
import pyfits as pf
import matplotlib.pyplot as plt
from numpy import mean, sqrt, square
from scipy import interpolate

iwavfile_full=pf.open('keck_iwav.fits')[0].data
iwavfile=iwavfile_full[0]

#open all the files and store data
#out of transit
spectrum1_full=pf.open('ij22.213.fits')[0].data
spectrum2_full=pf.open('ij22.214.fits')[0].data
spectrum3_full=pf.open('ij22.215.fits')[0].data
spectrum4_full=pf.open('ij22.216.fits')[0].data
spectrum5_full=pf.open('ij22.217.fits')[0].data
#in transit
spectrum6_full=pf.open('ij02.929.fits')[0].data
spectrum7_full=pf.open('ij02.930.fits')[0].data
spectrum8_full=pf.open('ij08.1997.fits')[0].data
spectrum9_full=pf.open('ij13.1014.fits')[0].data
spectrum10_full=pf.open('ij01.1033.fits')[0].data

#make a list of spectra to work with
spectra_full=[spectrum1_full,spectrum2_full,spectrum3_full,spectrum4_full,spectrum5_full,
spectrum6_full,spectrum7_full,spectrum8_full,spectrum9_full]
spectra=[]
for x in spectra_full:
	spectra.append(x[0])

#interpolate spectra
wavelength_range=iwavfile[-1]-iwavfile[0]
pixel_range=len(iwavfile)
interval=wavelength_range/(pixel_range-1)
splined_wavelength=iwavfile[0]+np.arange(pixel_range)*interval
interpolated_spectra=[]
for x in spectra:
	tck=interpolate.splrep(iwavfile,x)
	ynew=interpolate.splev(splined_wavelength,tck)
	interpolated_spectra.append(ynew)


#normalize spectra
wavelength=6563
window=5
#note: the median is calculated for wavelengths 1 angstrom from the absorption feature 
#wavelength. For wider or thinner features, adjust this number
spectra_medians=[]
normalized_spectra=[]
for x in interpolated_spectra:
	median_rangex=[]
	for n in range(len(x)-1):
			if iwavfile[n]<(wavelength+1)+window:
				if iwavfile[n]>(wavelength-1)-window:
					median_rangex.append(x[n])
	medianx=np.median(median_rangex)
	spectra_medians.append(medianx)
#	thresholdx=0.8*medianx
	def normalize(spectrum):
		return spectrum*(1-((medianx-spectra_medians[0])/spectra_medians[0]))/spectra_medians[0]
	array=np.array(normalize(x))
	normalized_spectra.append(array)
	

#shift spectra horizontally
#first define the parabola fitting function
vertices=[]
spectrum_min_indices=[]
def fit_parabola(spectrum,array):
	#calculate the minimum to know where to fit the parabola
	spectrum_min=[]
	vertex=0
	spectrum_min_final=0
	for n in range(len(spectrum)-1):
		if spectrum[n]<spectrum[n+1]:
			if spectrum[n]<spectrum[n-1]:
 				if splined_wavelength[n]>6562.4:
 					if splined_wavelength[n]<6563.8:
						spectrum_min.append(spectrum[n])
						if spectrum[n]==np.min(spectrum_min):
							spectrum_min_final=splined_wavelength[n]
#	spectrum_min_final=np.median(spectrum_min)
# 	#fitting parabola around minimum
 	wavfitarray=[]
 	spectrumfitarray=[]
 	for n in range(len(spectrum)-1):
 			if splined_wavelength[n]>spectrum_min_final-4:
 				if splined_wavelength[n]<spectrum_min_final+4:
 					wavfitarray.append(splined_wavelength[n])
 					spectrumfitarray.append(spectrum[n])
 	a,b,c=np.polyfit(wavfitarray,spectrumfitarray,2)
	quadratic_fit=[]
	for i in wavfitarray:
		quadratic_fit.append(a*i**2+b*i+c)
	x=-b/(2*a)
	point=quadratic_fit[np.abs(quadratic_fit - x).argmin()]
	index = quadratic_fit.index(point)
	spectrum_min_indices.append(index)
	array.append(-b/(2*a))
	#plt.plot(wavfitarray,quadratic_fit)
	
#define wavelength boundaries for clipped sections
wavelength_cut=[]
for n in range(len(splined_wavelength)):
	if splined_wavelength[n]<6565.2:
		if splined_wavelength[n]>6561.2:
			wavelength_cut.append(splined_wavelength[n])

section1_min=0
section1_max=0
for n in range(len(splined_wavelength)):
		if splined_wavelength[n]==wavelength_cut[0]:
			section1_min=n
		if splined_wavelength[n]==wavelength_cut[-1]:
			section1_max=n

	
#shift the spectra horizontally
clipped_spectra=[]
for x in normalized_spectra:
	fit_parabola(x,vertices)
	spectrum_shiftx=(vertices[-1]-vertices[0])*33.65
	clipped_spectra.append(x[section1_min+spectrum_shiftx:section1_max+1+spectrum_shiftx])

plt.subplot(211)
for x in clipped_spectra:
	plt.plot(wavelength_cut,x)
	plt.xlim(6561.5,6565)
	
diff=clipped_spectra[0]-clipped_spectra[8]
s1=diff[0:len(wavelength_cut)/3]
s2=diff[len(wavelength_cut)/3:2*len(wavelength_cut)/3]
s3=diff[2*len(wavelength_cut)/3:wavelength_cut[-1]]
rms1=np.around((np.mean(s1**2))**(1/2.),decimals=7)
rms2=np.around((np.mean(s2**2))**(1/2.),decimals=7)
rms3=np.around((np.mean(s3**2))**(1/2.),decimals=7)


plt.subplot(212)
for x in clipped_spectra:
	plt.plot(wavelength_cut,x-clipped_spectra[0])
	plt.xlim(6561.5,6565)
plt.figtext(0.17,0.40,str(rms1),fontsize=12)
plt.figtext(0.47,0.40,str(rms2),fontsize=12)
plt.figtext(0.77,0.40,str(rms3),fontsize=12)

plt.show()
