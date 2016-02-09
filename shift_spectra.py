import numpy as np
import pyfits as pf
import matplotlib.pyplot as plt
from numpy import mean, sqrt, square
import matplotlib.patches as mpatches

iwavfile=pf.open('keck_iwav.fits')[0].data

#open all the files and store data
#out of transit
spectrum1=pf.open('ij22.213.fits')[0].data
spectrum2=pf.open('ij22.214.fits')[0].data
spectrum3=pf.open('ij22.215.fits')[0].data
spectrum4=pf.open('ij22.216.fits')[0].data
spectrum5=pf.open('ij22.217.fits')[0].data

#in transit
spectrum6=pf.open('ij02.929.fits')[0].data
spectrum7=pf.open('ij02.930.fits')[0].data
spectrum8=pf.open('ij08.1997.fits')[0].data
spectrum9=pf.open('ij13.1014.fits')[0].data
spectrum10=pf.open('ij01.1033.fits')[0].data

#put the files into a list
spectra=[spectrum1,spectrum2,spectrum3,spectrum4,spectrum5,spectrum6,spectrum7,spectrum8,
spectrum9]

#CALCULATE AND PERFORM VERTICAL SHIFTS

#calculating threshold for minimum calculation
wavelength=6563
window=5
median_range1=[]
for n in range(len(spectrum1[:,0])):
	spectrum1_n= spectrum1[n,::]
	iwavfile_n=iwavfile[n,::]
	for i in range(len(spectrum1_n)):
		if iwavfile_n[i]<(wavelength+1)+window:
			if iwavfile_n[i]>(wavelength-1)-window:
				median_range1.append(spectrum1_n[i])
median1=np.median(median_range1)				
threshold=0.8*median1

#calculate median y value for each spectra
#calculating median of spectra
spectra_medians=[]
normalized_spectra=[]
for x in spectra:
	median_rangex=[]
	for n in range(len(x[:,0])-1):
		x_n=x[n,::]
		iwavfile_n=iwavfile[n,::]
		for i in range(len(x_n)-1):
			if iwavfile_n[i]<(wavelength+1)+window:
				if iwavfile_n[i]>(wavelength-1)-window:
					median_rangex.append(x_n[i])
	medianx=np.median(median_rangex)
	spectra_medians.append(medianx)
	thresholdx=0.8*medianx
	def normalize(spectrum):
		return spectrum*(1-((medianx-spectra_medians[0])/spectra_medians[0]))/spectra_medians[0]
	array=np.array(normalize(x))
	normalized_spectra.append(array)

#CALCULATE HORIZONTAL SHIFTS
vertices=[]
spectrum_min_indices=[]
def fit_parabola(spectrum,array):
	#calculate the minimum to know where to fit the parabola
	spectrum_min=[]
	for n in range(len(spectrum[:,0])-1):
		spectrum_n= spectrum[n,::]
		iwavfile_n=iwavfile[n,::]
		for i in range(len(spectrum_n)-1):
			if spectrum_n[i]<spectrum_n[i+1]:
				if spectrum_n[i]<spectrum_n[i-1]:
					if spectrum_n[i]<threshold:
						if iwavfile_n[i]>6562.4:
							if iwavfile_n[i]<6563.8:
								spectrum_min.append(iwavfile_n[i])
	spectrum_min_final=np.min(spectrum_min)
	
	
	#fitting parabola around minimum
	wavfitarray=[]
	spectrumfitarray=[]
	for n in range(len(spectrum[:,0])-1):
		spectrum_n=spectrum[n,::]
		iwavfile_n=iwavfile[n,::]
		for i in range(len(spectrum_n)-1):
			if iwavfile_n[i]>spectrum_min_final-0.7:
				if iwavfile_n[i]<spectrum_min_final+0.7:
					wavfitarray.append(iwavfile_n[i])
					spectrumfitarray.append(spectrum_n[i])
	a,b,c=np.polyfit(wavfitarray,spectrumfitarray,2)
	quadratic_fit=[]
	for i in wavfitarray:
		quadratic_fit.append(a*i**2+b*i+c)
	#x=np.min(quadratic_fit)
	x=-b/(2*a)
	point=quadratic_fit[np.abs(quadratic_fit - x).argmin()]
	index = quadratic_fit.index(point)
	spectrum_min_indices.append(index)
	array.append(-b/(2*a))
	#plt.plot(quadratic_fit)

iwavfile_cut=[]
for n in range(len(iwavfile[:,0])):
	iwavfile_n=iwavfile[n,::]
	for i in range(len(iwavfile_n)):
		if iwavfile_n[i]<6565.2:
			if iwavfile_n[i]>6561.2:
				iwavfile_cut.append(iwavfile_n[i])

#clip spectra around spectral feature
#first define the boundaries of the clipped section
section1_min=[]
section1_max=[]
for n in range(len(iwavfile[:,0])):
	iwavfile_n=iwavfile[n,::]
	for i in range(len(iwavfile_n)):
		if iwavfile_n[i]==iwavfile_cut[0]:
			section1_min.append(i)
		l1=len(iwavfile_cut)-1
		if iwavfile_n[i]==iwavfile_cut[l1]:
			section1_max.append(i)

clipped_spectra=[]
for x in normalized_spectra:
	fit_parabola(x,vertices)
	spectrum_shiftx=(vertices[-1]-vertices[0])*33.65
	clipped_spectra.append(x[0,section1_min[0]+spectrum_shiftx:section1_max[0]+1+spectrum_shiftx])

#Spectra Plot
plt.subplot(211)
for x in clipped_spectra:
	plt.plot(x)
plt.xlim(0,120)
plt.ylabel('Normalized Flux')
plt.title('Transit H-Alpha Spectra for HD189733')
	
#Difference Plot
diff=clipped_spectra[0]-clipped_spectra[1]
s1=diff[0:30]
s2=diff[30:90]
s3=diff[90:120]
rms1=np.around((np.mean(s1**2))**(1/2.),decimals=7)
rms2=np.around((np.mean(s2**2))**(1/2.),decimals=7)
rms3=np.around((np.mean(s3**2))**(1/2.),decimals=7)
	
plt.subplot(212)
for x in clipped_spectra:	
	plt.plot(x-clipped_spectra[0])
plt.xlim(0,120)
plt.ylim(-0.1,0.1)
plt.figtext(0.15,0.15,'section 1',fontsize=12)
plt.figtext(0.45,0.15,'section 2',fontsize=12)
plt.figtext(0.75,0.15,'section 3',fontsize=12)
plt.figtext(0.15,0.35,str(rms1),fontsize=12)
plt.figtext(0.45,0.35,str(rms2),fontsize=12)
plt.figtext(0.75,0.35,str(rms3),fontsize=12)
plt.ylabel('Normalized Flux')
plt.xlabel('Wavelength')
plt.title('H-Alpha Difference Spectra for HD189733')
	
plt.axvline(x=30,color='black',linestyle='--')
plt.axvline(x=90,color='black',linestyle='--')
plt.show()

#RMS Plot
RMS_array_s1=[]
RMS_array_s2=[]
RMS_array_s3=[]
#plotting the rms values for each section in each array
for x in clipped_spectra:
	diff=clipped_spectra[0]-x
	s1=diff[0:30]
	s2=diff[30:90]
	s3=diff[90:120]
	rms1=np.around((np.mean(s1**2))**(1/2.),decimals=7)
	rms2=np.around((np.mean(s2**2))**(1/2.),decimals=7)
	rms3=np.around((np.mean(s3**2))**(1/2.),decimals=7)
	RMS_array_s1.append(rms1)
	RMS_array_s2.append(rms2)
	RMS_array_s3.append(rms3)
	plt.plot(RMS_array_s1,'co')
	plt.plot(RMS_array_s2,'mo')
	plt.plot(RMS_array_s3,'yo')
	plt.xlabel('Spectrum')
	plt.ylabel('RMS')
	plt.title('RMS Values for Spectra Sections')
	plt.xlim(0,9)

plt.plot(RMS_array_s1,'co',label='section 1')
plt.plot(RMS_array_s2,'mo',label='section 2')
plt.plot(RMS_array_s3,'yo',label='section 3')
plt.legend(loc='upper left', numpoints=1)
#bbox_to_anchor=(1.1, 1) add this command to move legend

plt.show()
