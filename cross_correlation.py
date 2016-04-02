import numpy as np
import pyfits as pf
import matplotlib.pyplot as plt
from numpy import mean, sqrt, square
import scipy as sp
from scipy import ndimage
import time



plt.ion()

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

spectra_full=[spectrum1_full,spectrum2_full,spectrum3_full,spectrum4_full,spectrum5_full,
spectrum6_full,spectrum7_full,spectrum8_full,spectrum9_full]

spectra=[]
for x in spectra_full:
	y=x[0]
	spectra.append(y[500:850])
	
window=15
spectra_medians=[]
normalized_spectra=[]

def normalize(x):
	medianrange=[]
	p=np.where(x==np.min(x))
	for i in range(len(x)):
		if i < (p[0]-window) or i > (p[0]+window):
			medianrange.append(x[i])
	median=np.median(medianrange)
	spectra_medians.append(median)
	normspec=x*(1-((median-spectra_medians[0])/spectra_medians[0]))/spectra_medians[0]
	normalized_spectra.append(normspec)
	
for x in spectra:
	normalize(x)




x=np.linspace(0,len(normalized_spectra[0]), len(normalized_spectra[0]))
xvals=np.linspace(0,len(normalized_spectra[0]),(len(normalized_spectra[0]))*10)
s1interp=np.interp(xvals,x,normalized_spectra[0])


def shift(s):
	sinterp=np.interp(xvals,x,s)
	chi=[]
	for i in range(len(s1interp)):
		n=np.roll(sinterp,i)
		diff= s1interp-n
		chi.append(np.sum(diff**2))
	
	p=np.where(chi==np.min(chi))[0]
 	a=np.roll(sinterp,p)
	plt.subplot(311)
	plt.plot(a)
	plt.subplot(312)
	diff=a-s1interp
	rms=np.around((np.mean(diff[1000:2000]**2))**(1/2.),decimals=7)
	plt.plot(diff)
	plt.ylim(-0.05,0.05)
	plt.subplot(313)
	plt.plot(rms,'o')
