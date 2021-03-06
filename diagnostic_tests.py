import numpy as np
import pyfits as pf
from astropy import table
import matplotlib.pyplot as plt
plt.ion()

#open spectra and store the data and header in arrays
iwavfile_full=pf.open('keck_iwav.fits')[0].data[0]
iwavfile=pf.open('keck_iwav.fits')[0].data[0,500:850]
spectra_filenames=[]
with open('spectra2.txt') as text:
	for line in text:
		spectra_filenames.append(line)

spectra_full=[]
spectra_headers=[]
spectra_files=[]
for name in spectra_filenames:
	spectra_full.append(pf.open(name)[0].data)
	spectra_headers.append(pf.open(name)[0].header)
	spectra_files.append(pf.open(name))
	
spectra=[]
for x in spectra_full:
 	a,b,c=np.polyfit(iwavfile_full,x[0],2)
 	quadratic_fit=[]
 	for i in iwavfile_full:
 		quadratic_fit.append(a*i**2+b*i+c)
	y=x[0]/quadratic_fit
	spectra.append(y[500:850])
	
#define function to normalize spectra
def normalize(x,norm,meanarray,plot=False,window=30):
	meanrange=[]
	p=np.where(x==np.min(x))
	for i in range(len(x)):
		if i < (p[0]-window) or i > (p[0]+window):
			meanrange.append(x[i])
	mean=np.mean(meanrange)
	meanarray.append(mean)
	normspec=x/mean
	norm.append(normspec)
	if plot==True:
		plt.plot(normspec)

spectra_medians=[]
normalized_spectra=[]
for x in spectra:
	normalize(x,normalized_spectra,spectra_medians)
	
#define interpolated, cross-correlation shift
x=np.linspace(0,len(normalized_spectra[0]), len(normalized_spectra[0]))
xvals=np.linspace(0,len(normalized_spectra[0]),(len(normalized_spectra[0]))*10)
s1interp=np.interp(xvals,x,normalized_spectra[0])
iwavfileinterp=np.interp(xvals,x,iwavfile)

#define a shift function that plots shifted spectra and its difference spectrum
shifted_spectra=[]
shifts=[]
def shift(s,store=False, plot=True,shifted_specs=shifted_spectra,store_shifts=False,
shifts_array=shifts,plotreg=False):
	sinterp=np.interp(xvals,x,s)
	chi=[]
	for i in range(len(s1interp)):
		n=np.roll(sinterp,i)
		diff= s1interp-n
		chi.append(np.sum(diff**2))
	p=np.where(chi==np.min(chi))[0]
	if store_shifts==True:
		if p<(len(s1interp)/2.):
			shifts_array.append(p[0])
		if p>(len(s1interp)/2.):
			shifts_array.append(-1*(len(s1interp)/2.))
 	a=np.roll(sinterp,p)
 	if store == True:
 		shifted_specs.append(a)
 	if store == False:
 		None
	if plot == True:
		s1=np.where(s1interp==np.min(s1interp))[0]-500
		s2=np.where(s1interp==np.min(s1interp))[0]+1000
		p1=iwavfileinterp[s1]
		p2=iwavfileinterp[s2]
		plt.subplot(211)
		plt.plot(iwavfileinterp,a)
		plt.axvline(x=p1,linestyle='--',color='black')	
		plt.axvline(x=p2,linestyle='--',color='black')
		plt.title('Spectra')
		plt.ylabel('Normalized Flux')
		plt.axvline(x=iwavfileinterp[200],linestyle='--',color='black')
		plt.axvline(x=iwavfileinterp[3300],linestyle='--',color='black')
		
		plt.subplot(212)
		diff=a-s1interp
		plt.plot(iwavfileinterp,diff)
		#plt.axvline(x=s1,linestyle='--',color='black')	
		#plt.axvline(x=s2,linestyle='--',color='black')
		plt.ylim(-0.05,0.05)
		plt.title('Difference Spectra')
		plt.axvline(x=p1,linestyle='--',color='black')	
		plt.axvline(x=p2,linestyle='--',color='black')
		plt.axvline(x=iwavfileinterp[200],linestyle='--',color='black')
		plt.axvline(x=iwavfileinterp[3300],linestyle='--',color='black')
	else:
		None
	if plotreg==True:
		plt.plot(iwavfile,s)
		
#define a function to calculate RMS values
def calc_RMS(RMS1,RMS2,RMS3):
	for j in normalized_spectra:
		sinterp=np.interp(xvals,x,j)
		chi=[]
		for i in range(len(s1interp)):
			n=np.roll(sinterp,i)
			diff= s1interp-n
			chi.append(np.sum(diff**2))
		p=np.where(chi==np.min(chi))[0]
 		a=np.roll(sinterp,p)
 		diff=a-s1interp
 		point=np.where(s1interp==np.min(s1interp))[0]
 		rms1=(np.mean(diff[point-800:point-500]**2))**(1/2.)
		rms2=(np.mean(diff[point-500:point+1000]**2))**(1/2.)
		rms3=(np.mean(diff[point+1000:point+1300]**2))**(1/2.)
		RMS1.append(rms1)
		RMS2.append(rms2)
		RMS3.append(rms3)


#the next part of the code will plot 6 diagnostic plots

#store JD and phase in arrays
P=2.218575200
Tdur=0.075152
MP=2453955.5255511
JD=[]
phase=[]
for spec in spectra_headers:
	jd=float(spec[59])+2400000.5
	JD.append(jd)
	phase.append(((jd-(MP-P/2.))%P)/P)

#store signal:noise in array
shiftsarray1=[]
shifted_tablespecs=[]
for n in spectra:
	shift(n,store=True, plot=False,shifted_specs=shifted_tablespecs,store_shifts=True,
	shifts_array=shiftsarray1)
	shiftsarray=np.array(shiftsarray1)/10.0
SN=[]
for spec in shifted_tablespecs:
	min=np.where(spec==np.min(spec))[0]
	sn1=[]
	sn1.append(spec[min-1000:min-500])
	sn1.append(spec[min+1000:min+1500])
	sn=np.array(sn1).ravel()
	SN.append((np.median(2*sn)**(1/2.)))
	
#store RMS values in arrays
rms1=[]
rms2=[]
rms3=[]
calc_RMS(rms1,rms2,rms3)
rms1med=str(np.around(np.median(rms1),decimals=4))
rms3med=str(np.around(np.median(rms3),decimals=4))

#make plots :)
def rms_plot():
	plt.title('RMS for Spectrum Sections 1 and 3')	
	plt.plot(phase,rms1,'bo',label='section 1'+' median='+rms1med)
	plt.plot(phase,rms3,'ro',label='section 3'+' median='+rms3med)
	plt.xlabel('phase')
	plt.ylabel('Section RMS')
	plt.legend()
	plt.legend(bbox_to_anchor=(1.0,1.09))
	
def diagnostic_plot():
	plt.suptitle('Diagnostic Plots')
	plt.title('Spectra Diagnostic Plots')
	plt.subplot(231)
	plt.plot(JD,SN,'*')
	plt.xlabel('JD')
	plt.ylabel('S:N')
	plt.subplot(232)
	plt.plot(JD,shiftsarray,'*')
	plt.xlabel('JD')
	plt.ylabel('Pixel Shift')
	plt.subplot(233)
	plt.plot(JD,rms1,'bo',label='section 1')
	plt.plot(JD,rms3,'ro',label='section 3')
	plt.ylabel('Section RMS')
	plt.xlabel('JD')
	plt.legend()
	plt.legend(bbox_to_anchor=(1.05,1.05))
	plt.subplot(234)
	plt.plot(phase,SN,'*')
	plt.ylabel('S:N')
	plt.xlabel('phase')
	plt.subplot(235)
	plt.plot(phase,shiftsarray,'*')
	plt.ylabel('Pixel Shift')
	plt.xlabel('phase')
	plt.subplot(236)	
	plt.plot(phase,rms1,'bo',label='section 1')
	plt.plot(phase,rms3,'ro',label='section 3')
	plt.xlabel('phase')
	plt.ylabel('Section RMS')
	plt.legend()
	plt.legend(bbox_to_anchor=(1.15,1.15))
	
def histogram_plot():
	plt.subplot(121)
	ingress=(((MP-Tdur)-(MP-P/2.))%P)/P
	egress=(((MP+Tdur)-(MP-P/2.))%P)/P
	plt.axvspan(ingress,egress,color='red',alpha=0.4,hatch='xxx')
	plt.hist(phase,bins=np.linspace(0,1,50),alpha=0.5)
	plt.title('Orbital Phase Distribution')
	plt.xlabel('phase')
	plt.ylabel('Spectra')
	plt.subplot(122)
	plt.hist(SN,alpha=0.5)
	plt.title('S:N Distribution')
	plt.xlabel('S:N')
	plt.ylabel('Spectra')
	
def shift_plot():
	plt.plot(shiftsarray,rms1,'bo',label='section 1')
	plt.plot(shiftsarray,rms2,'go',label='section 2')
	plt.plot(shiftsarray,rms3,'ro',label='section 3')
	plt.title('Spectra Shifts Versus RMS')
	plt.xlabel('Pixel Shift')
	plt.ylabel('RMS')
	plt.legend(bbox_to_anchor=(1.05,1.05))
	
def plot_spectra():
	shifted_specs_plot=[]
	shifts_array_plot=[]
	for spec in normalized_spectra:
		shift(spec,store=False, plot=True,shifted_specs=shifted_specs_plot,store_shifts=False,
		shifts_array=shifts_array_plot)









