import numpy as np
import pyfits as pf
from astropy import table
import matplotlib.pyplot as plt
from astropy.table import Table, Column
plt.ion()

#open spectra and store the data and header in arrays
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
	y=x[0]
	spectra.append(y[500:850])
	
#normalize spectra and store them into an array "norm"
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
RMS=[]
def shift(s,store=False, plot=True,shifted_specs=shifted_spectra,store_shifts=False,shifts_array=shifts,RMS=RMS):
	sinterp=np.interp(xvals,x,s)
	chi=[]
	for i in range(len(s1interp)):
		n=np.roll(sinterp,i)
		diff= s1interp-n
		chi.append(np.sum(diff**2))
	p=np.where(chi==np.min(chi))[0]
	if store_shifts==True:
		shifts_array.append(p[0])
 	a=np.roll(sinterp,p)
 	if store == True:
 		shifted_specs.append(a)
 	if store == False:
 		None
	if plot == True:
		plt.subplot(211)
		plt.plot(a)
		plt.subplot(212)
		diff=a-s1interp
		rms=np.around((np.mean(diff[1000:2000]**2))**(1/2.),decimals=7)
		RMS.append(rms)
		plt.plot(diff)
		plt.ylim(-0.05,0.05)
	else:
		None
	
#PLOT RMS VALUES
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
 		rms1=(np.mean(diff[point-1000:point-500]**2))**(1/2.)
		rms2=(np.mean(diff[point-500:point+1000]**2))**(1/2.)
		rms3=(np.mean(diff[point+1000:point+1500]**2))**(1/2.)
		RMS1.append(rms1)
		RMS2.append(rms2)
		RMS3.append(rms3)
		

def plot_rms():
	RMS11=[]
	RMS21=[]
	RMS31=[]
	calc_RMS(RMS11, RMS21, RMS31)
	specx=np.arange(1,len(spectra_full)+1,1)
	plt.plot(specx,RMS11,'bo',label='section 1')
	plt.plot(specx,RMS21,'go',label='section 2')
	plt.plot(specx,RMS31,'ro',label='section 3')
	plt.xlim(0,len(normalized_spectra)+1)
	plt.ylim(0,0.05)
	plt.title('Spectra RMS Values')
	plt.xlabel('Spectrum Index')
	plt.ylabel('RMS')
	plt.legend()
	
#PLOT RMS AS FUNCTION OF PHASE	
def plot_rmsphase():
	RMS12=[]
	RMS22=[]
	RMS32=[]	
	calc_RMS(RMS12, RMS22, RMS32)
	JD=[]
	for x in spectra_headers:
		JD.append(float(x[59])+2400000.5)
	plt.plot(JD,RMS12,'bo',label='section 1')
	plt.plot(JD,RMS22,'go',label='section 2')
	plt.plot(JD,RMS32,'ro',label='section 3')
	plt.title('RMS versus Julian Date')
	plt.xlabel('Julian Date')
	plt.ylabel('RMS')
	#convert to phase: find distance to midpoint of transit, convert to radians
	MP=2453955.5255511
	d=[]
	for x in JD:
		d.append(x-MP)
	#plt.plot(d,RMS,'ro')
	
#plot difference images for averaged in transit and out of transit spectra

#sort spectra
in_transit_spectra=[]
out_of_transit_spectra=[]
P=2.218575200
Tdur=0.075152
MP=2453955.5255511
def sort(spectrum):
	JD=float(spectrum[0].header[59])+2400000.5
	for n in range(int(JD/100)):
		if JD+P*n>=(MP-Tdur/2):
			if JD+P*n<=(MP+Tdur/2):
				in_transit_spectra.append(spectrum)
				break
		elif JD-P*n>=(Tdur/2-MP):
			if JD-P*n<=(Tdur/2+MP):
				in_transit_spectra.append(spectrum)
				break
	else:
		out_of_transit_spectra.append(spectrum)
		
for spectrum in spectra_files:
	sort(spectrum)

def plot_averaged_spectra():
	specdata=[]
	for x in in_transit_spectra:
		specdata.append(x[0].data[0,500:850])
	normspec=[]
	medspec=[]
	for spec in specdata:
		normalize(spec,normspec,medspec)
	shifted_intransitspec=[]
	for n in normspec:
		shift(n,store=True,plot=False,shifted_specs=shifted_intransitspec)
		
	in_avg=[]
	for i in range(len(shifted_intransitspec[0])):
		vals=[]
		for spectrum in shifted_intransitspec:
			vals.append(spectrum[i])
		in_avg.append(np.mean(vals))

	specdata2=[]
	for x in out_of_transit_spectra:
		specdata2.append(x[0].data[0,500:850])
	normspec2=[]
	medspec2=[]
	for spec in specdata2:
		normalize(spec,normspec2,medspec2)
	shifted_out_of_transitspec=[]
	for n in normspec2:
		shift(n,store=True,plot=False,shifted_specs=shifted_out_of_transitspec)
	out_avg=[]
	for i in range(len(shifted_out_of_transitspec[0])):
		vals=[]
		for spectrum in shifted_out_of_transitspec:
			vals.append(spectrum[i])
		out_avg.append(np.mean(vals))
	
	plt.subplot(211)
	plt.plot(in_avg,'g-',label='In Transit')
	plt.plot(out_avg,'b-',label='Out of Tranist')
	plt.ylabel('Normalized Flux')
	plt.title('Averaged Spectra')
	
	plt.subplot(212)
	plt.plot(np.array(in_avg)-np.array(out_avg),'g-')
	plt.plot(np.array(out_avg)-np.array(out_avg),'b-')
	plt.ylim(-0.1,0.1)
	plt.xlabel('Wavelength')
	plt.ylabel('Normalized Flux')
	plt.title('Difference Spectra')
	plt.legend()
	
def plot_table():
	#create an array SN that calculates s:n for each spectrum
	shiftsarray1=[]
	shifted_tablespecs=[]
	for n in spectra:
		shift(n,store=True, plot=False,shifted_specs=shifted_tablespecs,store_shifts=True,shifts_array=shiftsarray1)
	shiftsarray=np.array(shiftsarray1)/10.0
	SN=[]
	for spec in shifted_tablespecs:
		min=np.where(spec==np.min(spec))[0]
		sn1=[]
		sn1.append(spec[min-1000:min-500])
		sn1.append(spec[min+1000:min+1500])
		sn=np.array(sn1).ravel()
		SN.append((np.median(2*sn)**(1/2.)))
	#create an array of orbital phases
	phase=[]
	for spec in spectra_headers:
		JD=float(spec[59])+2.44*10**6
		phase.append(((JD-MP)%P)/P)
	#create arrays for RMS for sections 1,2, and 3
	RMS1tab=[]
	RMS2tab=[]
	RMS3tab=[]
	calc_RMS(RMS1tab,RMS2tab,RMS3tab)
	#plot arrays in a table
	t = Table([spectra_filenames,np.arange(0,len(spectra),1),SN,shiftsarray,phase,RMS1tab,
	RMS2tab,RMS3tab],names=('File Name', 'Spectrum Index', 'S:N','Interpolated Pixel Shift','Phase',
	'RMS1','RMS2','RMS3'))
	t.show_in_browser(jsviewer=True)
	
	
		
		
	
	
	
	
	
	
	