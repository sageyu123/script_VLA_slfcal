import os
import numpy as np
import vla_prep
import shutil
import pdb
import time
from scipy.signal import medfilt2d



start_timestamp = time.time()
# specfile='SUN01_20141101TT164613.125-164629.525.50ms.slfcal.ms.spec.npz'
specfile='SUN01_20141101TT164613.125-164629.525.spw1.50ms.slfcal.ms.spec.npz'
specdata=np.load(specfile)
spec=specdata['spec']
npol=specdata['npol']
nbl=specdata['nbl']
ntim=specdata['ntim']
nfreq=specdata['nfreq']
tim=specdata['tim']
freq=specdata['freq']
bl=specdata['bl'].item()
b=0
pol='LL'
sz_spec=spec.shape
spec_plt=np.zeros((4,sz_spec[2],sz_spec[3]))
spec_plt[0,:,:]=spec[1,b,:,:]
# spec_plt[1,:,:]=spec[1,b,:,:]
# spec_plt[2,:,:]=spec[1,b,:,:]
spec_plt[3,:,:]=spec_med = medfilt2d(spec[1,b,:,:],kernel_size=[3,3])

threshold=250
# threshold=300
idx_thresh1 = spec_plt[3,:,:] > threshold
# idx_thresh2 = spec_plt[0,:,:] <= threshold
spec_plt[1,:,:][idx_thresh1] = 400
spec_plt[2,:,:][idx_thresh1] = spec_plt[0,:,:][idx_thresh1]
idx_selec = np.where(spec_plt[3,:,:] > threshold)
idx_tim=np.unique(idx_selec[1])


mspath='/srg/sjyu/20141101/'
os.chdir(mspath)
msfile = 'SUN01_20141101.T163940-164700.50ms.cal.ms'
ephemfile='horizons_sun_20141101.radecp'
msinfofile='SUN01_20141101.T163940-164700.50ms.cal.msinfo.npz'
slfcalms_timeran='16:46:13.125~16:46:29.525'
slfcalms_s=['1']
slfcalms_chan=['0~63']
structure_id_slfcalms = 'rU01'
spwchan = ','.join('%s:%s' % t for t in zip(slfcalms_s,slfcalms_chan))

''' ----- step 1 ----- '''


dir_ms_split='./ms_split/'
if not os.path.exists(dir_ms_split):
    os.mkdir(dir_ms_split)

strtmp=[t.replace(':','') for t in slfcalms_timeran.split('~')]
timestr='T'+strtmp[0]+'-'+strtmp[1]
slfcalms = dir_ms_split+'SUN01_20141101T'+timestr+'.50ms.slfcal.ms'
slfcal_dbkg_ms = dir_ms_split+'SUN01_20141101T'+timestr+'.50ms.slfcal.dbkg.ms'


debkg=True
prep=False
tofits=True


if debkg:
    slfcalms = dir_ms_split+'SUN01_20141101T'+timestr+'.50ms.slfcal.dbkg.ms'
else:
    slfcalms = dir_ms_split+'SUN01_20141101T'+timestr+'.50ms.slfcal.ms'   


tofits=True


print 'Script for calibrating 2014 Nov 1 data --- '+structure_id_slfcalms
print ''



structure_id = 'rU01-test2'
pol='LL'
refantenna='ea04'
antennas='0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,26'

fitsprefix='./slfcal/fits2/'
if not os.path.exists(fitsprefix):
    os.mkdir(fitsprefix)
imgprefix='slfcal/'
if not os.path.exists(imgprefix):
    os.mkdir(imgprefix)
imgprefix='slfcal/'+structure_id+'/'
print 'image files output to --- '+imgprefix
if not os.path.exists(imgprefix):
    os.mkdir(imgprefix)	    
imgprefix=imgprefix+structure_id+'_'

# pb = progressbar(idx_tim.size, "*")    
# for ii in range(0,len(idx_selec[0])):
# for ii in [ 4528, 4690, 4855, 5021, 5184,
#        5347, 5505, 5660, 5814]:	
# for ii in range(0,idx_tim.size):
for ii in range(0,100):
	t_int=0.05
	f_int=2.0 #MHz
	t0=tim[idx_selec[1][ii]]#-0.5*t_int
	t1=t0+t_int
	t0str=qa.time(qa.quantity(t0,'s'),prec=9)[0]
	t1str=qa.time(qa.quantity(t1,'s'),prec=9)[0]
	timestr=qa.time(qa.quantity(t0,'s'),form=['no_dms','fits'],prec=9)[0].translate(None, ':')
	timeran=t0str+'~'+t1str	

	# for chn in range(2,62):
	# for chn in range(32,35):
	chn=idx_selec[0][ii]
	print timestr,chn	
	f0=freq[chn]/1.e6
	f1=f0+f_int
	freqstr='{:d}MHz'.format(int(round(f0)))
	print imgprefix+timestr+'_'+freqstr+'.slfcal'
	slfcal_img = imgprefix+timestr+'_'+freqstr+'.slfcal'
	# slfcal_img_local = imgprefix+timestr+'_'+freqstr+'.local.slfcal'.format(chn)
	default('clean')
	vis=slfcalms
	imagename=slfcal_img
	spw='0:{:d}'.format(chn)
	timerange=timeran
	mode='mfs'
	# interpolation = 'nearest'
	imagermode='csclean'
	weighting='briggs'
	gain=0.05
	psfmode='hogbom'
	# mask = ['circle [ [ 147pix, 181pix ], 28pix]']
	# mask = 'box [ [ 119pix , 153pix] , [175pix, 211pix ] ]'
	# mask='rU01.mask0'
	########## global region ##########
	imsize=[256,256]
	cell = ['10arcsec', '10arcsec']
	phasecenter = 'J2000 14h26m22.7351 -14d29m29.801'
	########## global region ##########
	# imsize=[512,512]
	# cell = ['5.0arcsec', '5.0arcsec']
	# phasecenter = 'J2000 14h26m22.7351 -14d29m29.801'	
	########## local region	##########
	# imsize=[128,128]
	# cell = ['5.0arcsec', '5.0arcsec']
	# phasecenter = 'J2000 14h26m59.250 -14d35m44.681'
	stokes=pol
	uvtaper=True
	outertaper=['50arcsec']
	uvrange=''
	niter=1000
	npercycle=100
	usescratch=False
	clean()	
	os.system('rm -rf '+slfcal_img+'.flux')
	os.system('rm -rf '+slfcal_img+'.mask')
	os.system('rm -rf '+slfcal_img+'.model')
	os.system('rm -rf '+slfcal_img+'.psf')
	os.system('rm -rf '+slfcal_img+'.residual')    

	# if os.path.exists(slfcal_img+'.image'):	
	# 	ia.open(slfcal_img+'.image')
	# 	ROI='centerbox[[14h26m59.250,-14d35m44.681], [384.0arcsec, 384.0arcsec]]'
	# 	ia.fromimage(outfile=slfcal_img_local+'.image', infile=slfcal_img+'.image', region=ROI, overwrite=true)  
	# 	ia.close()  
# print("--- %s seconds ---" % (time.time() - start_timestamp))
	''' ----- step 4 ----- '''
	if tofits and os.path.exists(slfcal_img+'.image'):
		ephem=vla_prep.read_horizons(ephemfile=ephemfile)
		reftime=[timeran]*1
		helio=vla_prep.ephem_to_helio(msinfo=msinfofile,ephem=ephem,reftime=reftime)
		imagenames=[slfcal_img]
		imagefile=[img+'.image' for img in imagenames]
		vlafits=[fitsprefix+img.split('/')[-1]+'.image.fits' for img in imagenames]
		for files in vlafits:
			if os.path.exists(files):
				os.system('rm -rf '+files)
		vla_prep.imreg(imagefile=imagefile,fitsfile=vlafits,helio=helio,toTb=F,scl100=True,blc=[50,100],trc=[82,132])
	# pb.progress(ii)
print("--- %s seconds ---" % (time.time() - start_timestamp))
	# os.system('rm -rf slfcal/rU01/*'+'.flux')
	# os.system('rm -rf slfcal/rU01/*'+'.mask')
	# os.system('rm -rf slfcal/rU01/*'+'.model')
	# os.system('rm -rf slfcal/rU01/*'+'.psf')
	# os.system('rm -rf slfcal/rU01/*'+'.residual')    

