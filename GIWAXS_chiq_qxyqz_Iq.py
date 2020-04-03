import numpy as np
import glob
import os
import pyFAI
import fabio
import pygix
#import kb_gen
from scipy.ndimage.filters import median_filter
import matplotlib.pyplot as plt
import matplotlib.colors  as colors
import matplotlib
# import pygix.plotting as pp
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pandas as pd


def make_wedge_zero(image, qxy):
	"""
	sets pixel values in missing wedge in qxyqz plot to zero
	pygix sets it to some stupid values and source code is compiled
	ideally a mask would be created only once
	"""
	import numpy as np
	x = int(np.argmin(abs(qxy)))
	# print x
	for row in np.arange(image.shape[0]):
		c = image[row,x]
		j = 1
		try:
			while image[row,x-j] == image[row,x]: j +=1
		except: pass
		image[row,x-j:x+j] = .0
	return image


def crop_qxyqz(intensity, qxy, qz, qxyrange=[-10.,10.], qzrange=[0.,30.]):
	"""
	crops qxyqz image
	units are 1/nm by default
	returns cropped intensity, cropped qxy, cropped qz
	"""
	import numpy as np
	qxy0 = np.argmin(np.abs(qxy-qxyrange[0]))
	qxy1 = np.argmin(np.abs(qxy-qxyrange[1]))
	qz0 = np.argmin(np.abs(qz-qzrange[0]))
	qz1 = np.argmin(np.abs(qz-qzrange[1]))
	return intensity[qz0:qz1, qxy0:qxy1], qxy[qxy0:qxy1], qz[qz0:qz1]




ponifile = './poni.poni'
parentdir = './'
resultdir = os.path.join(parentdir, 'results')
qbins = 1000
incident_angle = 2.5
sample_orientation = 3




def replace_at_index1(tup, ix, val):
     lst = list(tup)
     for i in range(0,len(ix)):
         lst[ix[i]] = val[i]
     return tuple(lst)




tifs = sorted(glob.glob(os.path.join(parentdir, '*.tif')))[::-1]
print(tifs)
ai = pyFAI.AzimuthalIntegrator()
ai.load(ponifile)
r = np.zeros((qbins+1, len(tifs)+1), dtype = object)
try: 
	os.makedirs(os.path.join(resultdir, 'png'))
except: 
	pass
try:
	os.makedirs(os.path.join(resultdir, 'png_qxyqz'))
except:
	pass
try:
	os.makedirs(os.path.join(resultdir, 'png_chiq'))
except:
	pass
for tif, i in zip(tifs, list(range(len(tifs)))):
	filename = tif.split('/')[-1][:-4]
	img = fabio.open(tif).data
	print(('Now reading ', tif))
	q, intensity = ai.integrate1d(img, qbins, polarization_factor=-.95)
	if i == 0:
		r[0,0] = 'q'
		r[1:,0] = q*.1
	r[0,i+1] = os.path.basename(tif)
	r[1:,i+1] = intensity

	df = pd.DataFrame(r[1:], columns=r[0])
	df['q'] = df['q'].apply(lambda x: x*10)
	df.set_index(df['q'], inplace=True)
	df.drop(['q'], axis=1, inplace=True)
	plt.imshow(np.log(img))
	plt.savefig(os.path.join(parentdir, 'results', 'png', os.path.basename(tif).replace('.tif','')+'.png'), bbox_inches='tight')
	plt.close()

	pg = pygix.Transform()
	pg.load(ponifile)
	pg.incident_angle = incident_angle
	pg.sample_orientation = sample_orientation
	
	# qxyxz plot:
	intensity, qxy, qz = pg.transform_reciprocal(img, npt=(qbins,qbins), polarization_factor=-.95, method='nearest')
	with open('{}/png_qxyqz/{}_imageData.dat'.format(resultdir, filename), 'w') as f:
		np.savetxt(f, intensity)
	
	# default method is "splitpix", which for some reason fills the wedge with interpolated pixels. Function is defined in C:/Anaconda2/Lib/site-packages/pyFAI/ext/splitpixelFull.pyd, which is already compiled, so cannot be fixed.
	# Methods 'np' and 'cython' don't fill wedge, but returned array is much too big and mostly black, so valuable data fills only a few pixels. ip_range and op_range are defined to customize the range, but don't work
	# two options: a) postprocess to fill in wedge with black pixels or b) fix 'np' or 'cython' methods
	matplotlib.rcdefaults() # pygix messes around with tex in matplotlib which crashes subsequent figures; this restores the default 
	fig, axarr = plt.subplots(1, 1, figsize=(18,10))
	intensity = make_wedge_zero(intensity, qxy)
	# intensity, qxy, qz = crop_qxyqz(intensity, qxy, qz, qxyrange=[0,qxy.max()], qzrange=[0.,qz.max()]) # in 1/nm
	# pp.implot(intensity, qxy, qz, mode='rsm', newfig=False, show=True, filename=resultpath+sample+'/png_qxyqz/'+f.replace('.tiff','')+'_qxyqz.png') # not good due to tex issue
	
	intensity[intensity <= 0.0] = 'NaN'
	maskIntensity = np.ma.masked_where(np.isnan(intensity),intensity, copy=True)
	colornorm=colors.LogNorm(vmin=maskIntensity.min(), vmax=maskIntensity.max())
	plot1=axarr.imshow(maskIntensity, extent=np.array([qxy.min(), qxy.max(), qz.min(), qz.max()]), aspect = 1., norm=colornorm, cmap='inferno', origin='lower', interpolation='none') # 0.1 is for conversion from 1/nm to 1/A
	axarr.set_title(os.path.basename(tif))
	axarr.set_xlabel('q [1/nm]', size=24)
	axarr.set_ylabel('q [1/nm]', size=24)
	axarr.tick_params(direction='in', which='major', length=7, labelsize=20)
	cbar_coord = replace_at_index1(make_axes_locatable(axarr).get_position(), [0,2], [0.83,0.02])
	cbar_ax = fig.add_axes(cbar_coord)
	cbar = fig.colorbar(plot1, cax=cbar_ax)
	cbar_ax.tick_params(direction='out', which='major', length=7, labelsize=20)
	cbar_ax.tick_params(direction='out', which='minor', length=4)
	# plt.show()
	plt.savefig(os.path.join(resultdir, 'png_qxyqz', os.path.basename(tif).replace('.tif','')+'_qxyqz.pdf'), bbox_inches='tight')
	plt.close()

	intensity, q, chi = pg.transform_polar(np.log(img), npt=(qbins, qbins), q_range=(5,30), chi_range=(-90, 90), polarization_factor=-.95, method='nearest')
	intensity = np.rot90(np.rot90(np.rot90(make_wedge_zero(np.rot90(intensity), chi))))
	intensity = median_filter(intensity, size = 5)
	# intensity, qxy, qz = crop_qxyqz(intensity, q, chi, qxyrange=[0,q.max()], qzrange=[-45.,chi.max()])
	intensity[intensity <= 0.0] = 'NaN'
	maskIntensity = np.ma.masked_where(np.isnan(intensity),intensity, copy=True)
	colornorm=colors.LogNorm(vmin=maskIntensity.min(), vmax=maskIntensity.max())
	

	fig, axarr = plt.subplots(2, 1, figsize=(18,20))
	plot2 = axarr[0].imshow(maskIntensity, extent = [q.min(), q.max(), chi.min(), chi.max()], aspect = "auto", cmap='inferno', origin='lower', interpolation='none')
	axarr[0].set_title(os.path.basename(tif))
	axarr[0].set_xlabel('q [1/nm]', size=24)
	axarr[0].set_ylabel('q [1/nm]', size=24)
	axarr[0].tick_params(direction='in', which='major', length=7, labelsize=20)
	
	cbar_coord2 = replace_at_index1(make_axes_locatable(axarr[0]).get_position(), [0,2], [0.92,0.02])
	
	cbar_ax2 = fig.add_axes(cbar_coord2)
	cbar2 = fig.colorbar(plot2,cax=cbar_ax2)
	cbar_ax2.tick_params(direction='out', which='major', length=7, labelsize=20)
	cbar_ax2.tick_params(direction='out', which='minor', length=4, labelsize=20)
	
	plot3 = axarr[1].plot(df[os.path.basename(tif)])
	axarr[1].set_title(os.path.basename(tif))
	axarr[1].set_xlim(q.min(), q.max())
	axarr[1].set_xlabel('q [1/nm]', size=24)
	axarr[1].set_ylabel('Intensity', size=24)
	axarr[1].tick_params(direction='in', which='major', length=7, labelsize=20)

	# plt.show()
	plt.savefig(os.path.join(resultdir, 'png_chiq', os.path.basename(tif).replace('.tif','')+'_chiq.pdf'), bbox_inches='tight')
	plt.close()

	
	


np.savetxt(os.path.join(parentdir, 'Iq.csv'), r, delimiter = ",", fmt = "%s")


