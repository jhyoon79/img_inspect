###
### This is to inspect artifact in the image
###
import numpy as np
import matplotlib.pylab as plt
import pdb
import pyfits
import jhy_math
from matplotlib.backends.backend_pdf import PdfPages
fig_pdf = PdfPages('check_artifact.pdf')

###
### define region in the sky that I'm interested
###
def define_region(hdr):
	racen_ngc891 = 35.639224# + 0.0020
	deccen_ngc891 = 42.349146# - 0.0030
	Nx = hdr['NAXIS1']
	Ny = hdr['NAXIS2']
	xcen = hdr['CRPIX1']
	ycen = hdr['CRPIX2']
	racen = hdr['CRVAL1']
	deccen = hdr['CRVAL2']
#	PA = hdr['CROTA2']
	try:
		dra_scale = hdr['CDELT1']
		ddec = hdr['CDELT2']
	except KeyError:
		dra_scale = hdr['CD1_1']+hdr['CD2_1']
		ddec = hdr['CD2_2']+hdr['CD1_2']
	dra = jhy_math.angsep2radiff(dra_scale,deccen)
	ra = (np.arange(Nx)+1 - xcen) * dra + racen
	dec = (np.arange(Ny)+1 - ycen) * ddec + deccen
	ra_diff = abs(ra - racen_ngc891)
	dec_diff = abs(dec - deccen_ngc891)
	xcen_ngc891 = list(ra_diff).index(min(ra_diff)) + 1
	ycen_ngc891 = list(dec_diff).index(min(dec_diff)) + 1

	return xcen_ngc891,ycen_ngc891

###
### conduct photometry in a rectangular region
###
def rect_phot(data,err,hdr):

	err = data*0
	xcen,ycen = define_region(hdr)

	total_flux = []
	pscale = 0.000791667*60	# pacs160 pix scale in arcmin

	z_min = 0	# initial scale height in arcmin
	z_max = 3/pscale
	h_z = 0.1
	dz = h_z/pscale	
	a_major = 80#14/pscale
	x1 = xcen - a_major
	x2 = xcen + a_major

	summed_flux_all = []
	summed_err_all = []
	z_mean_all = []
	pix_mean_all = []
	z_step = np.copy(dz)
	for j in [-1,1]:
		count = 0
		for i in range(30):
			z_range = count*j*dz#z_step
			if abs(z_range) >= z_max: break
			if i==0: pre_y2,z_step = ycen,dz
			y1 = pre_y2#np.copy(ycen + z_range)
			y2 = y1 + z_step*j
			if j < 0:
				data_sub = data[y2:y1,x1:x2]
				err_sub = err[y2:y1,x1:x2]
			if j > 0:
				data_sub = data[y1:y2,x1:x2]
				err_sub = err[y1:y2,x1:x2]
	
			Ny,Nx = data_sub.shape
			area = 1.#(Nx*Ny*(pscale/60.)**2.)/(180/np.pi)**2.	# in steradian
			summed_flux = np.sum(data_sub)/area
			summed_err = np.sqrt(np.sum(err_sub**2.+\
										(data_sub*0.10)**2.))/area
			pix_mean = ((y1+y2)/2.-ycen)
			z_mean = ((y1+y2)/2.-ycen)*pscale	# in arcmin
			count += 1

			summed_flux_all.append(summed_flux)
			summed_err_all.append(summed_err)
			z_mean_all.append(z_mean)
			pix_mean_all.append(pix_mean)
			pre_y2 = y2
			z_step = dz
			print summed_flux,pix_mean
	return np.array(summed_flux_all),np.array(summed_err_all),np.array(z_mean_all),np.array(pix_mean_all)

###
### Read fits images
###
color,hdr = pyfits.getdata('color_70_160_rot_flat.fits',0,header=True)
color_emurph01,hdr_emurph01 = pyfits.getdata('color_70_160_emurph01_rot_flat.fits',0,header=True)
color_emurph01_new,hdr_emurph01_new = pyfits.getdata('color_70_160_emurph01_new_rot_flat.fits',0,header=True)

###
### get the photometry of the region of interest
###
z_color,z_color_err,z_mean,pix_mean = rect_phot(color,color,hdr)
z_color_emurph01,z_color_err_emurph01,z_mean_emurph01,pix_mean_emurph01 = rect_phot(color_emurph01,color_emurph01,hdr_emurph01)
z_color_emurph01_new,z_color_err_emurph01_new,z_mean_emurph01_new,pix_mean_emurph01_new = rect_phot(color_emurph01_new,color_emurph01_new,hdr_emurph01_new)

###
### plot for the inspection of the artifact
###
plt.clf()
plt.plot(pix_mean,z_color,marker='o',linestyle='none',color='black')
plt.plot(pix_mean_emurph01,z_color_emurph01,marker='d',linestyle='none',color='blue')
plt.plot(pix_mean_emurph01_new,z_color_emurph01_new,marker='^',linestyle='none',color='red')
plt.xlabel('pix')
plt.ylabel('color')
plt.ylim(-50,500)
fig_pdf.savefig()
fig_pdf.close()
