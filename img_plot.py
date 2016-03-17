#!/usr/bin/env python

#################################
### make a tile plot for NGC 891
#################################
import numpy as np
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pylab as plt
import jhy_math
import pdb
import pyfits
import scipy.ndimage
from astropy.io import ascii
from matplotlib.patches import Ellipse
matplotlib.rcParams.update({'font.size':12})
from astropy.stats import sigma_clip

###
### find RA&Dec matched region for different bands
###
def define_region(hdr,bar=False,CD=False):
	cen = np.loadtxt('image_center.dat',skiprows=1,dtype=float)
	xcen = 894
	ycen = 1350

	if CD == False:
		dra_scale = hdr['CDELT1']
		ddec = hdr['CDELT2']
	if CD == True:
		dra_scale = hdr['CD1_1']#/3600.
		ddec = hdr['CD2_2']#/3600.
	pscale = 0.0003458	# in arcmin

	dx = 8.5/60.	# 8 arcmin FOV
	xmin = round(xcen - abs(dx/pscale))
	xmax = round(xcen + abs(dx/pscale))
	ymin = round(ycen - dx/pscale)
	ymax = round(ycen + dx/pscale)

	return xmin,xmax,ymin,ymax

###
### plot ellipse photometry aperture
### 
def plot_ellipse(data,hdr,bandname,xcen=-999,ycen=-999):
	band_list = {'MIPS24':0,'PACS70':1,'PACS160':2,'SPIRE250':3,\
					'SPIRE350':4,'SPIRE500':5}
	pix_scale = np.array([0.000345833334658,0.000388889,0.000791667,0.00125,0.00173611,0.0025])

	plt.plot(xcen,ycen,linestyle='none',marker='+',mew=1,color='black')

	a_spire500 = 49
	a_b = 25/50.
	xcen_ellipse = xcen#+37.51
	ycen_ellipse = ycen#+24.22
	b_spire500 = a_spire500 * a_b
	pscale = pix_scale[band_list[bandname]]

	a_major = pix_scale[band_list['SPIRE500']]/pscale*a_spire500
	b_minor = pix_scale[band_list['SPIRE500']]/pscale*b_spire500

	Ny = len(data[:,0])+1
	Nx = len(data[0,:])+1
	x = np.arange(Nx*50)/50.
	y = np.arange(Ny*50)/50.
	PA = 0#np.radians(180. - 69)	# -72deg
	N_pixel = 0
	sum_in_ellipse = 0
	err_in_ellipse = []
	all_sky = []

	xx = x - xcen_ellipse
	y_upper = np.sqrt(1.-xx**2./a_major**2.)*b_minor
	y_lower = -1*np.sqrt(1.-xx**2./a_major**2.)*b_minor

	x_rot_upper = np.cos(PA)*xx - np.sin(PA)*y_upper
	y_rot_upper = np.sin(PA)*xx + np.cos(PA)*y_upper
	x_rot_lower = np.cos(PA)*xx - np.sin(PA)*y_lower
	y_rot_lower = np.sin(PA)*xx + np.cos(PA)*y_lower

	x_final_upper = x_rot_upper + xcen_ellipse
	y_final_upper = y_rot_upper + ycen_ellipse
	x_final_lower = x_rot_lower + xcen_ellipse
	y_final_lower = y_rot_lower + ycen_ellipse

	plt.plot(x_final_upper,y_final_upper,color='black',lw=1.5)
	plt.plot(x_final_lower,y_final_lower,color='black',lw=1.5)

	d_ngc891_Ibata = 9.73*1000	# kpc, Ibata+2009
	hz_star_kpc = 1.44 # 1.44 kpc Ibata+2009
	hz_star_deg = np.degrees(np.arctan(hz_star_kpc/d_ngc891_Ibata))
	hz_star_pix = hz_star_deg/pix_scale[0]
	Ny, Nx = data.shape
	xarr = np.arange(Nx)
	yarr = np.arange(Ny)
	slope_deg = 0
	slope = np.tan(np.radians(slope_deg))
	a_major = slope*(xarr-xcen)+ycen
	disk_upper = a_major + hz_star_pix/np.sin(np.radians(90+slope_deg))
	disk_lower = a_major - hz_star_pix/np.sin(np.radians(90+slope_deg))
	sub = x_final_upper > 0
	xsub1 = np.min(x_final_upper[sub])
	xsub2 = np.max(x_final_upper[sub])
	plt.plot(xarr[xsub1:xsub2],disk_upper[xsub1:xsub2],linestyle=':',color='black',lw=2)
	plt.plot(xarr[xsub1:xsub2],disk_lower[xsub1:xsub2],linestyle=':',color='black',lw=2)

	return 

###
### plot arrow pointing where the flare is
###
def plot_arrow(xmin,xmax,ymin,ymax):
	x_arrow = (xmax-xmin)*4.4/10.+xmin
	y_arrow = (ymax-ymin)*20./100.+ymin
	dxy_arrow = (xmax-xmin)*1./9.
	hw_arrow = (xmax-xmin)*1./20.
	hl_arrow = (ymax-ymin)*1./15.
	plt.arrow(x_arrow,y_arrow,dxy_arrow,-0*dxy_arrow,fc='black',ec='black',\
				linewidth=3,head_width=hw_arrow,head_length=hl_arrow)
	return 

###
### plot NSWE direction
###
def plot_direction(xmin,xmax,ymin,ymax):
	dx = (xmax-xmin)
	dy = (ymax-ymin)
	x_arrow = dx*1.5/10.+xmin
	y_arrow = dy*1.5/10.+ymin
	dxy_arrow = dx*1./35.
	x_N0 = x_arrow
	y_N0 = y_arrow
	dx_N1 = -dxy_arrow*np.tan(np.radians(69))
	dy_N1 = dxy_arrow
	hw_arrow = (xmax-xmin)*1./50.
	hl_arrow = (ymax-ymin)*1./50.
	plt.arrow(x_N0,y_N0,dx_N1,dy_N1,fc='black',ec='black',\
				linewidth=1,head_width=hw_arrow,head_length=hl_arrow)
	plt.text(x_N0+dx_N1*1.4,y_N0+dy_N1*1.3,'N',rotation=69,fontsize=10)
	dx_E1 = -dxy_arrow
	dy_E1 = -dxy_arrow*np.tan(np.radians(69))
	plt.arrow(x_N0,y_N0,dx_E1,dy_E1,fc='black',ec='black',\
				linewidth=1,head_width=hw_arrow,head_length=hl_arrow)
	plt.text(x_N0+dx_E1*0.9,y_N0+dy_E1*1.5,'E',rotation=69,fontsize=10)

	return 
###
### Body
###

###
### read an image
###
data_mips24,hdr_mips24 = pyfits.getdata('NGC891_mips24_masked_smooth_rot_flat.fits',0,header=True)
xmin_mips24,xmax_mips24,ymin_mips24,ymax_mips24 = \
												define_region(hdr_mips24,CD=True)

pix_scale = np.array([0.000345833334658,0.000388889,0.000791667,0.00125,0.00173611,0.0025])

###
### plot image
###
matplotlib.rcParams.update({'font.size':15})
plt.clf()
fig = plt.figure(figsize=(8,5))
ax = fig.add_subplot(111)
ax.set_title('MIPS 24$\\mu$m',fontsize=12)
sub_lowSN = data_mips24 < 1e-9
data_mips24[sub_lowSN] = -999.
dat = np.log10(data_mips24/(pix_scale[0]*np.pi/180.)**2.)
plt.imshow(dat,origin='lower',cmap='spectral',vmin=3,vmax=7.7)
cb = plt.colorbar(orientation='vertical')
cb.set_label('log SB (Jy/sr)',fontsize=12)

pscale = 0.0003458*60	# in arcmin
xtick = np.arange(-8,9,2)
ytick = np.arange(-8,9,2)
ax.set_xticks(xtick/pscale+894)
ax.set_xticklabels(map(str,xtick))
ax.set_yticks(ytick/pscale+1350)
ax.set_yticklabels(map(str,ytick))

ax.set_aspect('equal')
plt.xlim(xmin_mips24,xmax_mips24)
plt.ylim(ymin_mips24,ymax_mips24-150)
plot_arrow(xmin_mips24,xmax_mips24,ymin_mips24,ymax_mips24)
plot_ellipse(data_mips24,hdr_mips24,'MIPS24',\
				xcen=(xmin_mips24+xmax_mips24)/2,\
				ycen=(ymin_mips24+ymax_mips24)/2,)
plt.xlabel('$\\Delta$X (arcmin)')
plt.ylabel('$\\Delta$Y (arcmin)')
plot_direction(xmin_mips24,xmax_mips24,ymin_mips24,ymax_mips24)
plt.savefig('img_plot.pdf')

