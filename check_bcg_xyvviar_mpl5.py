"""
This script plots individual BCGs.

Usage: python check_bcg_xyvviar_mpl5.py bcg_number
"""
import os
import sys
import numpy as np
import pyfits
from matplotlib import pylab as plt
import sep
from PIL import Image

threshold = 0.25  # Threshold for r band image masking
bcg_cata = np.loadtxt('./bcg_list_alexie_5')  # n by 5 where 0: plateid, 1: designid, 4: dist_kpc
k = np.int32(sys.argv[1])

plateid, designid, dist_kpc = int(bcg_cata[k,0]), int(bcg_cata[k,1]), bcg_cata[k,4]*1000
print "plateid, designid, dist_kpc=", plateid, designid, dist_kpc

# Read the velocity map.
hdu = pyfits.open('./data/manga-'+str(plateid)+'-'+str(designid)+'-MAPS-VOR10-GAU-MILESHC.fits.gz')
xc, yc, xm, ym, dx, dy = hdu['STELLAR_VEL'].header['CRPIX1'], hdu['STELLAR_VEL'].header['CRPIX2'], \
                         hdu['STELLAR_VEL'].header['NAXIS1'], hdu['STELLAR_VEL'].header['NAXIS2'], \
                         hdu['STELLAR_VEL'].header['PC1_1'], hdu['STELLAR_VEL'].header['PC2_2']

# Read the velocity map excluding the satellites.
satellite_mask_hdu = pyfits.open('./sex/'+str(plateid)+'-'+str(designid)+'-MAPS-VOR10-GAU-MILESHC-nosatelite.fits')
assert satellite_mask_hdu[0].header['CRPIX1'] == xc, satellite_mask_hdu[0].header['CRPIX2'] == yc
satellite_mask = satellite_mask_hdu[0].data[:xm,:ym]
masked_image = np.ma.array(hdu['STELLAR_VEL'].data, mask=np.logical_or(hdu['STELLAR_VEL_MASK'].data,\
        satellite_mask))
masked_ivar = np.ma.array(hdu['STELLAR_VEL_IVAR'].data, mask=np.logical_or(hdu['STELLAR_VEL_MASK'].data,\
        satellite_mask))

# Normalize the pixels within 3 arcsec to the center to be 0.
xyrange = np.int32(3/.5)
vavgin3 = np.mean(masked_image[xm//2-xyrange:xm//2+xyrange,ym//2-xyrange:ym//2+xyrange])
print vavgin3
masked_image -= vavgin3

# Query the dababase for logcube files and apply a r band cut.
spectro = 'manga-'+str(plateid)+'-'+str(designid)+'-LOGCUBE.fits.gz'
# Type in userid and password (not allowed to distribute).
if not os.path.exists('./data/'+spectro):
    os.system(\
            'wget --user=xxx --password=xxx https://data.sdss.org/sas/mangawork/manga/spectro/redux/MPL-5/'+\
                    str(plateid)+'/stack/'+spectro+'; mv '+spectro+' data/')
hdusn = pyfits.open('./data/'+spectro)
rimg_sn = np.ma.array(hdusn['RIMG'].data, mask=hdusn['RIMG'].data<threshold)

# Check pixel info
assert dx == -0.000138889, dy == 0.000138889
assert hdu['STELLAR_VEL'].header['CRPIX1'] == hdusn['RIMG'].header['CRPIX1']
assert hdu['STELLAR_VEL'].header['CRPIX2'] == hdusn['RIMG'].header['CRPIX2']
assert hdu['STELLAR_VEL'].header['NAXIS1'] == hdusn['RIMG'].header['NAXIS1']
assert hdu['STELLAR_VEL'].header['NAXIS2'] == hdusn['RIMG'].header['NAXIS2']
assert hdu['STELLAR_VEL'].header['PC1_1'] == hdusn['RIMG'].header['CD1_1']
assert hdu['STELLAR_VEL'].header['PC2_2'] == hdusn['RIMG'].header['CD2_2']

# Extract pixel info for rebinning
x = np.arange(-xc, xm-xc)*dx*dist_kpc/180.*np.pi
y = np.arange(-yc, ym-yc)*dy*dist_kpc/180.*np.pi
for i in xrange(len(x)):
    for j in xrange(len(y)):
        if i == j == 0:
            a = [[x[i], y[j], masked_image[i,j], masked_ivar[i,j], rimg_sn[i,j]]]
        else:
            a.append([x[i], y[j], masked_image[i,j], masked_ivar[i,j], rimg_sn[i,j]])

a = np.ma.array(a)
nonnan_ind = np.where(np.logical_and(~np.isnan(a[:,2]),a[:,3]>1e-6))[0]
a = a[nonnan_ind]
nbins, resol = 100, 1
image_rebin, image_rebin_err, image_rebin_count = np.zeros((nbins,nbins)), np.zeros((nbins,nbins)), np.zeros((nbins,nbins))
rimg_rebin, rimg_rebin_count = np.ma.zeros((nbins,nbins)), np.ma.zeros((nbins,nbins))
reference = nbins // 2
x_range = (np.arange(nbins)-reference)*resol
y_range = (np.arange(nbins)-reference)*resol

xi = np.int32((a[:,0]-x_range[0])/resol)
yi = np.int32((a[:,1]-y_range[0])/resol)
for i in xrange(len(xi)):
    rimg_rebin_count[xi[i],yi[i]] += 1
    rimg_rebin[xi[i],yi[i]] += a[i,4]
    image_rebin_count[xi[i],yi[i]] += 1
    image_rebin[xi[i],yi[i]] += a[i,2]
    image_rebin_err[xi[i],yi[i]] += 1/a[i,3]

image_rebin /= image_rebin_count
image_rebin_err = np.sqrt(image_rebin_err) / image_rebin_count

nbins1d, dr_small, dr_large, pmass, perror2 = 10, 0, 30, a[:,2], 1/a[:,3]
feat_profile = np.zeros((nbins1d,4))
r = np.sqrt(a[:,0]*a[:,0]+a[:,1]*a[:,1])
dx = (dr_large - dr_small) * 1. / nbins1d
ind_array = np.floor((r - dr_small) / dx)
xa = np.linspace(dr_small, dr_large, nbins1d+1)
feat_profile[:,0] = xa[:-1] + np.diff(xa)/2
for i in xrange(nbins1d):
    feat_profile[i,1] = np.sum(ind_array==i)
    feat_profile[i,2] = np.sum(pmass[ind_array==i])/np.sum(ind_array==i)
    feat_profile[i,3] = np.sqrt(np.sum(perror2[ind_array==i]))/np.sum(ind_array==i)

# Add the segmentation map.
hdu = pyfits.open('./crop/'+str(plateid)+'-'+str(designid)+'-MAPS-VOR10-GAU-MILESHC-r-noconv.fits')
x_range0 = -(np.arange(hdu[0].header['NAXIS1'])-hdu[0].header['CRPIX1'])*hdu[0].header['CD1_1']*dist_kpc*np.pi/180
y_range0 = (np.arange(hdu[0].header['NAXIS2'])-hdu[0].header['CRPIX2'])*hdu[0].header['CD2_2']*dist_kpc*np.pi/180
data0 = hdu[0].data.byteswap().newbyteorder()
bkg = sep.Background(data0, bw=8, bh=8)
bkg.subfrom(data0)
rms = bkg.globalrms
objects, seg = sep.extract(data0, 1.0*rms, deblend_cont=0.01, segmentation_map=True)
hdu.close()

#-----
fig = plt.figure(figsize=(16,8))
ax = fig.add_subplot(231)
im = ax.imshow(image_rebin, origin='lower', vmin=-100, vmax=100, interpolation='nearest', extent=[x_range[0], x_range[-1], y_range[0], y_range[-1]])
ax.set_title('Velocities')
ax.set_xlim(-30,30)
ax.set_ylim(-30,30)
fig.colorbar(im)

ax = fig.add_subplot(234)
im = ax.imshow(np.sqrt(image_rebin_err), origin='lower', interpolation='nearest', extent=[x_range[0], x_range[-1], y_range[0], y_range[-1]])
ax.set_title('Errors on velocities')
ax.set_xlim(-30,30)
ax.set_ylim(-30,30)

fig.colorbar(im)
ax = fig.add_subplot(232)
im = ax.imshow(rimg_rebin/rimg_rebin_count, origin='lower', interpolation='nearest', extent=[x_range[0], x_range[-1], y_range[0], y_range[-1]])
ax.set_title('rimg')
ax.set_xlim(-30,30)
ax.set_ylim(-30,30)

fig.colorbar(im)
ax = fig.add_subplot(235)
ax.errorbar(feat_profile[:,0], feat_profile[:,2], yerr=feat_profile[:,3])
ax.set_title('1d profile')

ax = fig.add_subplot(233)
ax.imshow(hdu[0].data, interpolation='nearest', extent=[x_range0[0], x_range0[-1], y_range0[0], y_range0[-1]])
ax.set_xlim(-30,30)
ax.set_ylim(-30,30)
fig.colorbar(im)

ax = fig.add_subplot(236)
ax.imshow(seg, interpolation='nearest', extent=[x_range0[0], x_range0[-1], y_range0[0], y_range0[-1]])
ax.set_xlim(-30,30)
ax.set_ylim(-30,30)
fig.colorbar(im)

# Add the jpeg image
imjpeg = Image.open('./jpg/'+str(plateid)+'-'+str(designid)+'-MAPS-VOR10-GAU-MILESHC.png')
fig.figimage(imjpeg, xo=1450, yo=0)
plt.savefig(str(plateid)+'-'+str(designid)+'.png')

#if os.path.exists(spectro):
#    os.system('rm '+spectro)

