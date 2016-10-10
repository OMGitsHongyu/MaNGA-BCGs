import os
import sys
import numpy as np
import pyfits
from matplotlib import pylab as plt

def reject_outliers(data, m = 2.):
    d = np.abs(data - np.median(data))
    mdev = np.median(d)
    s = d/mdev if mdev else 0.
    return s < m


def stack_bcgs(bcg_cata, spirals, cut=0.25, weighted=0):
    """
    Stack BCGs in a directed/weighted fashion.

    Parameters
    ----------
    bcg_data : (nbcgs,5) array
        BCG info (plateid, designid, gmass, smass, dist_kpc).
    spirals : (nspirals,) list
        spirals plateid-designid in string.
    """

    p = 0
    for k in xrange(bcg_cata.shape[0]):

        plateid, designid, dist_kpc = int(bcg_cata[k,0]), int(bcg_cata[k,1]), bcg_cata[k,4]*1000
        if str(plateid)+'-'+str(designid) in spirals:
            print k, p, plateid, designid
            try:
                hdu = pyfits.open('./data/manga-'+str(plateid)+'-'+str(designid)+'-MAPS-VOR10-GAU-MILESHC.fits.gz')
                xc, yc, xm, ym, dx, dy = hdu['STELLAR_VEL'].header['CRPIX1'], hdu['STELLAR_VEL'].header['CRPIX2'], \
                                         hdu['STELLAR_VEL'].header['NAXIS1'], hdu['STELLAR_VEL'].header['NAXIS2'], \
                                         hdu['STELLAR_VEL'].header['PC1_1'], hdu['STELLAR_VEL'].header['PC2_2']

                hdusn = pyfits.open('./data/manga-'+str(plateid)+'-'+str(designid)+'-LOGCUBE.fits.gz')
                satellite_mask_hdu = pyfits.open('./sex/'+str(plateid)+'-'+str(designid)+'-MAPS-VOR10-GAU-MILESHC-nosatelite.fits')

                assert dx == -0.000138889, dy == 0.000138889
                assert hdu['STELLAR_VEL'].header['CRPIX1'] == hdusn['RIMG'].header['CRPIX1'] == satellite_mask_hdu[0].header['CRPIX1']
                assert hdu['STELLAR_VEL'].header['CRPIX2'] == hdusn['RIMG'].header['CRPIX2'] == satellite_mask_hdu[0].header['CRPIX2']
                assert hdu['STELLAR_VEL'].header['NAXIS1'] == hdusn['RIMG'].header['NAXIS1'] == satellite_mask_hdu[0].header['NAXIS1']-1
                assert hdu['STELLAR_VEL'].header['NAXIS2'] == hdusn['RIMG'].header['NAXIS2'] == satellite_mask_hdu[0].header['NAXIS2']-1
                assert hdu['STELLAR_VEL'].header['PC1_1'] == hdusn['RIMG'].header['CD1_1']
                assert abs(hdu['STELLAR_VEL'].header['PC1_1'] - satellite_mask_hdu[0].header['CD1_1']) < 1e-6
                assert hdu['STELLAR_VEL'].header['PC2_2'] == hdusn['RIMG'].header['CD2_2']
                assert abs(hdu['STELLAR_VEL'].header['PC2_2'] - satellite_mask_hdu[0].header['CD2_2']) < 1e-6

                masked_image = np.ma.array(hdu['STELLAR_VEL'].data, mask=np.logical_or(np.logical_or(\
                        hdu['STELLAR_VEL_MASK'].data, hdusn['RIMG'].data<cut), satellite_mask_hdu[0].data[0:xm,0:ym]))
                masked_ivar = np.ma.array(hdu['STELLAR_VEL_IVAR'].data, mask=np.logical_or(np.logical_or(\
                        hdu['STELLAR_VEL_MASK'].data, hdusn['RIMG'].data<cut), satellite_mask_hdu[0].data[0:xm,0:ym]))
#        masked_image = np.ma.array(hdu['STELLAR_VEL'].data, mask=hdu['STELLAR_VEL_MASK'].data)
#        masked_ivar = np.ma.array(hdu['STELLAR_VEL_IVAR'].data, mask=hdu['STELLAR_VEL_MASK'].data)

                xyrange = np.int32(3/.5)
                vavgin3 = np.mean(masked_image[xm//2-xyrange:xm//2+xyrange,ym//2-xyrange:ym//2+xyrange])
                print vavgin3
                masked_image -= vavgin3

                hdu.close()
                hdusn.close()
                satellite_mask_hdu.close()
                x = np.arange(-xc, xm-xc)*dx*dist_kpc/180.*np.pi
                y = np.arange(-yc, ym-yc)*dy*dist_kpc/180.*np.pi
                for i in xrange(len(x)):
                    for j in xrange(len(y)):
#                a = np.ma.array([x[i], y[j], masked_image[i,j], masked_ivar[i,j]])
                        if p == i == j == 0:
                            a = [[x[i], y[j], masked_image[i,j], masked_ivar[i,j]]]
                        else:
                            a.append([x[i], y[j], masked_image[i,j], masked_ivar[i,j]])
                p += 1
            except IOError as e:
                print e
                pass

    a = np.ma.array(a)
    nonnan_ind = np.where(np.logical_and(~np.isnan(a[:,2]),a[:,3]>1e-6))[0]
    a = a[nonnan_ind]
    nbins, resol = 100, 1
    image_rebin, image_rebin_err, image_rebin_count = np.zeros((nbins,nbins)), np.zeros((nbins,nbins)), np.zeros((nbins,nbins))
    reference = nbins // 2
    x_range = (np.arange(nbins)-reference)*resol
    y_range = (np.arange(nbins)-reference)*resol

    xi = np.int32((a[:,0]-x_range[0])/resol)
    yi = np.int32((a[:,1]-y_range[0])/resol)
    for i in xrange(len(xi)):
        image_rebin_count[xi[i],yi[i]] += 1
        image_rebin[xi[i],yi[i]] += a[i,2]
        image_rebin_err[xi[i],yi[i]] += 1/a[i,3]

    image_rebin /= image_rebin_count
    image_rebin_err = np.sqrt(image_rebin_err) / image_rebin_count

    nbins1d, dr_small, dr_large, pmass, perror2, pivar = 50, 0, 30, a[:,2], 1/a[:,3], a[:,3]
    feat_profile = np.ma.zeros((nbins1d,4))
    r = np.sqrt(a[:,0]*a[:,0]+a[:,1]*a[:,1])
    dx = (dr_large - dr_small) * 1. / nbins1d
    ind_array = np.floor((r - dr_small) / dx)
    xa = np.linspace(dr_small, dr_large, nbins1d+1)
    feat_profile[:,0] = xa[:-1] + np.diff(xa)/2
    print ind_array
    print pmass
    if weighted:
        for i in xrange(nbins1d):
            feat_profile[i,1] = np.sum(ind_array==i)
            if feat_profile[i,1] > 0:
                pmass_ind = reject_outliers(pmass[ind_array==i])
                pmass_weighed = pmass * pivar
                if pmass_weighed[ind_array==i].shape[0] > 1:
                    feat_profile[i,2] = np.sum(pmass_weighed[ind_array==i][pmass_ind])/np.sum(pivar[ind_array==i][pmass_ind])
                    feat_profile[i,3] = np.sqrt(np.sum(perror2[ind_array==i]))/np.sum(ind_array==i)
                else:
                    feat_profile[i,2] = np.sum(pmass_weighed[ind_array==i])/np.sum(pivar[ind_array==i])
                    feat_profile[i,3] = np.sqrt(np.sum(perror2[ind_array==i]))/np.sum(ind_array==i)
            else:
                feat_profile[i,2] = np.nan
                feat_profile[i,3] = np.nan
    else:
        for i in xrange(nbins1d):
            feat_profile[i,1] = np.sum(ind_array==i)
            feat_profile[i,2] = np.sum(pmass[ind_array==i])/np.sum(ind_array==i)
            if feat_profile[i,1] > 0:
                feat_profile[i,2] = np.median(pmass[ind_array==i])
                feat_profile[i,3] = np.sqrt(np.sum(perror2[ind_array==i]))/np.sum(ind_array==i)
            else:
                feat_profile[i,2] = np.nan
                feat_profile[i,3] = np.nan

    fig = plt.figure()
    ax = fig.add_subplot(221)
    im = ax.imshow(image_rebin, origin='lower', interpolation='nearest', vmin=-10, vmax=60, extent=[x_range[0], x_range[-1], y_range[0], y_range[-1]])
    ax.set_title('Velocities')
    fig.colorbar(im)
    ax = fig.add_subplot(223)
    im = ax.imshow(np.sqrt(image_rebin_err), origin='lower', interpolation='nearest', vmin=0, vmax=5, extent=[x_range[0], x_range[-1], y_range[0], y_range[-1]])
    ax.set_title('Errors on velocities')
    fig.colorbar(im)
    ax = fig.add_subplot(222)
    im = ax.imshow(image_rebin_count, origin='lower', interpolation='nearest', vmin=0, vmax=200, extent=[x_range[0], x_range[-1], y_range[0], y_range[-1]])
    ax.set_title('counts')
    fig.colorbar(im)
    ax = fig.add_subplot(224)
    ax.errorbar(feat_profile[:,0], feat_profile[:,2], yerr=feat_profile[:,3])
    ax.set_title('1d profile')
    ax.set_xlim((0,20))
    ax.set_ylim((-100,100))
    plt.show()

if __name__ == '__main__':

    bcg_cata = np.loadtxt('./bcg_list_alexie_5')
    spirals = ['8082-12704',
               '8135-12701',
               '8254-6103',
               '8318-6102',
               '8459-3703',
               '8482-9101',
               '8482-12701']

    stack_bcgs(bcg_cata, spirals)
