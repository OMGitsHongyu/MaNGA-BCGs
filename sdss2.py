import os
import pyfits as fits
import glob
import numpy as np
import string
import scipy.ndimage.filters as ss
import fileinput as f
import time
import pdb
import sep
from matplotlib import pyplot as plt
from scipy.stats import mode


def download_jpgs(infn):
    """
    Download jpgs of galaxies.
    """

    with open(infn, 'r') as f:
        filelines = f.read()
    filestring = filelines.splitlines()

    if os.path.exists('./jpg/'):
        os.system('rm -rf jpg/')
    os.system('mkdir jpg/')

    for i in xrange(len(filestring)):
        fileinfo = filestring[i].split(',')
        imwidth, pxscl, ra, dec = 400., 0.396, fileinfo[1], fileinfo[2]
        url = 'http://skyservice.pha.jhu.edu/DR10/ImgCutout/getjpeg.aspx?ra='+str(ra)+\
              '&dec='+str(dec)+'&scale='+str(pxscl)+'&opt=G&width='+str(imwidth)+'&height='+str(imwidth)
        imfilename = fileinfo[0]+'.png'
        os.system('wget -O '+ imfilename + " '" + url + "'")
        os.system('mv *png ./jpg/')


def extract_sources(infn):
    """
    Extract unrelated satellites.
    """

    with open(infn, 'r') as f:
        filelines = f.read()
    filestring = filelines.splitlines()

    for i in xrange(len(filestring)):
        fileinfo = filestring[i].split(',')
        hdu = fits.open('./crop/'+fileinfo[0]+'-r-noconv.fits')
        x_range = (np.arange(hdu[0].header['NAXIS1'])-hdu[0].header['CRPIX1'])*hdu[0].header['CD1_1']*np.float64(fileinfo[5])*np.pi/180
        y_range = (np.arange(hdu[0].header['NAXIS2'])-hdu[0].header['CRPIX2'])*hdu[0].header['CD2_2']*np.float64(fileinfo[5])*np.pi/180
        data = hdu[0].data.byteswap().newbyteorder()
        objects = sep.extract(data, 3.0)
        mask = np.zeros_like(data)
        for i in xrange(len(objects)):
            mask[objects['ymin'][i]:objects['ymax'][i],objects['xmin'][i]:objects['xmax'][i]] = i+1

        data = hdu[0].data.byteswap().newbyteorder()
        bkg = sep.Background(data, bw=8, bh=8)
        bkg.subfrom(data)

        rms = bkg.globalrms
        objects, seg = sep.extract(data, 2.0*rms, deblend_cont=0.01, segmentation_map=True)

        segmode = mode(seg[int(hdu[0].header['CRPIX1']-3):int(hdu[0].header['CRPIX1']+3),int(hdu[0].header['CRPIX2']-3)\
                        :int(hdu[0].header['CRPIX2']+3)], axis=None)[0][0]
        satellite_mask = np.zeros_like(seg)
        satellite_mask[seg!=segmode] = 1

        if os.path.exists('./a.fits'):
            os.system('rm a.fits')
        fits.writeto('./a.fits',satellite_mask,header=hdu[0].header)
        x_range = (np.arange(hdu[0].header['NAXIS1'])-hdu[0].header['CRPIX1'])*hdu[0].header['CD1_1']*3600
        y_range = (np.arange(hdu[0].header['NAXIS2'])-hdu[0].header['CRPIX2'])*hdu[0].header['CD2_2']*3600

        hduorigin = fits.open('../data/manga-'+fileinfo[0]+'.fits.gz')
        npix = hduorigin['STELLAR_VEL'].header['NAXIS1'] + 1
        if os.path.exists('./coadd.fits'):
            os.system('rm coadd.fits')
        os.system('swarp ./a.fits -center '+fileinfo[1]+','+fileinfo[2]+' -image_size '+str(npix)+' -subtract_back N -resampling_type nearest -pixel_scale 0.5')
        if not os.path.exists('./sex/'):
            os.system('mkdir sex/')
        os.system('mv coadd.fits sex/'+fileinfo[0]+'-nosatelite.fits')
        hdu2 = fits.open('./sex/'+fileinfo[0]+'-nosatelite.fits')

        hduorigin['STELLAR_VEL'].header['CRPIX1'] == hdu2[0].header['CRPIX1']
        hduorigin['STELLAR_VEL'].header['CRPIX2'] == hdu2[0].header['CRPIX2']
        abs(np.float32(hduorigin['STELLAR_VEL'].header['PC1_1']) - np.float32(hdu2[0].header['CD1_1'])) < np.finfo(np.float32).eps
        abs(np.float32(hduorigin['STELLAR_VEL'].header['PC2_2']) - np.float32(hdu2[0].header['CD2_2'])) < np.finfo(np.float32).eps

        hdu.close()
        hdu2.close()
        hduorigin.close()

if __name__ == '__main__':
    download_jpgs('./input.txt')
    extract_sources('./input.txt')
