"""
This file is modified from David V. Stark @ IPMU's code. It is used to query
sdss database to preprocess r-band image masks and jpegs and rebin them in the
designed format.
"""
import subprocess
import re
import os
import time
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

def parse_input(inputfile):
    """
    Parse the input file for query_sdss.py.
    """
    f = open(inputfile,'r')
    for line in f:
        if line.rfind('=') > -1 and line.rfind('#') != 0:

            s = string.split(line,'=')
            arg1 = string.strip(s[0])
            arg2 = string.split(s[1],'#')
            arg2 = string.strip(arg2[0])

            if arg1 == 'inputfile':
                inputfile = arg2
                print inputfile
            elif arg1 == 'infofile':
                infofile = arg2
                print infofile
            elif arg1 == 'tablesdir':
                tablesdir = arg2
                print tablesdir
            elif arg1 == 'cropdir':
                cropdir = arg2
                print cropdir
            elif arg1 == 'tempdir':
                tempdir = arg2
                print tempdir
            elif arg1 == 'calibdir':
                calibdir = arg2
                print calibdir
            elif arg1 == 'bands':
                bands = arg2
                print bands
            elif arg1 == 'rundownload':
                rundownload = arg2
                print rundownload
            elif arg1 == 'runsextractor':
                runsextractor = 1
                print runsextractor
            else:
                print 'Warning: input file line, '+ arg1 + ', unrecognized...skipping...'

    out = [inputfile, infofile, tablesdir, cropdir, tempdir, calibdir, bands, rundownload\
            , runsextractor]
    return out


def getsdssframeinfo(photofieldfile,field):
    """
    Get frame info.
    """
    photofield = fits.open(photofieldfile)
    info = photofield[1].data
    sel = np.where(info['field'] == field)
    dark_variance = info[sel]['dark_variance']
    gain = info[sel]['gain']
    nmgypercount = info[sel]['nMgyPerCount']
    sky = info[sel]['sky']
    # psfs=info[sel]['psf_width']
    psfs = info[sel]['PSF_SIGMA1_2G']
    score = info[sel]['SCORE']
    out = [dark_variance, gain, nmgypercount, sky, psfs, score]

    return out


def getavgsdsserrinfo(imfilelist):

    #create lists to hold everything
    pxsclarr = []
    gainarr = []
    nmgycntarr = []
    skyarr = []
    dvarr = []
    texparr = []
    psfarr = []

    for i in range(len(imfilelist)):

        im = fits.open(imfilelist[i])
        imhdr = im[0].header
        pxsclx = imhdr['CD1_1']
        pxscly = imhdr['CD1_2']
        pxscl = pow(pxsclx*pxsclx+pxscly*pxscly, 0.5) * 3600.
        gain = imhdr['GAIN']
        nmgycnt = imhdr['NMGYCNT']
        sky = imhdr['SKY']
        dv = imhdr['DV']
        texp = float(imhdr['EXPTIME'])
        psf = imhdr['PSF']

        print sky

        pxsclarr += [pxscl]
        texparr += [texp]
        gainarr += [gain]
        nmgycntarr += [nmgycnt]
        skyarr += [sky]
        dvarr += [dv]
        psfarr += [psf]

    pxscl = np.mean(pxsclarr)
    texp = np.mean(texparr)
    gain = np.mean(gainarr)
    nmgycnt = np.mean(nmgycntarr)
    sky = np.mean(skyarr)
    dv = np.mean(dvarr)
    maxpsf = np.amax(psfarr)
    minpsf = np.amin(psfarr)
    meanpsf = np.median(psfarr)

    out = [pxscl, gain, nmgycnt, sky, dv, texp, maxpsf, minpsf, meanpsf]

    return out


def copy_header(image,outpath):

    if os.path.exists(outpath) == 0:
        os.system('mkdir ' + outpath)

    im = fits.open(image)

    hdu = fits.PrimaryHDU()
    hdu.header = im[0].header

    # Get rid of any paths in the input file name
    outname = string.split(image,'/')
    outname = outname[-1]

    hdu.writeto(outpath+outname, clobber=1)


# Definition of sdss,UKIDSS,GALEX bands
def download_sdss(setupfile, sdssbands, runconv=0):
    setup = parse_input(setupfile)
    inputfile = setup[0]
    infofile = setup[1]
    tablesdir = setup[2]
    cropdir = setup[3]
    tempdir = setup[4]
    calibdir = setup[5]
    bands = setup[6]
    rundownload = setup[7]
    runsextractor = setup[8]

    print 'inputfile: ', inputfile
    print 'tablesdir: ', tablesdir
    print 'cropdir: ', cropdir
    print 'tempdir: ', tempdir
    print 'calibdir: ', calibdir
    print 'bands: ', bands

    print ''
    print 'reading input file'
    print ''

    # Make bands a list
    bands = string.split(bands,',')

    # Test that directories exist, make them if not
    if not os.path.exists(tablesdir):
        os.system('mkdir ' + tablesdir)
    if not os.path.exists(cropdir):
        os.system('mkdir ' + cropdir)
    if not os.path.exists(tempdir):
        os.system('mkdir ' + tempdir)
    if not os.path.exists(calibdir):
        os.system('mkdir ' + calibdir)

    # Open input file
    i = open(inputfile, 'r')

    # Open output info file
    o=open(infofile, 'w')

    pxscl = 0.396         # arcseconds per pixel

    # Cycle through targets now
    for line in i:

        s = string.split(line,',')
        print s
        name = string.strip(s[0])
        ra = float(s[1])
        dec = float(s[2])
        radius = float(s[3])
        print ra, dec

        worstpsf = 0. #will keep running account of worst psf as we dl ukidss and SDSS

        print ''
        print 'Searching for image data for '+name
        print ''

        range_deg = radius / 3600.   #needs to be in degrees
        range_px = radius / pxscl    #needs to be in pixels

        runsdss = ('u' in bands) or ('g' in bands) or ('r' in bands) or ('i' in bands) or ('z' in bands)

        if runsdss:

            print ''
            print 'Querying SDSS database'
            print ''

            queryURL = 'python ./sqlcl.py -s http://skyserver.sdss3.org/dr8/en/tools/search/x_sql.asp -f csv -q "SELECT p.run,p.rerun,p.camcol,p.field FROM PhotoObj as p WHERE (p.ra > ' + str(ra-range_deg) + ' AND p.ra < ' + str(ra+range_deg) + ') AND (p.dec > ' + str(dec-range_deg)+ ' AND p.dec < ' + str(dec+range_deg)+ ')"'
            # command line sql query of sloan database, returns output to
            # "query" which I then extract. extracted elements is a "tuple"
            # with 2 elements.  The first elements is the SQL output, so I
            # extact that and split it into 1 line per galaxy.  Each line
            # is now a string with run,rerun,camcol,field

            query = subprocess.Popen(queryURL, shell=True, stdout=subprocess.PIPE)
            sqlresult = query.communicate()
            sqlresult, junk = sqlresult
            sqlresult = re.split(', |\n', sqlresult)
            sqlresult = sqlresult[1:len(sqlresult)-1]

            # Remove duplicates
            sqlresult = list(set(sqlresult))

            sdssimnumber = -1
            print sqlresult
            for element in sqlresult:
                sdssimnumber = sdssimnumber + 1
                rrcf = element.split(',')
                run = rrcf[0]
                rerun = rrcf[1]
                camcol = rrcf[2]
                field = rrcf[3]

                # This just makes naming conventions easier
                if int(run) < 1000:
                    fff = '000'
                elif int(run) < 100:
                    fff = '0000'
                elif int(run) < 10:
                    fff = '00000'
                else:
                    fff = '00'
                if int(field) < 100:
                    aaa = '00'
                else:
                    aaa = '0'

                # Now write this to the output info file
                space = '   '
                outstring = name+space + str(ra) + space + str(dec) + space + fff \
                         + run + space + rerun + space + camcol + space + aaa \
                         + field + '\n'
                o.write(outstring)

                # Download photofield file
                photofieldfile = 'photoField-' + fff + run + '-' + camcol + '.fits'
                url = 'http://data.sdss3.org/sas/dr8/groups/boss/photoObj/301/' + \
                        run + '/' + photofieldfile
                os.system('wget ' + url)
                os.system('mv ' + photofieldfile + ' ' + tempdir)

                # Have run into issues where the file was not found, which causes this
                # to crash.  See if file exists, if not, then set "score" to zero so it
                # will skip over it
                photofieldfile = tempdir + photofieldfile

                if os.path.exists(photofieldfile):

                # Open photofield file and extract frame info

                    dark_variance, gain, nmgypercount, sky, psfs, score = \
                            getsdssframeinfo(photofieldfile, int(field))

                    print ''
                    print 'image score: ', score
                    print ''

                    # Don't need photofield file anymore, erase
                    os.system('rm ' + photofieldfile)

                else:
                    print 'cannot find file! (' + photofieldfile + ')'

                    score = 0

                if score < 0.6:
                    print 'bad score...skipping'

                for band in bands:

                    if (band in sdssbands) and (score >= 0.6): # Added requirement that "score" be photometric quality

                        print 'good score, downloading:'

                        # Download image file
                        tempfile = 'frame-' + band + '-' + fff + run + '-' + camcol\
                                + '-' + aaa + field + '.fits'
                        url = 'http://mirror.sdss3.org/sas/dr9/boss/photoObj/frames/301/'\
                                + run + '/' + camcol + '/' + tempfile + '.bz2'
                        os.system('wget ' + url)
                        os.system('bunzip2 ' + tempfile + '.bz2')
                        os.system('mv ' + tempfile + ' ' + tempdir)

                        # Read in first fits extension and save
                        im = fits.open(tempdir + tempfile)
                        hdr = im[0].header
                        data = im[0].data

                        # Add frame info to header
                        hdr.update('field', field, '')
                        hdr.update('dv', float(dark_variance[0, np.where(sdssbands==band)]), 'dark variance')
                        hdr.update('gain', float(gain[0, np.where(sdssbands==band)]), 'gain')
                        hdr.update('nmgycnt', float(nmgypercount[0, np.where(sdssbands==band)]), 'nanomaggies per count')
                        hdr.update('sky', float(sky[0, np.where(sdssbands==band)]), 'sky level')
                        hdr.update('psf', float(psfs[0, np.where(sdssbands==band)]), 'Effective PSF (arcseconds)')

                        fits.writeto(calibdir + band + '-' + str(sdssimnumber) + '-noconv.fits', data, hdr, clobber=1)

                        # Remove frame file now
                        os.system('rm ' + tempdir + tempfile)

                        #also create version convolved to the worst PSF
                        #sel=np.where(sdssbands == bands)
                        #thispsf=psfs[sel]
                        #filtersz=np.sqrt(worstpsf*worstpsf-thispsf*thispsf)
                        #convim=ss.gaussian_filter(data,filtersz/pxscl)
                        #fits.writeto(calibdir+band+'-'+str(sdssimnumber)+'-conv.fits',convim,hdr,clobber=1)

                        maxframepsf = np.amax(psfs)
                        if maxframepsf > worstpsf:
                            worstpsf = maxframepsf

        for band in bands:

            # Find all images in this band
            files = glob.glob(calibdir + band + '*-noconv.fits')
            print files
            if len(files) > 0:

                np.savetxt('filelist.txt', files, fmt='%s')

                # Account for different pixel sizes, whether to background subtract
                if band in sdssbands:
                    bgsub = 'N'
                    desired_pxscl = 0.5

                swarpfile = 'stripe82_alt.swarp'
                swarpcommand ='swarp @filelist.txt -c ' + swarpfile + ' -center '\
                              + str(ra) + ',' + str(dec) + ' -image_size ' + str(\
                              2*radius/desired_pxscl) + ' -subtract_back ' + bgsub\
                              + ' -pixel_scale ' + str(desired_pxscl)

                print 'swarp', ra, dec, 2*radius/desired_pxscl, bgsub, desired_pxscl
                print swarpcommand
                os.system(swarpcommand)

                finalim = cropdir + name + '-' + band + '-noconv.fits'
                os.system('mv coadd.fits ' + finalim)

                # Extract uncertainty and psf info

                if band in sdssbands:

                    errinfo = getavgsdsserrinfo(files)

                    # Write this info the cropped final image
                    cropim = fits.open(finalim, mode='update')
                    cropim[0].header.update('OPXSCL', errinfo[0], 'Original pixel scale (arcsec pixel^-1)')
                    cropim[0].header.update('GAIN', errinfo[1], 'Mean gain')
                    cropim[0].header.update('NMGYCNT', errinfo[2], 'Mean nanomaggies per count')
                    cropim[0].header.update('SKY', errinfo[3], 'Mean sky level')
                    cropim[0].header.update('DV', errinfo[4], 'Mean dark variance')
                    cropim[0].header.update('TEXP', errinfo[5], 'Mean exposure time (s)')
                    cropim[0].header.update('MAXPSF', errinfo[6], 'Maximum seeing included')
                    cropim[0].header.update('MINPSF', errinfo[7], 'Minimum seeing included')
                    cropim[0].header.update('MEANPSF', errinfo[8], 'Median seeing of all frames')

                    cropim.flush()

                print 'MAX PSF', worstpsf


    # Now convolve SDSS and UKIDSS data to worst SDSS resolution

    # Copy over the header for the fits images
    files = glob.glob(cropdir + name + '-*-noconv.fits') + glob.glob(cropdir + name + '*-conv.fits') +\
            glob.glob(cropdir + name + '*-galexconv.fits')
    outpath = tablesdir + '/headers/'

    for eachfile in files:
        copy_header(eachfile, outpath)

    o.close()



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
        imfilename = fileinfo[0] + '.png'
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
        x_range = (np.arange(hdu[0].header['NAXIS1'])-hdu[0].header['CRPIX1'])*hdu[0].header['CD1_1']\
                *np.float64(fileinfo[5])*np.pi/180
        y_range = (np.arange(hdu[0].header['NAXIS2'])-hdu[0].header['CRPIX2'])*hdu[0].header['CD2_2']\
                *np.float64(fileinfo[5])*np.pi/180
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
#    sdssbands = np.array(['u','g','r','i','z'])
#    setupfile = './setup.txt'
#    download_sdss(setupfile, sdssbands)

    inputfn = './input.txt'
    download_jpgs(inputfn)
    extract_sources(inputfn)
