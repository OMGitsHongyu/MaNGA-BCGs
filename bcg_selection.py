import os
import sys
import pyfits
import numpy as np
from astropy.cosmology import WMAP7
import functools
import glob
import sdss

def select_BCGs(cataname, bcgs=1, download=1, version='mpl5', outfn='bcg_list_alexie_5'):
    """
    Select BCGs from MaNGA, params are set below.

    cataname : string
        catalog name.
    """
    # Select BCGs according to a catalog fits file.
    bcghdu = pyfits.open(cataname)
    if bcgs == 1:
        a0 = functools.reduce(np.logical_and, (bcghdu[1].data['DRP_MATCH']==1,\
                                               bcghdu[1].data['GRP_MHALO_LEST']>13.5,\
                                               bcghdu[1].data['BCG_VIS']==1,\
                                               bcghdu[1].data['VIS_BCG_FLAG']!=0))
        a1 = functools.reduce(np.logical_or, (bcghdu[1].data['VIS_CLUSTER_FLAG']==1,\
                                              bcghdu[1].data['VIS_CLUSTER_FLAG']==2,\
                                              bcghdu[1].data['VIS_CLUSTER_FLAG']==5))
        cata = np.where(np.logical_and(a0,a1))[0]
    else:
        cata = np.where(bcghdu[1].data['PLATE']!=-1)[0]

    print "Number of BCGs", len(cata)

    # Save BCGs info to a list of lists.
    # Format: plate ifudsgn grp_mhalo_lest stellarmass z
    for j in range(len(cata)):
        plateid, designid, gmass, smass, z = bcghdu[1].data['PLATE'][cata[j]],\
                                             bcghdu[1].data['IFUDSGN'][cata[j]],\
                                             bcghdu[1].data['GRP_MHALO_LEST'][cata[j]],\
                                             bcghdu[1].data['STELLARMASS'][cata[j]],\
                                             bcghdu[1].data['Z'][cata[j]]
        print j, plateid, designid, gmass, smass, WMAP7.comoving_distance(z).value
        if j == 0:
            bcg_list = [[plateid, designid, gmass, smass, WMAP7.comoving_distance(z).value]]
        else:
            bcg_list.append([plateid, designid, gmass, smass, WMAP7.comoving_distance(z).value])

        # Download
        if download:
            if version == 'mpl4':
                os.system(\
                        'wget --user=sdss --password=2.5-meters https://data.sdss.org/sas/'+\
                                'mangawork/manga/sandbox/mangadap/MPL-4/full/'+\
                                str(plateid)+'/'+str(designid)+'/manga-'+str(plateid)+'-'+str(designid)\
                                +'-LOGCUBE_MAPS-NONE-023.fits.gz')
            if version == 'mpl5':
                os.system(\
                        'wget --user=sdss --password=2.5-meters https://data.sdss.org/sas/'+\
                                'mangawork/manga/spectro/analysis/MPL-5/VOR10-GAU-MILESHC/'+\
                                str(plateid)+'/'+str(designid)+'/manga-'+str(plateid)+'-'+str(designid)\
                                +'-MAPS-VOR10-GAU-MILESHC.fits.gz')
    bcghdu.close()

    # Save outputs
    if outfn is not None:
        if os.path.exists(outfn):
            os.system('rm '+outfn)
            np.savetxt(outfn, np.array(bcg_list), fmt='%d %d %.4f %.5f %.5f')

    return bcg_list


def write_input(infn, datadir='./', outfn='input.txt'):
    """
    Write the input files for the source extractor script.
    """
    bcg_cata = np.loadtxt(infn)
    if os.path.exists(outfn):
        os.system('rm '+outfn)
    with open(outfn, 'w') as f:
        for k in xrange(len(bcg_cata)):
            try:
                plateid, designid, dist_kpc = int(bcg_cata[k,0]), int(bcg_cata[k,1]), bcg_cata[k,4]*1000
                datafn = glob.glob(datadir+'*'+str(plateid)+'-'+str(designid)+'*.fits.gz')
                hdu = pyfits.open(datafn[0])
                ra, dec = hdu['STELLAR_VEL'].header['CRVAL1'], hdu['STELLAR_VEL'].header['CRVAL2']
                fnstart, fnend = datafn[0].find(str(plateid)+'-'+str(designid)), datafn[0].find('.fits.gz')
                f.write(datafn[0][fnstart:fnend]+','+str(ra)+','+str(dec)+',25,1,'+str(dist_kpc)+'\n')
                hdu.close()
                assert len(datafn) == 1, 'might be two versions, check!!!'
            # Some files may be not appropriately downloaded
            except IndexError:
                pass


if __name__ == '__main__':
    bcg_list = select_BCGs('../../data/combined_yang_cat_drpv1_5_1.fits')
    write_input('bcg_list_alexie_5')
