import os
import numpy as np
from legacypipe.ps1cat import HealpixedCatalog

class SMSSCatalog(HealpixedCatalog):

    delveband = dict(u = "umag")

    def __init__(self, file_prefix=None, indexing=None, ccdwcs=None, **kwargs):
        self.smssdir = os.getenv('SMSS_CAT_DIR')
        if self.smssdir is None:
            raise ValueError('You must have the SMSS_CAT_DIR environment variable set to point to healpixed SkyMapper catalogs')
        if indexing is None:
            indexing = os.getenv('SMSS_CAT_SCHEME', 'nested')
        if not indexing in ['nested', 'ring']:
            raise ValueError('Supported values for the SMSS_CAT_SCHEME environment variable or healpix indexing scheme are "nested" or "ring"')
        if file_prefix is None:
            file_prefix = os.getenv('SMSS_CAT_PREFIX', 'delve')
        #
        fnpattern = os.path.join(self.delvedir, file_prefix + '-%(hp)05d.fits')
        super(SMSSCatalog, self).__init__(fnpattern, indexing=indexing, **kwargs)
        if ccdwcs is not None:
            self.ccdwcs = ccdwcs

    def get_stars(self,magrange=None,band='g'):
        """Return the set of SkyMapper stars on a given CCD with well-measured gri
        magnitudes. Optionally trim the stars to a desired r-band magnitude
        range.
        """
        cat = self.get_catalog_in_wcs(self.ccdwcs)
        print('Found {} good SkyMapper stars'.format(len(cat)))
        if magrange is not None:
            keep = np.where((getattr(cat,SMSSCatalog.smssband[band])>magrange[0])*
                            (getattr(cat,SMSSCatalog.smssband[band])<magrange[1]))[0]
            cat = cat[keep]
            print('Trimming to {} stars with {}=[{},{}]'.
                  format(len(cat),band,magrange[0],magrange[1]))
        return cat

    def get_catalog_in_wcs(self, wcs, step=100., margin=10):
        # Grid the CCD in pixel space
        W,H = wcs.get_width(), wcs.get_height()
        xx,yy = np.meshgrid(
            np.linspace(1-margin, W+margin, 2+int((W+2*margin)/step)),
            np.linspace(1-margin, H+margin, 2+int((H+2*margin)/step)))
        # Convert to RA,Dec and then to unique healpixes
        ra,dec = wcs.pixelxy2radec(xx.ravel(), yy.ravel())
        healpixes = set()
        for r,d in zip(ra,dec):
            healpixes.add(self.healpix_for_radec(r, d))
        # Read catalog in those healpixes
        cat = self.get_healpix_catalogs(healpixes)
        # Cut to sources actually within the CCD.
        _,xx,yy = wcs.radec2pixelxy(cat.ra, cat.dec)
        cat.x = xx
        cat.y = yy
        onccd = np.flatnonzero((xx >= 1.-margin) * (xx <= W+margin) *
                               (yy >= 1.-margin) * (yy <= H+margin))
        cat.cut(onccd)
        return cat
