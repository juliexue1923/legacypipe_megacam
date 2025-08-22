import os
import numpy as np
from legacypipe.ps1cat import HealpixedCatalog

class GaiaCatalog(HealpixedCatalog):
    def __init__(self, file_prefix=None, indexing=None, ccdwcs=None, **kwargs):
        self.gaiadir = os.getenv('GAIA_CAT_DIR')
        if self.gaiadir is None:
            raise ValueError('You must have the GAIA_CAT_DIR environment variable set to point to healpixed Gaia catalogs')
        if indexing is None:
            indexing = os.getenv('GAIA_CAT_SCHEME', 'ring')
        if not indexing in ['nested', 'ring']:
            raise ValueError('Supported values for the GAIA_CAT_SCHEME environment variable or healpix indexing scheme are "nested" or "ring"')
        if file_prefix is None:
            file_prefix = os.getenv('GAIA_CAT_PREFIX', 'chunk')
        #
        fnpattern = os.path.join(self.gaiadir, file_prefix + '-%(hp)05d.fits')
        super(GaiaCatalog, self).__init__(fnpattern, indexing=indexing, **kwargs)
        if ccdwcs is not None:
            self.ccdwcs = ccdwcs

    def get_stars(self,magrange=None,band='r'):
        """Return the set of gaia stars on a given CCD with well-measured grz
        magnitudes. Optionally trim the stars to a desired r-band magnitude
        range.
        """
        cat = self.get_catalog_in_wcs(self.ccdwcs)
        print('Found {} good Gaia stars'.format(len(cat)))
        if magrange is not None:
            keep = np.where((cat.median[:,ps1cat.ps1band[band]]>magrange[0])*
                            (cat.median[:,ps1cat.ps1band[band]]<magrange[1]))[0]
            cat = cat[keep]
            print('Trimming to {} stars with {}=[{},{}]'.
                  format(len(cat),band,magrange[0],magrange[1]))
        return cat

    def get_catalog_radec_box(self, ralo, rahi, declo, dechi):
        import numpy as np

        wrap = False
        if rahi < ralo:
            # wrap-around?
            rahi += 360.
            wrap = True

        # Prepare RA,Dec grid to pick up overlapping healpixes
        rr,dd = np.meshgrid(np.linspace(ralo,  rahi,  2+int(( rahi- ralo)/0.1)),
                            np.linspace(declo, dechi, 2+int((dechi-declo)/0.1)))
        healpixes = set()
        for r,d in zip(rr.ravel(), dd.ravel()):
            healpixes.add(self.healpix_for_radec(r, d))
        # Read catalog in those healpixes
        cat = self.get_healpix_catalogs(healpixes)
        #print('Read', len(cat), 'Gaia catalog entries.  RA range', cat.ra.min(), cat.ra.max(),
        #      'Dec range', cat.dec.min(), cat.dec.max())
        cat.cut((cat.dec >= declo) * (cat.dec <= dechi))
        if wrap:
            cat.cut(np.logical_or(cat.ra <= ralo, cat.ra >= (rahi - 360.)))
        else:
            cat.cut((cat.ra  >= ralo ) * (cat.ra  <= rahi))
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
        #self.indexing = 'nested'
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

    def get_healpix_catalog(self, healpix):
        from astrometry.util.fits import fits_table
        fname = self.fnpattern % dict(hp=healpix)
        #fname = fname.replace('chunk', 'gaia')
        print('Reading', fname)
        return fits_table(fname)

    @staticmethod
    def catalog_nantozero(gaia):
        gaia.pmra = nantozero(gaia.pmra)
        gaia.pmdec = nantozero(gaia.pmdec)
        gaia.parallax = nantozero(gaia.parallax)
        return gaia

def nantozero(x):
    import numpy as np
    x = x.copy()
    x[np.logical_not(np.isfinite(x))] = 0.
    return x
