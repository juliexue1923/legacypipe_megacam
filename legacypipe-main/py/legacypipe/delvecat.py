import os
import numpy as np
from legacypipe.ps1cat import HealpixedCatalog

class DELVECatalog(HealpixedCatalog):
    def __init__(self, file_prefix=None, indexing=None, ccdwcs=None, **kwargs):
        self.delvedir = os.getenv('DELVE_CAT_DIR')
        if self.delvedir is None:
            raise ValueError('You must have the DELVE_CAT_DIR environment variable set to point to healpixed DELVE catalogs')
        if indexing is None:
            indexing = os.getenv('DELVE_CAT_SCHEME', 'ring')
        if not indexing in ['nested', 'ring']:
            raise ValueError('Supported values for the DELVE_CAT_SCHEME environment variable or healpix indexing scheme are "nested" or "ring"')
        if file_prefix is None:
            file_prefix = os.getenv('DELVE_CAT_PREFIX', 'delve')
        #
        fnpattern = os.path.join(self.delvedir, file_prefix + '-%(hp)05d.fits')
        super(DELVECatalog, self).__init__(fnpattern, indexing=indexing, **kwargs)
        if ccdwcs is not None:
            self.ccdwcs = ccdwcs
