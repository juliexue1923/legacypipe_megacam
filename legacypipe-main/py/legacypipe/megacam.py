import os
from datetime import datetime
import fitsio

from legacypipe.image import LegacySurveyImage

class MegaCamImage(LegacySurveyImage):

    zp0 = dict(
        g = 25.0,
        u = 25.0,
        r = 25.0,
    )
       
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.scamp_wcs = None

    def read_image_primary_header(self, **kwargs):
        '''
        Reads the FITS data (HDU 1) header from self.imgfn.

        Returns
        -------
        primary_header : fitsio header
            The FITS header
        '''
       # if self._first_header is not None:
       #     return self._first_header
        self._primary_header = self.read_image_fits()[1].read_header()
        return self._primary_header

    def get_ccdname(self, primhdr, hdr):
        if 'EXTNAME' in hdr:
            return hdr['EXTNAME'].strip().upper()
        elif 'MEXTNO' in hdr:
            ext_num = hdr['MEXTNO'].strip().upper()
            return 'IM' + ext_num
        else: 
            return ''

    def get_expnum(self, primhdr):
        date = primhdr['DATE-OBS']
        # DATE-OBS= '2022-01-29T04:21:01' / UT date (yyyy-mm-ddThh:mm:ss) of observation
        d = datetime.strptime(date[:19], "%Y-%m-%dT%H:%M:%S")
        expnum = d.second + 100*(d.minute + 100*(d.hour + 100*(d.day + 100*(d.month + 100*d.year))))
        return expnum

    def get_date(self, primhdr):
        date = primhdr['DATE-OBS']
        # DATE-OBS= '2022-10-04T05:20:19.335'
        return datetime.strptime(date[:19], "%Y-%m-%dT%H:%M:%S")

    def get_mjd(self, primhdr):
        from astrometry.util.starutil_numpy import datetomjd
        d = self.get_date(primhdr)
        return datetomjd(d)

    def compute_filenames(self):
        # Rename to find masks (irafmask) and weight-maps (swarpmask)
        self.dqfn = self.imgfn.replace('.fits', '.irafmask.fits.gz')
        self.wtfn = self.imgfn.replace('.fits', '.swarpmask.fits.gz')
        assert(self.dqfn != self.imgfn)
        assert(self.wtfn != self.imgfn)

    def get_base_name(self):
        # Returns the base name to use for this Image object.  This is
        # used for calib paths, and is joined with the CCD name to
        # form the name of this Image object and for calib filenames.
        basename = os.path.basename(self.image_filename)
        ### HACK -- keep only the first dotted component of the base filename.
        # This allows, eg, create-testcase.py to use image filenames like BASE.N3.fits
        # with only a single HDU.
        basename_list = basename.split('.')[:2]
        basename = basename_list[0] + basename_list[1]
        return basename

    def set_calib_filenames(self):
        super().set_calib_filenames()
        calibdir = self.survey.get_calib_dir()
        imgdir = os.path.dirname(self.image_filename)
        basename = self.get_base_name()
        if len(self.ccdname):
            calname = basename + '-' + self.ccdname
        else:
            calname = basename
        self.scamp_fn = self.imgfn.replace('.fits', '.head')
        self.sefn         = os.path.join(calibdir, 'se',           imgdir, basename, calname + '-se.fits')
        self.psffn        = os.path.join(calibdir, 'psfex-single', imgdir, basename, calname + '-psfex.fits')
        self.skyfn        = os.path.join(calibdir, 'sky-single',   imgdir, basename, calname + '-splinesky.fits')
        self.merged_psffn = os.path.join(calibdir, 'psfex',        imgdir, basename + '-psfex.fits')
        self.merged_skyfn = os.path.join(calibdir, 'sky',          imgdir, basename + '-splinesky.fits')
        self.wcs_initial_fn = os.path.join(imgdir, 'test.wcs')

    def read_scamp_wcs(self, hdr=None):
        from astrometry.util.util import wcs_pv2sip_hdr
        import tempfile

        print('Reading Scamp file', self.scamp_fn, 'HDU', self.hdu)
        lines = open(self.scamp_fn,'rb').readlines()
        lines = [line.strip() for line in lines]
        iline = 0
        header = []
        # find my HDU in the header
        for i in range(1, self.hdu+1):
            header = []
            while True:
                if iline >= len(lines):
                    raise RuntimeError('Failed to find HDU %i in Scamp header file %s' %
                                       (self.hdu, self.scamp_fn))
                line = lines[iline]
                header.append(line)
                iline += 1
                if line == b'END':
                    break

        # print('Keeping Scamp header:')
        # for line in header:
        #     print(line)

        # Write to a temp file and then read w/ fitsio!
        tmp = tempfile.NamedTemporaryFile(delete=False)
        preamble = [b'SIMPLE  =                    T / file does conform to FITS standard',
                    b'BITPIX  =                   16 / number of bits per data pixel',
                    b'NAXIS   =                    0 / number of data axes',
                    b'EXTEND  =                    T / FITS dataset may contain extensions'
                    ]
        hdrstr = b''.join([s + b' '*(80-len(s)) for s in preamble + header])
        tmp.write(hdrstr + b' '*(2880-(len(hdrstr)%2880)))
        tmp.close()
        scamp_hdr = fitsio.read_header(tmp.name)
        del tmp
        # print('Parsed scamp header:', scamp_hdr)

        # Read original WCS header
        if hdr is None:
            hdr = self.read_image_header()

        # Copy Scamp header cards in...
        for key in ['EQUINOX', 'RADESYS', 'CTYPE1', 'CTYPE2', 'CUNIT1', 'CUNIT2',
                    'CRVAL1', 'CRVAL2', 'CRPIX1', 'CRPIX2',
                    'CD1_1', 'CD1_2', 'CD2_1', 'CD2_2',
                    'PV1_0',  'PV1_1', 'PV1_2', 'PV1_4', 'PV1_5', 'PV1_6',
                   # 'PV1_7', 'PV1_8', 'PV1_9', 'PV1_10' ,
                    'PV2_0',  'PV2_1', 'PV2_2', 'PV2_4', 'PV2_5', 'PV2_6',
                   # 'PV2_7', 'PV2_8', 'PV2_9', 'PV2_10'
                    ]:
            hdr[key] = scamp_hdr[key]
        hdr['PV2_3'] = 0
        wcs = wcs_pv2sip_hdr(hdr)
        return wcs
  
    def get_wcs(self, hdr=None):
        if self.scamp_wcs is not None:
            return self.scamp_wcs
        if not os.path.exists(self.scamp_fn):
            print(self.scamp_fn)
        # Look for Scamp "head" file
        if os.path.exists(self.scamp_fn):
            # Load Scamp WCS
            self.scamp_wcs = self.read_scamp_wcs(hdr=hdr)
        if self.scamp_wcs is not None:
            # print("Reading scamp wcs...")
            return self.scamp_wcs

    def get_fwhm(self, primhdr, imghdr):
        seeing = float(primhdr['SEEING'])
        # print("Seeing is..", seeing)
        # If PsfEx file exists, read FWHM from there
        if seeing == 0.0:
            if not hasattr(self, 'merged_psffn'):
                return super().get_fwhm(primhdr, imghdr)
            psf = self.read_psf_model(0, 0, pixPsf=True)
            fwhm = psf.fwhm
            return fwhm
        fwhm = seeing / self.pixscale
        return fwhm

    def colorterm_gaia_to_observed(self, cat, band):
        # See, eg, ps1cat.py's ps1_to_decam.
        # "cat" is a table of PS1 stars;
        # Grab the g-i color:
       # g_index = ps1cat.ps1band['g']
       # i_index = ps1cat.ps1band['i']
       # gmag = cat[:,g_index]
       # imag = cat[:,i_index]
       # gi = gmag - imag

        coeffs = dict(
            g = [ 0. ]
            )[band]
        colorterm = np.zeros(len(gi))
        for power,coeff in enumerate(coeffs):
            colorterm += coeff * 1 **power
        return colorterm
