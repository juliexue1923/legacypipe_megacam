import os
from datetime import datetime
import fitsio
import numpy as np

from legacypipe.image import LegacySurveyImage

class MegaCamImage(LegacySurveyImage):

    zp0 = dict(
        u = 25.0,
        g = 25.0,
        r = 25.0
    )
      
    k_ext = dict(# copied from decam
                 g = 0.173,
                 r = 0.090,
                 # From Arjun 2021-03-17 based on DECosmos (calib against SDSS)
                 u = 0.63
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
    
    def get_nominal_pixscale():
        # from McLeod 2015
        return 0.16

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

    def get_band(self, primhdr):
        band = primhdr['FILTER']
        band = band.split()[0]
        if band == 'open,g':
            return 'g'
        if band == 'open,r':
            return 'r'
        if band == 'open,u':
            return 'u'
        return band
    
    def set_ccdzpt(self, ccdzpt):
        # Adjust zeropoint for exposure time
        self.ccdzpt = ccdzpt + 2.5 * np.log10(self.exptime)

    def compute_filenames(self):
        # Rename to find masks (irafmask) and weight-maps (swarpmask)
        self.dqfn = self.imgfn.replace('.fits', '.irafmask.fits.gz')
        self.wtfn = self.imgfn.replace('.fits', '.swarpmask.fits.gz')
        self.dqfn = self.dqfn.replace('working', 'working/../data_quality_dir')
        self.wtfn = self.wtfn.replace('working', 'working/../data_quality_dir')
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
            if key in scamp_hdr.keys():
                hdr[key] = scamp_hdr[key]
            else:
                print("Warning: no " + str(key) + " key. Set to 0.")
                hdr[key] = 0.
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
        if 'SEEING' in primhdr: 
            seeing = float(primhdr['SEEING'])
        # If PsfEx file exists, read FWHM from there
        else:
            if not hasattr(self, 'merged_psffn'):
                return super().get_fwhm(primhdr, imghdr)
            psf = self.read_psf_model(0, 0, pixPsf=True)
            fwhm = psf.fwhm
            return fwhm
        import numpy as np
        if seeing == 0.0 or seeing < 0.0 or np.isnan(seeing):
            if not hasattr(self, 'merged_psffn'):
                return super().get_fwhm(primhdr, imghdr)
            psf = self.read_psf_model(0, 0, pixPsf=True)
            fwhm = psf.fwhm 
            return fwhm
        fwhm = seeing / self.pixscale
        return fwhm

    def get_photometric_calibrator_cuts(self, name, cat):
        '''Returns whether to keep sources in the *cat* of photometric calibration
        stars from, eg, Pan-STARRS1 or SDSS.
        '''
        if name == 'delve':
            return np.ones(len(cat), bool)
        if name == 'smss':
            return np.ones(len(cat), bool)
        if name == 'sdss':
            return np.ones(len(cat), bool)
        if name == 'gxp':
            return np.ones(len(cat), bool)
        raise RuntimeError('Unknown photometric calibration set: %s' % name)

    def photometric_calibrator_to_observed(self, name, cat):
        if name == 'delve':
            band = self.get_delve_band()
            return getattr(cat, band)
        elif name == 'smss':
            band = self.get_smss_band()
        #    return getattr(cat, band)
            colorterm = self.colorterm_smss_to_observed(cat, self.band)
            return getattr(cat, band) + np.clip(colorterm, -1., +1.)
        elif name == 'sdss':
            band = self.get_sdss_band()
            return getattr(cat, band)
        elif name == 'gxp':
            band = self.get_gxp_band()
            colorterm = self.colorterm_gxp_to_observed(cat, self.band)
            return getattr(cat, band) + colorterm
        else:
            raise RuntimeError('No photometric conversion from %s to camera' % name)
    
    def get_smss_band(self):
        from legacypipe.smsscat import SMSSCatalog
        # A known filter?
        if self.band in SMSSCatalog.smssband:
            return SMSSCatalog.smssband[self.band]

    def get_sdss_band(self):
        from legacypipe.ps1cat import sdsscat
        # A known filter?
        if self.band in sdsscat.sdssband:
            return sdsscat.sdssband[self.band]

    def get_delve_band(self):
        from legacypipe.delvecat import DELVECatalog
        # A known filter?
        if self.band in DELVECatalog.delveband:
            return DELVECatalog.delveband[self.band]
        else:
            raise RuntimeError('This band is not in DELVE. Choose another photometric calibrator catalog.') 

    def get_gxp_band(self):
        from legacypipe.gaiaxp_cat import GaiaXPCatalog
        # A known filter?
        if self.band in GaiaXPCatalog.gxpband:
            return GaiaXPCatalog.gxpband[self.band]

    def colorterm_gxp_to_observed(self, cat, band):
        from legacypipe.gaiaxp_cat import GaiaXPCatalog
        # See, eg, ps1cat.py's ps1_to_decam.
        # "cat" is a table of reference stars;
        # Grab the g-r and u-r color:
        gmag = getattr(cat, 'synth_gmag')
        rmag = getattr(cat, 'synth_rmag')
        umag = getattr(cat, 'synth_umag')
        gr = gmag - rmag
        ur = umag - rmag

        # np.polyfit returns coefficients counting from the highest order, but enumerate() pairs the first entry with 0.
        # As such, reverse the order from np.polyfit
        u_coeffs = [
                (-1.0,  1.85, [-0.03997513, 0.0448199, 0.02125202]),
               # ( 0.5,  2.5, [ -0.0687748, 0.09175484]),
               # ( 2.5, 4.15, [ 1.93843559, -1.35614424, 0.2572292])
                ]
        
        if band in ['g','r']:
            colorterm = np.zeros(len(gr))
            return colorterm
        if band in ['u']:
            colorterm = np.zeros(len(ur))
            for ur_min, ur_max, coeffs in u_coeffs:
                mask = (ur >= ur_min) & (ur < ur_max)
                for power, coeff in enumerate(coeffs):
                    colorterm[mask] += coeff * ur[mask]**power
                    print(colorterm)
                    return colorterm

    def colorterm_smss_to_observed(self, cat, band):
        from legacypipe.smsscat import SMSSCatalog
        # See, eg, ps1cat.py's ps1_to_decam.
        # "cat" is a table of reference stars;
        # Grab the g-r, u-g and u-r color:
        g_index = 7
        r_index = 11
        u_index = 3
        gmag = getattr(cat, 'gmag')
        rmag = getattr(cat, 'rmag')
        umag = getattr(cat, 'umag')
        gr = gmag - rmag
        ug = umag - gmag
        ur = umag - rmag
        print(gmag, umag)

        # np.polyfit returns coefficients counting from the highest order, but enumerate() pairs the first entry with 0.
        # As such, reverse the order from np.polyfit
        coeffs = dict(
            g = [-0.26456767,  0.60549870 ],
            r = [ 0.06496791, -0.09319692 ],
            u = [ 0.02854417, -0.02239842 ]
            )[band]

        # Most are with respect to g-r, some are u-r...
        color = gr
        if band in ['u']:
            color = ur

        colorterm = np.zeros(len(color))
        for power,coeff in enumerate(coeffs):
            colorterm += coeff * color**power
        return colorterm

    def check_image_header(self, imghdr):
        # check consistency between the CCDs table and the image header 
        if 'EXTNAME' in imghdr:
            e =  imghdr['EXTNAME'].strip().upper()
        elif 'MEXTNO' in imghdr:
            ext_num = imghdr['MEXTNO'].strip().upper()
            e = 'IM' + ext_num
        if e.strip() != self.ccdname.strip():
            warnings.warn('Expected header EXTNAME="%s" to match self.ccdname="%s", self.imgfn=%s' % (e.strip(), self.ccdname,self.imgfn))
