import astropy.io.fits as fits
from specutils import SpectrumList
from specutils.io.registers import data_loader

from .loaders import FITS_FILE_EXTS, SINGLE_SPLIT_LABEL

GALAH_CONFIG = {
    "hdus": {
        "0": {"purpose": "science"},
        "1": {"purpose": "error_stdev"},
        "2": {"purpose": "unreduced_science"},
        "3": {"purpose": "unreduced_error_stdev"},
        "4": {"purpose": "skip"},
    },
    "wcs": {
        "pixel_reference_point_keyword": "CRPIX1",
        "pixel_reference_point_value_keyword": "CRVAL1",
        "pixel_width_keyword": "CDELT1",
        "wavelength_unit": "Angstrom",
    },
    "units": {"flux_unit": "count"},
    "all_standard_units": False,
    "all_keywords": False,
    "valid_wcs": False,
}


def identify_galah(origin, *args, **kwargs):
    """
    Identify if the current file is a GALAH file
    """
    file_obj = args[0]
    if isinstance(file_obj, fits.hdu.hdulist.HDUList):
        hdulist = file_obj
    else:
        hdulist = fits.open(file_obj, **kwargs)

    if "galah" in hdulist[0].header.get("REFERENC"):
        if not isinstance(file_obj, fits.hdu.hdulist.HDUList):
            hdulist.close()
        return True

    if not isinstance(file_obj, fits.hdu.hdulist.HDUList):
        hdulist.close()
    return False


@data_loader(
    label="GALAH", extensions=FITS_FILE_EXTS, dtype=SpectrumList,
    identifier=identify_galah,
)
def galah_loader(fname):
    spectra = SpectrumList.read(
        fname, format=SINGLE_SPLIT_LABEL, **GALAH_CONFIG
    )
    return spectra
