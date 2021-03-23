import astropy.io.fits as fits
from specutils import SpectrumList
from specutils.io.registers import data_loader

from .ssvloaders import FITS_FILE_EXTS, SINGLE_SPLIT_LABEL, MULTILINE_SINGLE_LABEL

MARZ_CONFIG = {
    "hdus": {
        "0": {
            "purpose": "science",
            "units": {"flux_unit": "10^-17 erg/s/cm^2/A"},
        },
        "1": {"purpose": "error_variance"},
        "2": {"purpose": "sky"},
        "3": {"purpose": "skip"},
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
# Try this one first (needs a patch to specutils until they except the change)
MARZ_CONFIG_NEXT = {
    "hdus": {
        "0": {
            "purpose": "science",
            "units": {"flux_unit": "10^-17 erg/s/cm^2/A"},
        },
        "1": {"purpose": "error_variance"},
        "2": {"purpose": "sky"},
        "3": {"purpose": "skip"},
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
    "fallback_header": True,
}


def identify_marz(origin, *args, **kwargs):
    """
    Identify if the current file is a OzDES file
    """
    file_obj = args[0]
    if isinstance(file_obj, fits.hdu.hdulist.HDUList):
        hdulist = file_obj
    else:
        hdulist = fits.open(file_obj, **kwargs)

    header = hdulist[0].header

    if "AAOMEGA-2dF" in header.get("INSTRUME", "") and header.get("NAXIS", 0) == 1:
        if not isinstance(file_obj, fits.hdu.hdulist.HDUList):
            hdulist.close()
        return True

    if "Combined" in header.get("SOURCE", ""):
        if not isinstance(file_obj, fits.hdu.hdulist.HDUList):
            hdulist.close()
        return True

    if not isinstance(file_obj, fits.hdu.hdulist.HDUList):
        hdulist.close()
    return False


@data_loader(
    label="MARZ", extensions=FITS_FILE_EXTS, dtype=SpectrumList,
    identifier=identify_marz,
)
def marz_loader(fname):
    try:
        spectra = SpectrumList.read(
            fname, format=SINGLE_SPLIT_LABEL, **MARZ_CONFIG_NEXT
        )
    except:
        # fallback (to not having support for "fallback_header"!)
        spectra = SpectrumList.read(
            fname, format=SINGLE_SPLIT_LABEL, **MARZ_CONFIG
        )
    return spectra