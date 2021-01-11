import astropy.io.fits as fits
from specutils import SpectrumList
from specutils.io.registers import data_loader

from .loaders import FITS_FILE_EXTS, SINGLE_SPLIT_LABEL

OZDES_CONFIG = {
    "hdus": {
        "0": {"purpose": "combined_science"},
        "1": {"purpose": "combined_error_variance"},
        "2": {"purpose": "skip"},
        "cycle": {
            "0": {"purpose": "science"},
            "1": {"purpose": "error_variance"},
            "2": {"purpose": "skip"},
        },
    },
    "units": None,
    "wcs": None,
    "all_standard_units": True,
    "all_keywords": False,
    "valid_wcs": True,
}


def identify_ozdes(origin, *args, **kwargs):
    """
    Identify if the current file is a OzDES file
    """
    file_obj = args[0]
    if isinstance(file_obj, fits.hdu.hdulist.HDUList):
        hdulist = file_obj
    else:
        hdulist = fits.open(file_obj, **kwargs)

    if "ozdes" in hdulist[0].header.get("REFERENC"):
        if not isinstance(file_obj, fits.hdu.hdulist.HDUList):
            hdulist.close()
        return True

    if not isinstance(file_obj, fits.hdu.hdulist.HDUList):
        hdulist.close()
    return False


@data_loader(
    label="OzDES", extensions=FITS_FILE_EXTS, dtype=SpectrumList,
    identifier=identify_ozdes,
)
def ozdes_loader(fname):
    spectra = SpectrumList.read(
        fname, format=SINGLE_SPLIT_LABEL, **OZDES_CONFIG
    )
    return spectra
