import astropy.io.fits as fits
from specutils import SpectrumList
from specutils.io.registers import data_loader
import ssv
import json

from .ssvloaders import FITS_FILE_EXTS, SINGLE_SPLIT_LABEL, MULTILINE_SINGLE_LABEL

JSON_FILE_EXTS = ["json", "JSON"]


def identify_marzjson(origin, *args, **kwargs):
    """
    Identify if the current file is a OzDES file
    """
    file_obj = args[0]
    if file_obj.endswith(".json"):
        return True
    
    return False


@data_loader(
    label="MARZJSON", extensions=JSON_FILE_EXTS, dtype=SpectrumList,
    identifier=identify_marzjson,
)
def marzjson_loader(fname):
    spectra = None
    with open(fname) as f:
        spec_data = json.load(f)
        spectra = ssv.utils.fromMarzJSON(spec_data)
    return spectra