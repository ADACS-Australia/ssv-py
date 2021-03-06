from copy import deepcopy
from enum import Enum
from astropy.io import registry
import pathlib
import os
import sys

import astropy.io.fits as fits
from astropy.nddata import (
    VarianceUncertainty,
    StdDevUncertainty,
    InverseVariance,
)
import astropy.units as u
from astropy.wcs import WCS
from astropy.wcs.utils import pixel_to_pixel

from specutils import Spectrum1D, SpectrumList
from specutils.io.registers import data_loader

_unregistered_loaders = {}
_unregistered_identifiers = {}
_low_priority_loaders = {}

def whatformat(*args, format=None, **kwargs):
    """
    Read in data.

    The arguments passed to this method depend on the format.
    """
    cls = SpectrumList
    formats = []
    ctx = None
    try:
        if format is None:
            path = None
            fileobj = None

            if len(args):
                if isinstance(args[0], (str, pathlib.Path)) and not os.path.isdir(args[0]):
                    from astropy.utils.data import get_readable_fileobj
                    # path might be a pathlib.Path object
                    if isinstance(args[0], pathlib.Path):
                        args = (str(args[0]),) + args[1:]
                    path = args[0]
                    try:
                        ctx = get_readable_fileobj(args[0], encoding='binary')
                        fileobj = ctx.__enter__()
                    except OSError:
                        raise
                    except Exception:
                        fileobj = None
                    else:
                        args = [fileobj] + list(args[1:])
                elif hasattr(args[0], 'read'):
                    path = None
                    fileobj = args[0]

                formats = whatformat_get_valid_format(
                    'read', cls, path, fileobj, args, kwargs)
        else:
            formats = [format]

    finally:
        if ctx is not None:
            ctx.__exit__(*sys.exc_info())

    return formats

def unregister(format):
    """
    The unregistered format also stays the less preferred format even once its restored
    """
    _low_priority_loaders[format] = True
    if format not in _unregistered_loaders.keys():
        _unregistered_loaders[format] = SpectrumList
        _unregistered_identifiers[format] = registry._identifiers[(format, SpectrumList)]
        registry.unregister_identifier(format, SpectrumList)

def restore_registered_loaders():
    """
    Restore all unregistered formats
    """
    for format in list(_unregistered_loaders):
        registry.register_identifier(format, SpectrumList, _unregistered_identifiers[format])
        _unregistered_loaders.pop(format)
        _unregistered_identifiers.pop(format)
    pass

def whatformat_get_valid_format(mode, cls, path, fileobj, args, kwargs):
    """
    Returns the first valid format that can be used to read/write the data in
    question.  Mode can be either 'read' or 'write'.
    """

    valid_formats = registry.identify_format(mode, cls, path, fileobj, args, kwargs)

    # A low priority loader is only intended to be used if its the only one
    reordered = []
    for format in valid_formats:
        if format in _low_priority_loaders:
            reordered.append(format)

    for format in valid_formats:
        if not format in _low_priority_loaders:
            reordered.append(format)
    return reordered

HEADER_PUPOSE_KEYWORDS = ["EXTNAME", "HDUNAME"]
HEADER_INDEX_PUPOSE_KEYWORDS = ["ROW", "ARRAY"]
FITS_FILE_EXTS = ["fit", "fits", "fts"]
SINGLE_SPLIT_LABEL = "Data Central Single-Split"
MULTILINE_SINGLE_LABEL = "Data Central Multiline-Single"
UNKNOWN_LABEL = "Unable to find a sensible label for spectrum"
# These are order in a best guess of priority, ideally the loader would know
# which label to use.
HEADER_LABEL_KEYWORDS = [
    "OBJECT",
    "OBJNAME",
    "OBS_ID",
    "EXTNAME",
    "HDUNAME",
    "TITLE",
    "ORIGIN",
    "ROOTNAME",
    "FILENAME",
    "AUTHOR",
    "OBSERVER",
    "CREATOR",
    "INSTRUME",
    "PROGRAM",
]


def guess_label_from_header(header):
    """
    Guess the label from `header`, which is assumed to be some mapping with
    FITS-like keys.
    """
    for header_key in HEADER_LABEL_KEYWORDS:
        label = header.get(header_key)
        if label is not None:
            return str(label)
    raise ValueError(UNKNOWN_LABEL)


class Purpose(Enum):
    SKIP = "skip"
    SCIENCE = "science"
    ERROR_STDEV = "error_stdev"
    ERROR_VARIANCE = "error_variance"
    ERROR_INVERSEVARIANCE = "error_inversevariance"
    SKY = "sky"
    COMBINED_SCIENCE = "combined_science"
    COMBINED_ERROR_STDEV = "combined_error_stdev"
    COMBINED_ERROR_VARIANCE = "combined_error_variance"
    COMBINED_ERROR_INVERSEVARIANCE = "combined_error_inversevariance"
    UNREDUCED_SCIENCE = "unreduced_science"
    UNREDUCED_ERROR_STDEV = "unreduced_error_stdev"
    UNREDUCED_ERROR_VARIANCE = "unreduced_error_variance"
    UNREDUCED_ERROR_INVERSEVARIANCE = "unreduced_error_inversevariance"


CREATE_SPECTRA = {
    Purpose.SCIENCE,
    Purpose.SKY,
    Purpose.COMBINED_SCIENCE,
    Purpose.UNREDUCED_SCIENCE,
}
ERROR_PURPOSES = {
    Purpose.ERROR_STDEV,
    Purpose.ERROR_VARIANCE,
    Purpose.ERROR_INVERSEVARIANCE,
    Purpose.COMBINED_ERROR_STDEV,
    Purpose.COMBINED_ERROR_VARIANCE,
    Purpose.COMBINED_ERROR_INVERSEVARIANCE,
    Purpose.UNREDUCED_ERROR_STDEV,
    Purpose.UNREDUCED_ERROR_VARIANCE,
    Purpose.UNREDUCED_ERROR_INVERSEVARIANCE,
}
PURPOSE_SPECTRA_MAP = {
    Purpose.SCIENCE: "reduced",
    Purpose.ERROR_STDEV: "reduced",
    Purpose.ERROR_VARIANCE: "reduced",
    Purpose.ERROR_INVERSEVARIANCE: "reduced",
    Purpose.SKY: "sky",
    Purpose.COMBINED_SCIENCE: "combined",
    Purpose.COMBINED_ERROR_STDEV: "combined",
    Purpose.COMBINED_ERROR_VARIANCE: "combined",
    Purpose.COMBINED_ERROR_INVERSEVARIANCE: "combined",
    Purpose.UNREDUCED_SCIENCE: "unreduced",
    Purpose.UNREDUCED_ERROR_STDEV: "unreduced",
    Purpose.UNREDUCED_ERROR_VARIANCE: "unreduced",
    Purpose.UNREDUCED_ERROR_INVERSEVARIANCE: "unreduced",
}
UNCERTAINTY_MAP = {
    Purpose.ERROR_STDEV: StdDevUncertainty,
    Purpose.ERROR_VARIANCE: VarianceUncertainty,
    Purpose.ERROR_INVERSEVARIANCE: InverseVariance,
    Purpose.COMBINED_ERROR_STDEV: StdDevUncertainty,
    Purpose.COMBINED_ERROR_VARIANCE: VarianceUncertainty,
    Purpose.COMBINED_ERROR_INVERSEVARIANCE: InverseVariance,
    Purpose.UNREDUCED_ERROR_STDEV: StdDevUncertainty,
    Purpose.UNREDUCED_ERROR_VARIANCE: VarianceUncertainty,
    Purpose.UNREDUCED_ERROR_INVERSEVARIANCE: InverseVariance,
}
GUESS_TO_PURPOSE = {
    "badpix": Purpose.SKIP,
    "": Purpose.SKIP,
    "sky": Purpose.SKY,
    "stdev": Purpose.ERROR_STDEV,
    "sigma": Purpose.ERROR_STDEV,
    "variance": Purpose.ERROR_VARIANCE,
    "spectrum": Purpose.SCIENCE,
}


def add_labels(spec_list, use_purpose=False):
    not_labeled = 0
    label_set = set()
    for spec in spec_list:
        meta = spec.meta
        purpose = meta.get("purpose")
        if use_purpose:
            tail = " (" + str(purpose) + ")"
        else:
            tail = ""
        try:
            meta["label"] = guess_label_from_header(meta["header"]) + tail
        except ValueError:
            not_labeled += 1
        else:
            label_set.add(meta["label"])

    if len(label_set) + not_labeled < len(spec_list):
        # This implies there are duplicates
        for i, spec in enumerate(spec_list, start=1):
            label = spec.meta.get("label")
            if label is not None:
                spec.meta["label"] = label + " #" + str(i)


def compute_wcs_from_keys_and_values(
    header=None,
    *,
    wavelength_unit_keyword=None,
    wavelength_unit=None,
    pixel_reference_point_keyword=None,
    pixel_reference_point=None,
    pixel_reference_point_value_keyword=None,
    pixel_reference_point_value=None,
    pixel_width_keyword=None,
    pixel_width=None,
):
    if wavelength_unit is None:
        if wavelength_unit_keyword is None:
            raise ValueError(
                "Either wavelength_unit or wavelength_unit_keyword must be "
                "provided"
            )
        wavelength_unit = u.Unit(header[wavelength_unit_keyword])
    if pixel_reference_point is None:
        if pixel_reference_point_keyword is None:
            raise ValueError(
                "Either pixel_reference_point or "
                "pixel_reference_point_keyword must be provided"
            )
        pixel_reference_point = header[pixel_reference_point_keyword]
    if pixel_reference_point_value is None:
        if pixel_reference_point_value_keyword is None:
            raise ValueError(
                "Either pixel_reference_point_value or "
                "pixel_reference_point_value_keyword must be provided"
            )
        pixel_reference_point_value = header[
            pixel_reference_point_value_keyword
        ]
    if pixel_width is None:
        if pixel_width_keyword is None:
            raise ValueError(
                "Either pixel_width or pixel_width_keyword must be provided"
            )
        # RS: The javascript drivers allow for the keyword CDELT1 (TODO: discuss with James)
        if header.get(pixel_width_keyword) is None:
            if header.get("CDELT1") is not None:
                pixel_width_keyword = "CDELT1"
        pixel_width = header[pixel_width_keyword]

    w = WCS(naxis=1)
    w.wcs.crpix[0] = pixel_reference_point
    w.wcs.crval[0] = pixel_reference_point_value
    w.wcs.cdelt[0] = pixel_width
    w.wcs.cunit[0] = wavelength_unit
    return w


def get_flux_units_from_keys_and_values(
    header,
    *,
    flux_unit_keyword=None,
    flux_unit=None,
    flux_scale_keyword=None,
    flux_scale=None,
):
    if flux_unit is None:
        if flux_unit_keyword is None:
            raise ValueError(
                "Either flux_unit or flux_unit_keyword must be provided"
            )
        flux_unit = header[flux_unit_keyword]
    flux_unit = u.Unit(flux_unit)

    if flux_scale is None:
        if flux_scale_keyword is None:
            flux_scale = 1
        else:
            flux_scale = header[flux_scale_keyword]
    return flux_scale * flux_unit


def add_single_spectra_to_map(
    spectra_map,
    *,
    header,
    data,
    spec_info=None,
    wcs_info=None,
    units_info=None,
    purpose_prefix=None,
    all_standard_units,
    all_keywords,
    valid_wcs,
    index=None,
):
    spec_wcs_info = {}
    spec_units_info = {}
    if wcs_info is not None:
        spec_wcs_info.update(wcs_info)
    if units_info is not None:
        spec_units_info.update(units_info)

    if spec_info is not None:
        spec_wcs_info.update(spec_info.get("wcs", {}))
        spec_units_info.update(spec_info.get("units", {}))
        purpose = spec_info.get("purpose")
    else:
        purpose = None

    purpose = get_purpose(
        header,
        purpose=purpose,
        purpose_prefix=purpose_prefix,
        all_keywords=all_keywords,
        index=index,
    )

    if purpose == Purpose.SKIP:
        return None

    if valid_wcs or not spec_wcs_info:
        wcs = WCS(header)
    else:
        wcs = compute_wcs_from_keys_and_values(header, **spec_wcs_info)

    if all_standard_units:
        spec_units_info = {"flux_unit_keyword": "BUNIT"}
    flux_unit = get_flux_units_from_keys_and_values(header, **spec_units_info)
    flux = data * flux_unit

    meta = {"header": header, "purpose": PURPOSE_SPECTRA_MAP[purpose]}

    if purpose in CREATE_SPECTRA:
        spectrum = Spectrum1D(wcs=wcs, flux=flux, meta=meta)
        spectra_map[PURPOSE_SPECTRA_MAP[purpose]].append(spectrum)
    elif purpose in ERROR_PURPOSES:
        try:
            spectrum = spectra_map[PURPOSE_SPECTRA_MAP[purpose]][-1]
        except IndexError:
            raise ValueError(f"No spectra to associate with {purpose}")
        aligned_flux = pixel_to_pixel(wcs, spectrum.wcs, flux)
        spectrum.uncertainty = UNCERTAINTY_MAP[purpose](aligned_flux)
        spectrum.meta["uncertainty_header"] = header

    # We never actually want to return something, this just flags it to pylint
    # that we know we're breaking out of the function when skip is selected
    return None


def get_purpose(
    header, *, purpose=None, purpose_prefix=None, all_keywords, index=None
):
    def guess_purpose(header):
        for keyword in HEADER_PUPOSE_KEYWORDS:
            guess = header.get(keyword)
            if guess is not None:
                return GUESS_TO_PURPOSE[guess.strip().lower()]
        return None

    def guess_index_purpose(header, index):
        for keyword in HEADER_INDEX_PUPOSE_KEYWORDS:
            guess = header.get(keyword + str(index))
            if guess is not None:
                return GUESS_TO_PURPOSE[guess.strip().lower()]
        return None

    if all_keywords:
        if index is None:
            guessed_purpose = guess_purpose(header)
            if guessed_purpose is not None:
                return guessed_purpose
            if "XTENSION" not in header:
                # we have a primary HDU, assume science
                return Purpose.SCIENCE
            raise ValueError(
                "Cannot identify purpose, cannot use all_keywords"
            )
        guessed_purpose = guess_index_purpose(header, index)
        if guessed_purpose is not None:
            return guessed_purpose
        raise ValueError("Cannot identify purpose, cannot use all_keywords")
    if purpose is not None:
        return Purpose(purpose)
    if purpose_prefix is not None:
        if index is None:
            return Purpose(header.get(purpose_prefix))
        return Purpose(header.get(purpose_prefix + str(index)))
    raise ValueError(
        "Either all_keywords must be True, or one of purpose or "
        "purpose_prefix must not be None."
    )


def no_auto_identify(*args, **kwargs):
    return False
