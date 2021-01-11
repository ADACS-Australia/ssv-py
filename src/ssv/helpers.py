from warnings import warn

from astropy.convolution import Box1DKernel, Gaussian1DKernel
from astropy.table import QTable
from astropy.nddata import (
    VarianceUncertainty, StdDevUncertainty, InverseVariance,
)
from astropy.units import Quantity
from numpy import isfinite, power, zeros, empty, nan, max as np_max
from specutils import Spectrum1D, SpectrumList, SpectrumCollection
from specutils.manipulation import convolution_smooth

from .plotting import (
    DEFAULT_WAVELENGTH_COLUMN_NAME, DEFAULT_FLUX_COLUMN_NAME, DEFAULT_UNCERTAINTY_COLUMN_NAME
)

UNKNOWN_LABEL = "Unable to find a sensible label for spectrum"
# These are order in a best guess of priority, ideally the loader would know
# which label to use.
HEADER_LABEL_KEYWORDS = [
    "OBJECT",
    "OBJNAME",
    "OBS_ID",
    "SOURCE",
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
DEFAULT_BOX_FILTER_WIDTH = 3
NAMED_FILTERS = {
    "box": Box1DKernel,
    "gaussian": Gaussian1DKernel,
}


def guess_label_from_header(header):
    """
    Guess the label from `header`, which is assumed to be some mapping with
    FITS-like keys.
    """
    for header_key in HEADER_LABEL_KEYWORDS:
        label = header.get(header_key)
        if label is not None:
            return label
    raise ValueError(UNKNOWN_LABEL)


def specutils_spectrum_to_table_spectrum(spectrum):
    """
    Convert a specutils spectrum to an astropy QTable ready for plotting.
    """
    if isinstance(spectrum, Spectrum1D):
        return QTable({
            DEFAULT_WAVELENGTH_COLUMN_NAME: spectrum.wavelength,
            DEFAULT_FLUX_COLUMN_NAME: spectrum.flux,
        })
    if isinstance(spectrum, SpectrumList):
        return QTable({
            DEFAULT_WAVELENGTH_COLUMN_NAME: spectrum[0].wavelength,
            DEFAULT_FLUX_COLUMN_NAME: spectrum[0].flux,
        })
    if isinstance(spectrum, SpectrumCollection):
        return QTable({
            DEFAULT_WAVELENGTH_COLUMN_NAME: spectrum.wavelength,
            DEFAULT_FLUX_COLUMN_NAME: spectrum.flux,
        })
    raise TypeError("Must be a specutils spectrum object")


def _filter_spectra(spectra, prefer_combined):
    for spectrum in spectra:
        if prefer_combined and spectrum.meta.get("purpose") == "combined":
            yield spectrum
            break
        if prefer_combined and spectrum.meta.get("purpose") in [
            "sky", "unreduced"
        ]:
            continue
        yield spectrum


def _get_spec_label(spectrum, prefix):
    label = spectrum.meta.get("label")
    if label is not None:
        label = str(label)
        if prefix is not None:
            label = str(prefix) + ": " + label
        return label
    header = spectrum.meta.get("header")
    if header is None:
        if prefix is not None:
            warn(UNKNOWN_LABEL)
            return str(prefix)
        raise ValueError(UNKNOWN_LABEL)
    try:
        label = guess_label_from_header(header)
    except ValueError as e:
        if prefix is not None:
            warn(str(e))
            return str(prefix)
        raise e
    label = str(label)
    if prefix is not None:
        label = str(prefix) + ": " + label
    return label


def specutils_spectra_to_table_spectra(
    spec_list, *, prefix=None, prefer_combined=False
):
    """
    Convert a list of specutils spectrum to an astropy QTable ready for
    plotting.
    """
    if isinstance(spec_list, (Spectrum1D, SpectrumCollection)):
        spec_list = [spec_list]

    table_spectra = {}
    for spectrum in _filter_spectra(spec_list, prefer_combined):
        table_dict = {
            DEFAULT_WAVELENGTH_COLUMN_NAME: spectrum.wavelength,
            DEFAULT_FLUX_COLUMN_NAME: spectrum.flux
        }

        if isinstance(spectrum.uncertainty, StdDevUncertainty):
            table_dict[DEFAULT_UNCERTAINTY_COLUMN_NAME] = power(spectrum.uncertainty.quantity, 2)
        elif isinstance(spectrum.uncertainty, InverseVariance):
            table_dict[DEFAULT_UNCERTAINTY_COLUMN_NAME] = power(spectrum.uncertainty.quantity, -1)
        elif isinstance(spectrum.uncertainty, VarianceUncertainty):
            table_dict[DEFAULT_UNCERTAINTY_COLUMN_NAME] = spectrum.uncertainty.quantity
        else:
            blank_uncertainty = empty(len(spectrum.spectral_axis))
            blank_uncertainty.fill(nan)
            table_dict[DEFAULT_UNCERTAINTY_COLUMN_NAME] = Quantity(blank_uncertainty, unit=spectrum.flux.unit**2)

        # label = spectrum.meta.get('purpose') + ' : ' + _get_spec_label(spectrum, prefix)
        label = spectrum.meta.get('purpose')
        table_spectra[label] = QTable(table_dict)

    return table_spectra


def smooth_spectra(spectra, *, filter="box", **kwargs):
    """
    Smooth spectra via a given filter.
    """
    if filter == "box" and not kwargs:
        kwargs["width"] = DEFAULT_BOX_FILTER_WIDTH

    if filter in NAMED_FILTERS:
        filter = NAMED_FILTERS[filter]

    if isinstance(spectra, Spectrum1D):
        return convolution_smooth(spectra, filter(**kwargs))
    if isinstance(spectra, SpectrumList):
        return SpectrumList([
            convolution_smooth(spectrum, filter(**kwargs))
            for spectrum in spectra
        ])
    if isinstance(spectra, SpectrumCollection):
        return SpectrumCollection([
            convolution_smooth(spectrum, filter(**kwargs))
            for spectrum in spectra
        ])

    try:
        spectra_iter = iter(spectra)
    except TypeError:
        raise TypeError(
            "Must be a specutils spectrum object or an iterable object"
        )
    return [
        convolution_smooth(spectrum, filter(**kwargs))
        for spectrum in spectra_iter
    ]
