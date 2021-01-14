import numpy as np
from functools import reduce
from specutils import SpectrumList, Spectrum1D
from specutils.fitting import fit_generic_continuum, fit_continuum
from specutils.manipulation import snr_threshold, box_smooth, gaussian_smooth, trapezoid_smooth, convolution_smooth, median_smooth
import astropy.units as u
from astropy.io import fits, registry
from pathlib import Path
import json, os

# From SO post https://stackoverflow.com/a/39130019
# Allows applying multiple functions in sequence to data
def id(x):
    return x

def comp(f,g):
    return lambda x: f(g(x))

def compose(*fs):
    return reduce(comp, fs, id)

# Transform functions
def fit_continuum(spectrum):
    """Fits continuum to a spectrum

    Parameters
    ----------
    spectrum : specutils.Spectrum1D
        Spectrum to which the continuum is to be fitted

    Returns
    -------
    numpy.NDArray
        Array of values for the fitted continuum flux
    """

    fitted_continuum = fit_generic_continuum(spectrum)
    return fitted_continuum(spectrum.spectral_axis)

def subtract_continuum(spectrum):
    """Fit continuum and subtract from spectrum

    Parameters
    ----------
    spectrum : specutils.Spectrum1D
        Spectrum from which the continuum is to be subtracted

    Returns
    -------
    specutils.Spectrum1D
        Spectrum with the continuum subtracted
    """

    spectrum = spectrum - fit_continuum(spectrum)
    return spectrum

def remove_spurious_points(spectrum):
    """Remove points from a spectrum that have a signal-to-noise ratio of above 50

    Parameters
    ----------
    spectrum : specutils.Spectrum1D
        Spectrum from which the spurious points are to be removed

    Returns
    -------
    specutils.Spectrum1D
        Spectrum with the spurious points removed
    """

    #uncertainty should be stddev for this function
    if spectrum.uncertainty is None or all(np.isnan(spectrum.uncertainty.array)):
        return spectrum
    else:
        spectrum_masked = snr_threshold(spectrum, 50)
        flux_masked = spectrum_masked.flux
        flux_masked[np.logical_not(spectrum_masked.mask)] = np.nan
        return Spectrum1D(spectral_axis=spectrum_masked.spectral_axis, flux=flux_masked, uncertainty=spectrum_masked.uncertainty, mask=spectrum_masked.mask)

def apply_scaling(max_value=1):
    """Returns a function used to scale the flux values of a spectrum to a maximum value

    Parameters
    ----------
    max_value : int, optional
        Maximum value of the scaled spectrum, by default 1 (normalisation)

    Returns
    -------
    function
        Takes a spectrum and returns spectrum scaled to the `max_value`
    """    

    def scale_flux(spectrum):
        min_flux = np.nanmin(spectrum.flux)
        max_flux = np.nanmax(spectrum.flux)
        spectrum = ((spectrum - min_flux) / (max_flux - min_flux)) * max_value
        return spectrum
    return scale_flux

def apply_smoothing(smoothing_func=median_smooth, smoothing_width=5):
    """Return a function used to smooth an input spectrum

    Parameters
    ----------
    smoothing_func : function, optional
        Function which takes a Spectrum1D and a smoothing width, by default median_smooth
    smoothing_width : int, optional
        Width of the smoothing kernel, by default 5

    Returns
    -------
    function
        Takes a spectrum and returns a spectrum smoothed with the specified smoothing function
    """    
    def smooth_flux(spectrum):
        return smoothing_func(spectrum, smoothing_width)
    return smooth_flux

def offset_wavelength(offset, spectrum):
    """Offset the wavelength of a spectrum by a given amount

    Parameters
    ----------
    offset : int, float or astropy.units.Quantity
        The amount by which to offset the wavelength. If an `int` or `float`, then will use the current units of the wavelength.
    spectrum : specutils.Spectrum1D
        Spectrum to offset

    Returns
    -------
    specutils.Spectrum1D
        Spectrum with wavelength offset by desired amount
    """
    if not isinstance(offset, u.Quantity):
        offset *= spectrum.spectral_axis.unit
    return Spectrum1D(spectral_axis=spectrum.spectral_axis + offset, flux=spectrum.flux, uncertainty=spectrum.uncertainty, mask=spectrum.mask)

def redshift_wavelength(redshift_scale_factor, spectrum):
    """Redshift the spectrum by a given scale factor

    Parameters
    ----------
    redshift_scale_factor : int or float
        Multiplies the wavelength values of the spectrum
    spectrum : specutils.Spectrum1D
        Spectrum to offset

    Returns
    -------
    specutils.Spectrum1D
        Spectrum with wavelength offset by desired amount
    """
    return Spectrum1D(spectral_axis=spectrum.spectral_axis * redshift_scale_factor, flux=spectrum.flux, uncertainty=spectrum.uncertainty, mask=spectrum.mask)

def offset_flux(offset, spectrum):
    """Offset the flux of a spectrum by a given amount

    Parameters
    ----------
    offset : int, float or astropy.units.Quantity
        The amount by which to offset the flux. If an `int` or `float`, then will use the current units of the flux.
    spectrum : specutils.Spectrum1D
        Spectrum to offset

    Returns
    -------
    specutils.Spectrum1D
        Spectrum with flux offset by desired amount
    """
    if not isinstance(offset, u.Quantity):
        offset *= spectrum.flux.unit
    return spectrum + offset

def spectrum_range(spectrum):
    """Returns the minimum and maximum of the flux in the input spectrum

    Parameters
    ----------
    spectrum : specutils.Spectrum1D
        Spectrum for which to return the range

    Returns
    -------
    float
        Minimum of the flux
    float
        Maximum of the flux
    """    
    flux = spectrum.flux
    return np.nanmin(flux).to_value(), np.nanmax(flux).to_value()

def toMarzJSON(simplespectrum):
    """Returns the simplespectrum to a Marz compatible JSON object
    """
    result = {}
    keys = simplespectrum.spectra.keys()
    reduced_keyword = None
    if 'reduced' in keys:
        reduced_keyword = 'reduced'
    sky_keyword = None
    if 'sky' in keys:
        sky_keyword = 'sky'
    #df = simplespectrum.to_dataframe_all()

    #wavelength = df["wavelength"].to_numpy()*10.0   # TODO: Better way than this using Units?
    #print("WAVELENGTH UNITS "+str(simplespectrum.spectra['reduced'].object.data.spectral_axis.unit))
    wavelength = simplespectrum.spectra[reduced_keyword].object.data.spectral_axis.to(u.Unit('Angstrom'))
    wavelength = np.array(wavelength, ndmin=1, copy=False)
    reduced = np.array(simplespectrum.spectra[reduced_keyword].object.data.flux, ndmin=1, copy=False)
    if sky_keyword is not None:
        sky = np.array(simplespectrum.spectra['sky'].object.data.flux, ndmin=1, copy=False)
    else:
        sky = np.full(wavelength.size, None, wavelength.dtype)
    variance = simplespectrum.spectra[reduced_keyword].object.data.uncertainty.array

    reduced = np.where(np.isnan(reduced), None, reduced)
    sky = np.where(np.isnan(sky), None, sky)
    variance = np.where(np.isnan(variance), None, variance)

    last = len(reduced)
    result["wavelength"] = wavelength[0:last].tolist()
    result["intensity"] = reduced[0:last].tolist()
    result["sky"] = sky[0:last].tolist()
    result["variance"] = variance[0:last].tolist()
    return result


def determine_format(data_or_file, *args, **kwargs):
    """Determine which FITS format the input file or HDUList is in.

    Parameters
    ----------
    data_or_file : HDUList, Path or file object
        The FITS data for which to determine the format

    Returns
    -------
    string
        Name of the format
    """
    
    if isinstance(data_or_file, fits.HDUList):
        return format_from_hdulist(data_or_file, *args, **kwargs)

    elif isinstance(data_or_file, (str, Path)) and not os.path.isdir(data_or_file):
        with fits.open(data_or_file) as hdulist:
            return format_from_hdulist(hdulist, *args, **kwargs)


def format_from_hdulist(hdulist, *args, **kwargs):
    info = hdulist.fileinfo(0)
    formats = registry.identify_format('read', SpectrumList, None, None, [hdulist], kwargs)
    return formats[-1] if formats else None # We only want one format

def read_spectra_file(data_or_file, format=None, config_dict=None):
    """Read a FITS file of a spectrum into a SpectrumList

    Parameters
    ----------
    data_or_file : HDUList, Path or file object
        The FITS data to load
    format : str, optional
        The format of the FITS file, by default None
    config_dict : dict, optional
        Used to specify the format of the FITS file, if the format string doesn't work, by default None

    Returns
    -------
    specutils.SpectrumList
        Data from the FITS file in a SpectrumList object
    """    
    if data_or_file is None:
        return None
    
    if format:
        return SpectrumList.read(
            data_or_file,
            format=format
        )
    elif config_dict:
        return SpectrumList.read(
            data_or_file,
            **config_dict
        )
    else:
        return SpectrumList.read(
            data_or_file,
            format=determine_format(data_or_file)
        )

def read_template_file(path_to_file):
    """Read a file containing template spectra

    Parameters
    ----------
    path_to_file : str or Path
        The path pointing to the template file

    Returns
    -------
    list
        List of Spectrum1D objects containing the template spectra
    """    
    template_list = []
    with open(path_to_file) as template_file:
        template_json = json.load(template_file)
        for template in template_json:
            flux = template['spec']
            spectral_axis = np.linspace(template['start_lambda'], template['end_lambda'], len(template['spec']))

            if template['log_linear']:
                spectral_axis = 10. ** spectral_axis

            template_list.append(Spectrum1D(spectral_axis=spectral_axis * u.AA, flux=flux * u.ct, meta={'purpose': template['name']}))

    return template_list