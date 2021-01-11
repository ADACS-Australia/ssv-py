from .. import helpers
from .. import plotting
from .. import utils
import altair as alt
from astropy.table import QTable
import astropy.units as u
from astropy.io import fits
import numpy as np
from specutils import SpectrumList, Spectrum1D
import json
import math
from types import SimpleNamespace
import pandas as pd
from copy import deepcopy
from astropy.nddata import (
    NDData, VarianceUncertainty, StdDevUncertainty, InverseVariance,
)
from astropy.utils.exceptions import AstropyWarning
import warnings

warnings.simplefilter('ignore', category=AstropyWarning)

def read_spectra_file_simple(path_to_file, format=None):
    if format:
        return SpectrumList.read(
            path_to_file,
            format
        )
    else:
        return SpectrumList.read(
            path_to_file,
        )

def read_spectra_file(path_to_file, config_dict):
    return SpectrumList.read(
            path_to_file,
            **config_dict
        )


def read_template_file(path_to_file):
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

class SpectrumIndividual():
    """
    The container class to actually hold the spectrum data

    """
    def __init__(self, purpose, data):
        self.purpose = purpose

        self.data = self._format_spectrum(data)

        self.wavelength_unit = self.data.spectral_axis.unit
        self.flux_unit = self.data.unit

        self.wavelength_offset = 0
        self.flux_offset = 0

        self._wavelength_redshift = 0.0

        self._transform_functions = [utils.id]
        self._show_variance = False

    def _format_spectrum(self, spectrum):
        if isinstance(spectrum.uncertainty, StdDevUncertainty):
            uncertainty = spectrum.uncertainty.quantity
        elif isinstance(spectrum.uncertainty, InverseVariance):
            uncertainty = np.power(spectrum.uncertainty.quantity, -0.5)
        elif isinstance(spectrum.uncertainty, VarianceUncertainty):
            uncertainty = np.power(spectrum.uncertainty.quantity, 0.5)
        else:
            blank_uncertainty = np.full(len(spectrum.spectral_axis), np.nan)
            uncertainty = u.Quantity(blank_uncertainty, unit=spectrum.flux.unit)

        return Spectrum1D(spectral_axis=spectrum.spectral_axis, flux=spectrum.flux, uncertainty=StdDevUncertainty(uncertainty), mask=np.isnan(spectrum.flux), meta=spectrum.meta)

    def show_variance(self, show=True):
        """
        Toggles whether or not the variance should be included when converting to a dataframe

        Parameters
        ---------

        show
            If True, the variance will be included in the dataframe
        """
        self._show_variance = show

    def set_transform_functions(self, *functions):
        """
        Set the transform functions to be applied to every spectrum

        Parameters
        ---------

        *functions
            Transform functions
        """

        self._transform_functions = [utils.id, *functions]

    @property
    def wavelength_redshift(self):
        return self._wavelength_redshift

    @wavelength_redshift.setter
    def wavelength_redshift(self, z):
        self._wavelength_redshift = z

    def to_dataframe(self):
        """
        Applies transform functions and converts the spectrum data to a pandas DataFrame

        Returns
        -------
        pandas.DataFrame
            DataFrame columns for wavelength, flux and (optionally) the variance
        """

        data = utils.compose(*self._transform_functions)(self.data)
        data = utils.offset_flux(self.flux_offset, data)
        data = utils.offset_wavelength(self.wavelength_offset, data)
        original_redshift = 0.0 # TODO: Is this always zero, does it adjust with _wavelength_redshift? the current redshift of this spectrum
        z = self.wavelength_redshift
        fact = (1 + z) / (1 + original_redshift)
        data = utils.redshift_wavelength(fact, data)
        return plotting.convert_spectrum_to_dataframe(
                self.purpose,
                data,
                output_wavelength_unit=u.nm,
                show_variance=self._show_variance
            )

class SimpleSpectrum:
    """
    Controller class to more easily deal with the spectrum data
    Internally uses SpectrumIndividual objects to store pairs of wavelength and flux

    Attributes
    ---------

    name
        The name of the SimpleSpectrum object, used as part of the label in the final chart
    """
    def __init__(self, name, data):
        self.name = name
        self.spectra = self._map_data(data)

    def __str__(self):
        return f"SimpleSpectrum: {self.name}"

    def _map_data(self, data):
        data = deepcopy(data)
        data_dict = {}

        for datum in data:
            data_dict[datum.meta.get('purpose')] = SimpleNamespace(
                object=SpectrumIndividual(purpose=datum.meta.get('purpose'), data=datum),
                visible=True
            )

        return data_dict

    def to_dataframe(self, wavelength_unit=u.nm):
        """
        Create a pandas DataFrame containing all data from the spectra that are set to be visible in the chart

        Parameters
        ---------

        wavelength_unit
            Set the wavelength unit for the resulting DataFrame

        Returns
        -------

        None or pandas.DataFrame
        """
        dfs = [spectrum.object.to_dataframe() for spectrum in self.spectra.values() if spectrum.visible]
        if not dfs:
            return None
        df = pd.concat(dfs).melt(id_vars=[plotting.DEFAULT_SPECTRA_LABEL_COLUMN_NAME, plotting.DEFAULT_WAVELENGTH_COLUMN_NAME], var_name='type', value_name=plotting.DEFAULT_FLUX_COLUMN_NAME)
        df[plotting.DEFAULT_SPECTRA_LABEL_COLUMN_NAME] = df[plotting.DEFAULT_SPECTRA_LABEL_COLUMN_NAME].apply(lambda label: f"{self.name}: {label}")
        return df

    def _toggle_trace(self, trace_key, show):
        spectrum = self.spectra.get(trace_key)
        if spectrum is not None:
            spectrum.visible = show
        else:
            raise KeyError(f"{trace_key} key is not in the spectrum file")

    def set_variance_visible(self, trace_key, visible=True):
        """
        Toggle the visibility of the variance for a spectrum


        Parameters
        ---------

        trace_key
            The key for the SpectrumIndividual object to turn on variance

        visible
            If True, variance will be visible, else it will not be visible
        """

        if trace_key in self.spectra.keys():
            self.spectra[trace_key].object.show_variance(visible)

    def set_visible_traces(self, *trace_keys):
        """
        Toggle which spectra are visible in the final dataframe


        Parameters
        ---------

        *trace_keys
            Keys to set visible. Any key not listed here will be set to not visible
        """

        for key, spectrum in self.spectra.items():
            if key in trace_keys:
                self._toggle_trace(key, True)
            else:
                self._toggle_trace(key, False)

    def set_trace_visible(self, trace_key, visible=True):
        """
        Toggle the visibility of an individual spectrum


        Parameters
        ---------

        trace_key
            The key for the SpectrumIndividual object to toggle visibility

        visible
            If True, the spectrum will be visible, else it will not be visible
        """

        if trace_key in self.spectra.keys():
            self._toggle_trace(trace_key, visible)

    def set_transform_functions(self, trace_keys=None, functions=None):
        """
        Set the list of transform functions to be applied to spectra before the DataFrame is returned

        Parameters
        ---------

        trace_keys
            List of spectra keys to apply the transform functions to
        
        functions
            List of transform functions to apply
        """

        for key, spectrum in self.spectra.items():
            if key in trace_keys:
                self.spectra[key].object.set_transform_functions(*functions)

    def offset_wavelength(self, offset, *trace_keys):
        """
        Set an amount by which to offset the wavelength values of specific spectra

        Parameters
        ---------

        offset
            The amount by which to offset the wavelengths

        *trace_keys
            The keys of the spectra to offset
        """

        if not trace_keys:
            trace_keys = self.spectra.keys()
        
        for trace_key in trace_keys:
            spectrum = self.spectra.get(trace_key)
            if spectrum is not None:
                spectrum.object.wavelength_offset = offset

    def offset_flux(self, offset, *trace_keys):
        """
        Set an amount by which to offset the flux values of specific spectra

        Parameters
        ---------

        offset
            The amount by which to offset the flux

        *trace_keys
            The keys of the spectra to offset
        """

        if not trace_keys:
            trace_keys = self.spectra.keys()
        
        for trace_key in trace_keys:
            spectrum = self.spectra.get(trace_key)
            if spectrum is not None:
                spectrum.object.flux_offset = offset

    def redshift_wavelength(self, z, *trace_keys):
        """
        Set a redshift z to offset the wavelength values of specific spectra

        Parameters
        ---------

        z
            The redshift value used to offest wavelengths

        *trace_keys
            The keys of the spectra to offset
        """

        if not trace_keys:
            trace_keys = self.spectra.keys()
        
        for trace_key in trace_keys:
            spectrum = self.spectra.get(trace_key)
            if spectrum is not None:
                spectrum.object.wavelength_redshift = z

    def flux_range(self, *trace_keys):
        """
        Get dictionary of the maxima and minima of the flux for the specified spectra

        Parameters
        ---------

        *trace_keys
            Spectra keys for which to return the flux range

        Returns
        ------

        Dictionary of the input keys with the flux ranges of the associated spectra
        """

        result = {}
        if not trace_keys:
            trace_keys = self.spectra.keys()
        
        for trace_key in trace_keys:
            spectrum = self.spectra.get(trace_key)
            if spectrum is not None:
                result[trace_key] = utils.spectrum_range(spectrum.object.data)

        return result

    def build_chart(self, wavelength_unit=u.nm, **kwargs):
        """
        Build a altair Chart with the data from the visible spectra

        Parameters
        ---------

        wavelength_unit
            Set the wavelength unit for the data used in the Chart

        **kwargs
            Passed to plotting.plot_spectra

        Returns
        -------

        If no data in the Chart, returns None, else returns an altair Chart of the visible spectra
        """

        chart_data = self.to_dataframe(wavelength_unit=wavelength_unit)
        if chart_data is None:
            return None
        return plotting.plot_spectra(chart_data, **kwargs)