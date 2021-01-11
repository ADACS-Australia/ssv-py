import altair as alt
from astropy.table import QTable
import astropy.units as u
from astropy.io import fits
import numpy as np
from specutils import SpectrumList
import json
import math
from collections import OrderedDict
import pandas as pd

class SimpleSpectrumViewer:
    """
    SimpleSpectrumViewer class description

    """

    def __init__(self, name, wavelength_limits=None):
        self.name = name
        self.x_axis_title = "Wavelength"
        self.y_axis_title = "Flux"
        self.wavelength_unit = u.nm
        self.flux_unit = u.ct
        self.spectrum_dict = OrderedDict()
        self.lines = None
        self._show_grid = True
        self._show_legend = True
        self._chart_properties = {
            "height": 400,
            "width": 700,
            "encoding": {
                "x": {
                    "scale": {"domain": [200,1000]},
                    "type": "quantitative"
                }
            }
        }

    def __str__(self):
        return f'Simple Spectrum Viewer [{self.name}]'

    def set_x_range(self, x_min, x_max):
        """
        Sets the range of the x-axis in the final chart

        Parameters
        ---------

        x_min
            The minimum of the x-axis extent

        x_max
            The maximum of the x-axis extent
        """

        self._chart_properties['encoding']['x'] = {
                    "scale": {"domain": [x_min, x_max]},
                    "type": "quantitative"
                }

    def set_y_range(self, y_min, y_max):
        """
        Sets the range of the y-axis in the final chart

        Parameters
        ---------

        y_min
            The minimum of the y-axis extent

        y_max
            The maximum of the y-axis extent
        """

        self._chart_properties['encoding']['y'] = {
                    "scale": {"domain": [y_min, y_max]},
                    "type": "quantitative"
                }
    
    def add_spectrum(self, spectrum):
        """
        Add a SimpleSpectrum to the viewer

        Parameters
        ---------

        spectrum
            The SimpleSpectrum object to be added to the viewer
        """

        self.spectrum_dict[spectrum.name] = spectrum

    def remove_spectrum(self, spectrum):
        """
        Remove a SimpleSpectrum object from the viewer

        Parameters
        ---------

        spectrum
            Name of the SimpleSpectrum object to be removed from the viewer
        """

        del self.spectrum_dict[spectrum.name]

    def get_all_spectrum_chart_data(self):
        """
        Return a dataframe of wavelengths and spectra for all SimpleSpectrum objects in the viewer

        Returns
        -------

        all_spectrum_chart_data
            Pandas Dataframe containing the wavelength and spectra values for all SimpleSpectrum objects in the viewer

        wavelength_min
            The smallest wavelength value in the dataframe

        wavelength_max
            The largest wavelength value in the dataframe
        """

        all_spectrum_chart_data = pd.concat([spectrum.to_dataframe(wavelength_unit=self.wavelength_unit) for spectrum in self.spectrum_dict.values()])
        wavelength_min = all_spectrum_chart_data['wavelength'].min()
        wavelength_max = all_spectrum_chart_data['wavelength'].max()
        return all_spectrum_chart_data, wavelength_min, wavelength_max

    def add_lines(self, lines):
        """
        Add a SimpleSpectralLines object to the viewer

        Parameters
        ---------

        lines
            The SimpleSpectralLines object to be added to the viewer
        """

        self.lines = lines

    def remove_lines(self):
        """
        Removes SimpleSpectralLines object from the viewer
        """

        self.lines = None

    def set_x_axis_title(self, title):
        """
        Set the x-axis title

        Parameters
        ---------

        title
            Title of the x-axis
        """

        self.x_axis_title = title
    
    def set_y_axis_title(self, title):
        """
        Set the y-axis title

        Parameters
        ---------

        title
            Title of the y-axis
        """

        self.y_axis_title = title

    def _get_wavelength_title(self):
        return f"{self.x_axis_title} ({self.wavelength_unit.to_string()})"

    def _get_flux_title(self):
        return f"{self.y_axis_title} ({self.flux_unit.to_string()})"

    def set_wavelength_unit(self, unit):
        """
        Set the wavelength unit

        Parameters
        ---------

        unit
            Unit to be used for wavelength data
        """

        self.wavelength_unit = unit

    def set_flux_unit(self, unit):
        """
        Set the flux unit

        Parameters
        ---------

        unit
            Unit to be used for flux data
        """

        self.flux_unit = unit

    def show_grid(self, show=True):
        """
        Turn the chart grid on or off

        Parameters
        ---------

        show
            If True, grid will be visible in the chart
        """

        self._show_grid = show

    def show_legend(self, show=True):
        """
        Turn the chart legend on or off

        Parameters
        ---------

        show
            If True, legend will be visible in the chart
        """

        self._show_legend = show

    def set_chart_width_height(self, width=None, height=None):
        """
        Set the width and height of the chart

        Parameters
        ---------

        width
            Width of the chart in pixels

        height
            Height of the chart in pixels
        """

        if width is not None:
            self._chart_properties['width'] = width
        if height is not None:
            self._chart_properties['height'] = height

    def build_chart(self):
        """
        Altair Chart constructed using the selected data from SimpleSpectrum and SimpleSpectralLines objects, along with the set chart parameters

        Returns
        -------

        The Altair Chart
        """
        layer_list = []

        for spectrum in self.spectrum_dict.values():
            layer = spectrum.build_chart(
                wavelength_axis_label=self._get_wavelength_title(),
                flux_axis_label=self._get_flux_title(),
                wavelength_unit=self.wavelength_unit
            )
            if layer is not None:
                layer_list.append(layer)

        if self.lines is not None:
            _, wavelength_min, wavelength_max = self.get_all_spectrum_chart_data()
            self.lines.set_wavelength_limits(wavelength_min=wavelength_min, wavelength_max=wavelength_max)
            layer_list.extend(
                self.lines.build_chart(
                    wavelength_axis_label=self._get_wavelength_title(),
                    wavelength_unit=self.wavelength_unit
                )
            )
        
        base_chart = alt.layer(
                *layer_list,
            ).configure(
                axis=alt.AxisConfig(grid=self._show_grid),
                legend=alt.LegendConfig(disable=not(self._show_legend)),
            ).properties(
                **self._chart_properties
            ).interactive()

        return base_chart