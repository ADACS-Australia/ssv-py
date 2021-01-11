from .. import plotting
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


class SimpleSpectralLines:
    """
    Class to control spectral lines in the final Chart
    """
    def __init__(self, lines=None):
        """
        Parameters
        ---------

        lines
            Sets which spectral lines are contained in the class, should be an astropy QTable with columns of line name and wavelength
        """

        self.lines = self._get_default_lines() if lines is None else lines

    def __str__(self):
        return f"SimpleSpectralLines"

    def _get_default_lines(self):
        # from adacs spectrum viewer for now, should produce a list that works
        # across more wavelengths
        return QTable({
            "name": [
                'Lyβ',
                'Lyα',
                '[NV]',
                'Si4',
                'CIV',
                'CIII',
                'MgII',
                '[OII]',
                '[NeIII]',
                'K',
                'H',
                'Hδ',
                'G-band',
                'Hγ',
                'Hβ',
                '[OIII]',
                '[OIII]',
                'Mg',
                'Na',
                '[NII]',
                'Hα',
                '[NII]',
                '[SII]',
                '[SII]',
            ],
            "wavelength": [
                102.5722,
                121.5670,
                124.014,
                140.00,
                154.906,
                190.873,
                279.875,
                372.8485,
                386.981,
                393.3663,
                396.8468,
                410.292,
                430.44,
                434.169,
                486.1325,
                495.8911,
                500.6843,
                517.53,
                589.40,
                654.984,
                656.280,
                658.523,
                671.832,
                673.271,
            ] * u.nm,
        })

    def set_wavelength_limits(self, wavelength_min, wavelength_max):
        """
        Set the range within which the spectral lines should be rendered

        Parameters
        ---------

        wavelength_min
            Minimum wavelength value
            
        wavelength_max
            Maximum wavelength value
        """

        self.wavelength_min = wavelength_min
        self.wavelength_max = wavelength_max

    def offset_wavelength(self, offset):
        """
        Set an amount by which to offset the wavelength values of specific spectra

        Parameters
        ---------

        offset
            The amount by which to offset the wavelengths
        """

        offset_axis = self.lines['wavelength']

        if not isinstance(offset, u.Quantity):
            offset *= offset_axis.unit
        
        offset_axis += offset

    def redshift_wavelength(self, z):
        """
        Set a redshift z to offset the wavelength values of specific spectra

        Parameters
        ---------

        z
            The redshift value used to offest wavelengths
        """

        offset_axis = self.lines['wavelength']

        redshift = 0.0 # the current redshift of this spectrum
        fact = (1 + z) / (1 + redshift)
        
        offset_axis *= fact

    def to_dataframe(self, wavelength_unit=u.nm):
        """
        Create a pandas DataFrame containing all data from the spectral lines that are within the wavelength limits

        Parameters
        ---------

        wavelength_unit
            Set the wavelength unit for the resulting DataFrame

        Returns
        -------

        pandas.DataFrame
        """

        lines_df = plotting.convert_lines_to_dataframe(self.lines, output_wavelength_unit=wavelength_unit, wavelength_min=self.wavelength_min, wavelength_max=self.wavelength_max)
        return lines_df

    def build_chart(self, wavelength_unit=u.nm, **kwargs):
        """
        Build an altair Chart with the data from the visible spectra

        Parameters
        ---------

        wavelength_unit
            Set the wavelength unit for the data used in the Chart

        **kwargs
            Passed to plotting.plot_lines

        Returns
        -------

        altair.Chart
        """

        return plotting.plot_lines(
                self.to_dataframe(wavelength_unit=wavelength_unit),
                style={
                    "strokeDash": [5,3],
                },
                **kwargs
            )