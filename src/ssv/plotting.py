import textwrap

import altair as alt
import astropy.units as u
from astropy.table import QTable
import pandas as pd
import numpy as np


# altair limits by default the number of rows to 5000, switch this off
alt.data_transformers.disable_max_rows()

DEFAULT_WAVELENGTH_AXIS_LABEL = "Wavelength"
DEFAULT_FLUX_AXIS_LABEL = "Flux"
DEFAULT_WAVELENGTH_UNIT = u.nm
DEFAULT_LINE_SIZE = 1
DEFAULT_LINE_OPACITY = 0.8
DEFAULT_PLOT_HEIGHT = 400
DEFAULT_PLOT_WIDTH = 400
DEFAULT_SECONDARY_COLOUR = "grey"
DEFAULT_SECONDARY_OPACITY = 0.5

DEFAULT_WAVELENGTH_COLUMN_NAME = "wavelength"
DEFAULT_FLUX_COLUMN_NAME = "flux"
DEFAULT_SPECTRA_LABEL_COLUMN_NAME = "label"
DEFAULT_LINE_LABEL_COLUMN_NAME = "name"
DEFAULT_UNCERTAINTY_COLUMN_NAME = "variance"


def _custom_axis(*, axis, column, title):
    return axis(
        column,
        axis=alt.Axis(
            title=title
        )
    )

def custom_wavelength_axis(
    *,
    wavelength_column=DEFAULT_WAVELENGTH_COLUMN_NAME,
    wavelength_title=DEFAULT_WAVELENGTH_AXIS_LABEL,
):
    """
    Handles some additional shared logic on the wavelength axis

    Parameters
    ---------

    wavelength_column
        Column name for the wavelengths in the spectrum DataFrame
    
    wavelength_title
        Title to be shown on the wavelength axis of the Chart
    """
    return _custom_axis(axis=alt.X, column=wavelength_column, title=wavelength_title)


def custom_flux_axis(
    *,
    flux_column=DEFAULT_FLUX_COLUMN_NAME,
    flux_title=DEFAULT_FLUX_AXIS_LABEL,
):
    """
    Handles some additional shared logic on the spectra/flux axis

    Parameters
    ---------

    flux_column
        Column name for the flux in the spectrum DataFrame
    
    flux_title
        Title to be shown on the flux axis of the Chart
    """
    return _custom_axis(axis=alt.Y, column=flux_column, title=flux_title)

def convert_spectrum_to_dataframe(
    label, spectrum, *,
    output_wavelength_unit=u.nm,
    show_variance=False
):
    """
    Converts the astropy QTable in the required format to the pandas DataFrame
    format used by altair.

    label
        Label used to group the data
    
    spectrum | specutils.Spectrum1D
        Spectrum1D object used to fill out the returned DataFrame
    
    output_wavelength_unit | astropy.Unit
        Converts the wavelength values to be in these units in the returned DataFrame

    show_variance
        Include the column for variance in the output DataFrame

    Returns
    ------

    pandas.DataFrame
        DataFrame with columns for wavelength, flux, variance (optionally) and the label
    """

    wavelength = spectrum.spectral_axis.to_value(
            output_wavelength_unit, u.equivalencies.spectral(),
        )
    flux = spectrum.flux.to_value()

    spectrum_df = pd.DataFrame({
        DEFAULT_WAVELENGTH_COLUMN_NAME: wavelength,
        'spectrum': flux,
    })

    if show_variance and spectrum.uncertainty and not np.all(np.isnan(spectrum.uncertainty.array)): 
        spectrum_df[DEFAULT_UNCERTAINTY_COLUMN_NAME] = spectrum.uncertainty.quantity.to_value()
    
    spectrum_df[DEFAULT_SPECTRA_LABEL_COLUMN_NAME] = label
    return spectrum_df

def plot_spectra(
    spectrum_df, *,
    wavelength_axis_label=DEFAULT_WAVELENGTH_AXIS_LABEL,
    flux_axis_label=DEFAULT_FLUX_AXIS_LABEL,
    wavelength_column=DEFAULT_WAVELENGTH_COLUMN_NAME,
    flux_column=DEFAULT_FLUX_COLUMN_NAME
):
    """
    Plots a spectrum

    Parameters
    ---------

    spectrum_df
        Pandas DataFrame containing the data necessary to plot the spectrum

    wavelength_axis_label
        Title for the wavelength axis of the Chart

    flux_axis_label
        Title for the flux axis of the Chart

    wavelength_column
        Title for the wavelength column in `spectrum_df`

    flux_column
        Title for the flux column in `spectrum_df`

    Returns
    -------

    altair.Chart
        Plot of the spectrum in `spectrum_df`

    """

    spectra_plot = alt.Chart(spectrum_df).mark_line().encode(
        x=custom_wavelength_axis(
            wavelength_column=wavelength_column,
            wavelength_title=wavelength_axis_label
        ),
        y=custom_flux_axis(
            flux_column=flux_column,
            flux_title=flux_axis_label
        ),
        color=DEFAULT_SPECTRA_LABEL_COLUMN_NAME,
        strokeDash='type'
    )

    return spectra_plot

def convert_lines_to_dataframe(lines, wavelength_column=DEFAULT_WAVELENGTH_COLUMN_NAME,
    name_column=DEFAULT_LINE_LABEL_COLUMN_NAME, output_wavelength_unit=DEFAULT_WAVELENGTH_UNIT,
    wavelength_min=0, wavelength_max=1e99):
    """
    Converts the spectral lines astropy QTable in the required format to the pandas DataFrame
    format used by altair.

    lines | astropy.QTable
        QTable must have columns for spectral line name and wavelength
    
    wavelength_column
        Column title for the line wavelengths in `lines`

    name_column
        Column title for the line names in `lines`
    
    output_wavelength_unit | astropy.Unit
        Converts the wavelength values to be in these units in the returned DataFrame

    wavelength_min
        Defines a minimum bound for which spectral lines will be included in the returned DataFrame

    wavelength_max
        Defines a maximum bound for which spectral lines will be included in the returned DataFrame

    Returns
    ------

    pandas.DataFrame
        DataFrame with columns for spectral line names and wavelengths
    """

    lines_df = pd.DataFrame({
        name_column: lines[name_column],
        wavelength_column: lines[wavelength_column].to_value(
            output_wavelength_unit, u.equivalencies.spectral(),
        ),
    })
    lines_df = lines_df[lines_df[wavelength_column] > wavelength_min]
    lines_df = lines_df[lines_df[wavelength_column] < wavelength_max]
    return lines_df

def plot_lines(
    lines_df, *,
    style=None,
    output_wavelength_unit=DEFAULT_WAVELENGTH_UNIT,
    wavelength_column=DEFAULT_WAVELENGTH_COLUMN_NAME,
    name_column=DEFAULT_LINE_LABEL_COLUMN_NAME,
    wavelength_axis_label=DEFAULT_WAVELENGTH_AXIS_LABEL,
):
    """
    Plots spectral lines

    Parameters
    ---------

    lines_df
        Pandas DataFrame containing the data necessary to plot the spectral lines

    style
        Used to style the marks in the Chart

    wavelength_column
        Title for the wavelength column in `lines_df`

    name_column
        Title for the name column in `lines_df`

    wavelength_axis_label
        Title for the wavelength axis of the Chart

    Returns
    -------

    altair.Chart
        Plot of the spectral lines in `lines_df`

    altair.Chart
        Plot used to label the lines with their names

    """
    if style is None:
        style = {}
    style.setdefault("size", DEFAULT_LINE_SIZE)
    style.setdefault("opacity", DEFAULT_LINE_OPACITY)


    plot = alt.Chart(
        lines_df
    ).mark_rule(
        **style,
        clip=True
    ).encode(
        x=custom_wavelength_axis(
            wavelength_column=wavelength_column,
            wavelength_title=wavelength_axis_label
        ),
        tooltip=[name_column, wavelength_column],
    )

    line_names = plot.mark_text(clip=True).encode(text='name:O', y=alt.value(10))

    # Workaround for https://github.com/altair-viz/altair/issues/2009
    # plot = plot.add_selection(alt.selection_single())

    return plot, line_names