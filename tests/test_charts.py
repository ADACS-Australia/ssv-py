import pytest
from specutils import Spectrum1D, SpectrumList

import ssv.loaders as loaders
from ssv.helpers import (
    smooth_spectra, specutils_spectra_to_table_spectra
)
from ssv.viewer import read_spectra_file_simple, read_template_file, \
                            SimpleSpectrum, SimpleSpectralLines, SimpleSpectrumViewer

from specutils.manipulation import box_smooth, gaussian_smooth, trapezoid_smooth, median_smooth

from ssv import utils
import ssv

# Uncomment if running test locally
from pathlib import Path
@pytest.fixture
def shared_datadir():
    return Path("./tests/data")

OZDES_TEST_FILENAME = "OzDES-DR2_04720.fits"
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

def transform_function_array(scaled=False, scaling_max=1, processed=False, subtract_continuum=False, smoothing_function=None, smoothing_width=5):
    array = []
    if scaled:
        array.append(utils.apply_scaling(scaling_max))
    if processed:
        array.append(utils.remove_spurious_points)
    if subtract_continuum:
        array.append(utils.subtract_continuum)
    if smoothing_function:
        array.append(utils.apply_smoothing(smoothing_function, smoothing_width))
    return array

class TestChart1:

    def test_ozdes_smooth(self, shared_datadir):
        spectrum_file = shared_datadir / OZDES_TEST_FILENAME
        formats = ssv.loaders.whatformat(spectrum_file)
        if len(formats) > 1:
            ssv.loaders.unregister(formats[0])
        spectrum_data = read_spectra_file_simple(spectrum_file)
        ssv.loaders.restore_registered_loaders()

        show_sky = True
        show_templates = True
        show_variance = True
        show_processed_data = True
        show_continuum_subtracted = False
        show_scaled_spectra = True
        scaling_max = 1.0
        smoothing_choice = box_smooth
        smoothing_width = 3

        spectrum = SimpleSpectrum('Test SSV', spectrum_data)
        spectrum.set_visible_traces('reduced')
        flux_range = spectrum.flux_range('reduced')
        spectrum.set_variance_visible('reduced', show_variance)
        spectrum.set_trace_visible('sky', show_sky)
        spectrum.set_transform_functions('reduced', transform_function_array(scaled=show_scaled_spectra, scaling_max=scaling_max, processed=show_processed_data, subtract_continuum=show_continuum_subtracted, smoothing_function=smoothing_choice, smoothing_width=smoothing_width))
        spectrum.set_transform_functions('sky', transform_function_array(scaled=show_scaled_spectra, scaling_max=scaling_max, smoothing_function=smoothing_choice, smoothing_width=smoothing_width))


        best_fit_redshift = 0
        best_fit_template = 'Quasar'
        best_fit_template_index = 12

        lines = SimpleSpectralLines()
        lines.redshift_wavelength(best_fit_redshift)

        templatedir = Path('./tests/data/marz/')
        TEMPLATE_FILENAME = 'MarzTemplates.json'
        template_data = read_template_file(templatedir / TEMPLATE_FILENAME)
        templates = SimpleSpectrum('Templates', template_data)
        template_choice = 'Quasar'
        templates.set_visible_traces(template_choice)
        templates.set_trace_visible(template_choice, show_templates)
        templates.redshift_wavelength(best_fit_redshift, best_fit_template)
        templates.set_transform_functions(template_choice, transform_function_array(scaled=show_scaled_spectra, scaling_max=scaling_max, processed=show_processed_data, 
            subtract_continuum=show_continuum_subtracted, smoothing_function=smoothing_choice, smoothing_width=smoothing_width))

        viewer = SimpleSpectrumViewer('Simple')
        viewer.add_spectrum(spectrum)
        viewer.add_spectrum(templates)
        viewer.add_lines(lines)
        viewer.show_grid(True)
        viewer.show_legend(True)
        viewer.set_chart_width_height(height=500)

        viewer.build_chart()