import pytest
from specutils import Spectrum1D, SpectrumList

import ssv.loaders as loaders
from ssv.helpers import (
    smooth_spectra, specutils_spectra_to_table_spectra
)

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


class TestSmooth:

    def test_ozdes_smooth(self, shared_datadir):
        spectra = SpectrumList.read(
            shared_datadir / OZDES_TEST_FILENAME,
            format=loaders.SINGLE_SPLIT_LABEL,
            **OZDES_CONFIG
        )
        smooth_spectra(spectra)


class TestConvert:

    def test_ozdes_no_prefer_combined(self, shared_datadir):
        spectra = SpectrumList.read(
            shared_datadir / OZDES_TEST_FILENAME,
            format=loaders.SINGLE_SPLIT_LABEL,
            **OZDES_CONFIG
        )
        tables = specutils_spectra_to_table_spectra(spectra)
        len(tables) == 5

    def test_ozdes_prefer_combined(self, shared_datadir):
        spectra = SpectrumList.read(
            shared_datadir / OZDES_TEST_FILENAME,
            format=loaders.SINGLE_SPLIT_LABEL,
            **OZDES_CONFIG
        )
        tables = specutils_spectra_to_table_spectra(
            spectra, prefer_combined=True
        )
        len(tables) == 1

    def test_fits2json_1(self, shared_datadir):
        from ssv.viewer import read_spectra_file, read_spectra_file_simple, SimpleSpectrum
        from ssv import utils
        spectrum_file = shared_datadir / "marz/emlLinearVacuumNoHelio.fits"
        formats = loaders.whatformat(spectrum_file)
        if len(formats) > 1:
            loaders.unregister(formats[0])
        spectrum_data = read_spectra_file_simple(spectrum_file)
        loaders.restore_registered_loaders()
        spectrum = SimpleSpectrum('fits2JSON', spectrum_data)
        asjson = utils.toMarzJSON(spectrum)

    def test_fits2json_2(self, shared_datadir):
        from ssv.viewer import read_spectra_file, read_spectra_file_simple, SimpleSpectrum
        from ssv import utils
        spectrum_file = shared_datadir / "marz/spec-4444-55538-1000.fits"
        formats = loaders.whatformat(spectrum_file)
        if len(formats) > 1:
            loaders.unregister(formats[0])
        spectrum_data = read_spectra_file_simple(spectrum_file)
        loaders.restore_registered_loaders()
        spectrum = SimpleSpectrum('fits2JSON', spectrum_data)
        asjson = utils.toMarzJSON(spectrum)
