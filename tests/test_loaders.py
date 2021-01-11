import pytest

# Uncomment if running test locally
from pathlib import Path
@pytest.fixture
def shared_datadir():
    return Path("./tests/data")

from astropy.nddata import (
    VarianceUncertainty, StdDevUncertainty, InverseVariance,
)
import astropy.units as u
from specutils import Spectrum1D, SpectrumList

import ssv.loaders as loaders

DC_6DFGS_TEST_FILENAMES = [
    "all-c0022498-344732r_spectrum0.fits",
    "all-c0022498-344732v_spectrum0.fits",
]
GALAH_TEST_FILENAME = "1704070052012671.fits"
GAMA_2QZ_TEST_FILENAME = "J113606.3+001155a.fit"
GAMA_2SLAQ_QSO_TEST_FILENAME = "J091726.21+003424.0_a14_040423.fit"
GAMA_2DFGRS_TEST_FILENAME = "299869.fit"
GAMA_GAMA_LT_TEST_FILENAME = "LTF_09_1128_0555.fit"
GAMA_GAMA_TEST_FILENAME = "G12_Y2_009_044.fit"
GAMA_MGC_TEST_FILENAME = "MGC23320.fit"
GAMA_SDSS_TEST_FILENAME = "spSpec-51609-0304-087.fit"
GAMA_WIGGLEZ_TEST_FILENAME = "Spectrum-195254.fit"
OZDES_TEST_FILENAME = "OzDES-DR2_04720.fits"
SIXDFGS_TABLE_TEST_FILENAME = "1D-c0022498-344732.fits"
SIXDFGS_COMBINED_TEST_FILENAME = "all-g2359555-392832.fits"
TWODFGRS_TEST_FILENAME = "000001.fits"
AAOMEGA_WITH_RWSS = "OBJ0039red.fits"
AAOMEGA_WITHOUT_RWSS = "OBJ0032red.fits"

GALAH_CONFIG = {
    "hdus": {
        "0": {"purpose": "science"},
        "1": {"purpose": "error_stdev"},
        "2": {"purpose": "unreduced_science"},
        "3": {"purpose": "unreduced_error_stdev"},
        "4": {
            "purpose": "science",
            "units": {"flux_unit": ""},
        },
    },
    "wcs": {
        "pixel_reference_point_keyword": "CRPIX1",
        "pixel_reference_point_value_keyword": "CRVAL1",
        "pixel_width_keyword": "CDELT1",
        "wavelength_unit": "Angstrom",
    },
    "units": {"flux_unit": "count"},
    "all_standard_units": False,
    "all_keywords": False,
    "valid_wcs": False,
}
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
GAMA_2QZ_CONFIG = {
    "hdus": None,
    "wcs": {
        "pixel_reference_point_keyword": "CRPIX1",
        "pixel_reference_point_value_keyword": "CRVAL1",
        "pixel_width_keyword": "CD1_1",
        "wavelength_unit": "Angstrom",
    },
    "units": {"flux_unit": "count"},
    "all_standard_units": False,
    "all_keywords": True,
    "valid_wcs": False,
}
GAMA_2SLAQ_QSO_CONFIG = {
    "hdus": None,
    "wcs": {
        "pixel_reference_point_keyword": "CRPIX1",
        "pixel_reference_point_value_keyword": "CRVAL1",
        "pixel_width_keyword": "CDELT1",
        "wavelength_unit": "Angstrom",
    },
    "units": {"flux_unit": "count"},
    "all_standard_units": False,
    "all_keywords": True,
    "valid_wcs": False,
}
GAMA_LT_CONFIG = {
    "hdus": {"0": {"purpose": "science"},},
    "wcs": {
        "pixel_reference_point_keyword": "CRPIX",
        "pixel_reference_point_value_keyword": "CRVAL",
        "pixel_width_keyword": "CDELT",
        "wavelength_unit": "Angstrom",
    },
    "units": {"flux_unit": "count"},
    "all_standard_units": False,
    "all_keywords": False,
    "valid_wcs": False,
}
GAMA_WIGGLEZ_CONFIG = {
    "hdus": {
        "0": {"purpose": "science"},
        "1": {"purpose": "error_variance"},
        "2": {"purpose": "skip"},
    },
    "wcs": {
        "pixel_reference_point_keyword": "CRPIX1",
        "pixel_reference_point_value_keyword": "CRVAL1",
        "pixel_width_keyword": "CDELT1",
        "wavelength_unit": "Angstrom",
    },
    "units": {"flux_unit": "count"},
    "all_standard_units": False,
    "all_keywords": False,
    "valid_wcs": False,
}
GAMA_2DFGRS_CONFIG = {
    "hdu": None,
    "wcs": {
        "pixel_reference_point_keyword": "CRPIX1",
        "pixel_reference_point_value_keyword": "CRVAL1",
        "pixel_width_keyword": "CDELT1",
        "wavelength_unit": "Angstrom",
    },
    "units": {"flux_unit": "count"},
    "all_standard_units": False,
    "all_keywords": True,
    "valid_wcs": False,
}
GAMA_GAMA_CONFIG = {
    "hdu": {
        "1": {
            "purpose": "science",
            "units": {"flux_unit": "10^-17 erg/s/cm^2/A"},
        },
        "2": {"purpose": "error_stdev"},
        "3": {"purpose": "unreduced_science"},
        "4": {"purpose": "unreduced_error_stdev"},
        "5": {"purpose": "sky"},
    },
    "wcs": {
        "pixel_reference_point_keyword": "CRPIX1",
        "pixel_reference_point_value_keyword": "CRVAL1",
        "pixel_width_keyword": "CD1_1",
        "wavelength_unit": "Angstrom",
    },
    "units": {"flux_unit": "count"},
    "all_standard_units": False,
    "all_keywords": False,
    "valid_wcs": False,
}
GAMA_MGC_CONFIG = {
    "hdu": None,
    "wcs": {
        "pixel_reference_point_keyword": "CRPIX1",
        "pixel_reference_point_value_keyword": "CRVAL1",
        "pixel_width_keyword": "CD1_1",
        "wavelength_unit": "Angstrom",
    },
    "units": {"flux_unit": "count"},
    "all_standard_units": False,
    "all_keywords": True,
    "valid_wcs": False,
}

class TestSingleSplit:
    def test_galah(self, shared_datadir):
        spectra = SpectrumList.read(
            shared_datadir / GALAH_TEST_FILENAME,
            format=loaders.SINGLE_SPLIT_LABEL,
            **GALAH_CONFIG
        )
        # Should be main spectra, without sky, and normalised
        assert len(spectra) == 3

        assert spectra[0].flux.unit == u.count
        assert spectra[1].flux.unit == u.Unit('') # dimensionless
        assert spectra[2].flux.unit == u.count

        assert spectra[0].spectral_axis.unit == u.Angstrom
        assert spectra[1].spectral_axis.unit == u.Angstrom
        assert spectra[2].spectral_axis.unit == u.Angstrom

        assert isinstance(spectra[0].uncertainty, StdDevUncertainty)
        assert spectra[1].uncertainty is None
        assert isinstance(spectra[2].uncertainty, StdDevUncertainty)

        assert spectra[0].meta.get("label") is not None
        assert spectra[0].meta.get("header") is not None

        assert spectra[1].meta.get("label") is None
        assert spectra[1].meta.get("header") is not None

        assert spectra[2].meta.get("label") is not None
        assert spectra[2].meta.get("header") is not None

    @pytest.mark.skip(reason="Galah needs fixing")
    def test_galah_guess(self, shared_datadir):
        spectra = SpectrumList.read(shared_datadir / GALAH_TEST_FILENAME)
        # Should be main spectra, without sky, and normalised
        assert len(spectra) == 3

        assert spectra[0].flux.unit == u.count
        assert spectra[1].flux.unit == u.Unit('') # dimensionless
        assert spectra[2].flux.unit == u.count

        assert spectra[0].spectral_axis.unit == u.Angstrom
        assert spectra[1].spectral_axis.unit == u.Angstrom
        assert spectra[2].spectral_axis.unit == u.Angstrom

        assert isinstance(spectra[0].uncertainty, StdDevUncertainty)
        assert spectra[1].uncertainty is None
        assert isinstance(spectra[2].uncertainty, StdDevUncertainty)

        assert spectra[0].meta.get("label") is not None
        assert spectra[0].meta.get("header") is not None

        assert spectra[1].meta.get("label") is None
        assert spectra[1].meta.get("header") is not None

        assert spectra[2].meta.get("label") is not None
        assert spectra[2].meta.get("header") is not None

    def test_ozdes(self, shared_datadir):
        spectra = SpectrumList.read(
            shared_datadir / OZDES_TEST_FILENAME,
            format=loaders.SINGLE_SPLIT_LABEL,
            **OZDES_CONFIG
        )

        # The test file has the combined obs, and 4 other sets
        assert len(spectra) == 5

        assert spectra[0].flux.unit == u.count / u.Angstrom
        assert spectra[0].spectral_axis.unit == u.Angstrom
        assert isinstance(spectra[0].uncertainty, VarianceUncertainty)
        assert spectra[0].meta.get("label") is None
        assert spectra[0].meta.get("header") is not None

        assert spectra[1].flux.unit == u.count / u.Angstrom
        assert spectra[1].spectral_axis.unit == u.Angstrom
        assert isinstance(spectra[1].uncertainty, VarianceUncertainty)
        assert spectra[1].meta.get("label") is not None
        assert spectra[1].meta.get("header") is not None

        assert spectra[2].flux.unit == u.count / u.Angstrom
        assert spectra[2].spectral_axis.unit == u.Angstrom
        assert isinstance(spectra[2].uncertainty, VarianceUncertainty)
        assert spectra[2].meta.get("label") is not None
        assert spectra[2].meta.get("header") is not None

        assert spectra[3].flux.unit == u.count / u.Angstrom
        assert spectra[3].spectral_axis.unit == u.Angstrom
        assert isinstance(spectra[3].uncertainty, VarianceUncertainty)
        assert spectra[3].meta.get("label") is not None
        assert spectra[3].meta.get("header") is not None

        assert spectra[4].flux.unit == u.count / u.Angstrom
        assert spectra[4].spectral_axis.unit == u.Angstrom
        assert isinstance(spectra[4].uncertainty, VarianceUncertainty)
        assert spectra[4].meta.get("label") is not None
        assert spectra[4].meta.get("header") is not None

    def test_ozdes_by_label(self, shared_datadir):
        spectra = SpectrumList.read(
            shared_datadir / OZDES_TEST_FILENAME,
            format="OzDES"
        )

        # The test file has the combined obs, and 4 other sets
        assert len(spectra) == 5

        assert spectra[0].flux.unit == u.count / u.Angstrom
        assert spectra[0].spectral_axis.unit == u.Angstrom
        assert isinstance(spectra[0].uncertainty, VarianceUncertainty)
        assert spectra[0].meta.get("label") is None
        assert spectra[0].meta.get("header") is not None

        assert spectra[1].flux.unit == u.count / u.Angstrom
        assert spectra[1].spectral_axis.unit == u.Angstrom
        assert isinstance(spectra[1].uncertainty, VarianceUncertainty)
        assert spectra[1].meta.get("label") is not None
        assert spectra[1].meta.get("header") is not None

        assert spectra[2].flux.unit == u.count / u.Angstrom
        assert spectra[2].spectral_axis.unit == u.Angstrom
        assert isinstance(spectra[2].uncertainty, VarianceUncertainty)
        assert spectra[2].meta.get("label") is not None
        assert spectra[2].meta.get("header") is not None

        assert spectra[3].flux.unit == u.count / u.Angstrom
        assert spectra[3].spectral_axis.unit == u.Angstrom
        assert isinstance(spectra[3].uncertainty, VarianceUncertainty)
        assert spectra[3].meta.get("label") is not None
        assert spectra[3].meta.get("header") is not None

        assert spectra[4].flux.unit == u.count / u.Angstrom
        assert spectra[4].spectral_axis.unit == u.Angstrom
        assert isinstance(spectra[4].uncertainty, VarianceUncertainty)
        assert spectra[4].meta.get("label") is not None
        assert spectra[4].meta.get("header") is not None

    # This picks up the wcs1d-fits format
    @pytest.mark.xfail(reason="Format is ambiguous")
    def test_ozdes_guess(self, shared_datadir):
        spectra = SpectrumList.read(shared_datadir / OZDES_TEST_FILENAME)

        # The test file has the combined obs, and 4 other sets
        assert len(spectra) == 5

        assert spectra[0].flux.unit == u.count / u.Angstrom
        assert spectra[0].spectral_axis.unit == u.Angstrom
        assert isinstance(spectra[0].uncertainty, VarianceUncertainty)
        assert spectra[0].meta.get("label") is None
        assert spectra[0].meta.get("header") is not None

        assert spectra[1].flux.unit == u.count / u.Angstrom
        assert spectra[1].spectral_axis.unit == u.Angstrom
        assert isinstance(spectra[1].uncertainty, VarianceUncertainty)
        assert spectra[1].meta.get("label") is not None
        assert spectra[1].meta.get("header") is not None

        assert spectra[2].flux.unit == u.count / u.Angstrom
        assert spectra[2].spectral_axis.unit == u.Angstrom
        assert isinstance(spectra[2].uncertainty, VarianceUncertainty)
        assert spectra[2].meta.get("label") is not None
        assert spectra[2].meta.get("header") is not None

        assert spectra[3].flux.unit == u.count / u.Angstrom
        assert spectra[3].spectral_axis.unit == u.Angstrom
        assert isinstance(spectra[3].uncertainty, VarianceUncertainty)
        assert spectra[3].meta.get("label") is not None
        assert spectra[3].meta.get("header") is not None

        assert spectra[4].flux.unit == u.count / u.Angstrom
        assert spectra[4].spectral_axis.unit == u.Angstrom
        assert isinstance(spectra[4].uncertainty, VarianceUncertainty)
        assert spectra[4].meta.get("label") is not None
        assert spectra[4].meta.get("header") is not None

    def test_gama_2qz(self, shared_datadir):
        spectra = SpectrumList.read(
            shared_datadir / GAMA_2QZ_TEST_FILENAME,
            format=loaders.SINGLE_SPLIT_LABEL,
            **GAMA_2QZ_CONFIG
        )
        assert len(spectra) == 1

        assert spectra[0].flux.unit == u.count
        assert spectra[0].spectral_axis.unit == u.Angstrom

        assert isinstance(spectra[0].uncertainty, VarianceUncertainty)

        assert spectra[0].meta.get("label") is not None
        assert spectra[0].meta.get("header") is not None

    @pytest.mark.xfail(reason="Format is ambiguous")
    def test_gama_2qz_guess(self, shared_datadir):
        spectra = SpectrumList.read(
            shared_datadir / GAMA_2QZ_TEST_FILENAME,
        )
        assert len(spectra) == 1

        assert spectra[0].flux.unit == u.count
        assert spectra[0].spectral_axis.unit == u.Angstrom

        assert isinstance(spectra[0].uncertainty, VarianceUncertainty)

        assert spectra[0].meta.get("label") is not None
        assert spectra[0].meta.get("header") is not None

    def test_2slaq_qso(self, shared_datadir):
        spectra = SpectrumList.read(
            shared_datadir / GAMA_2SLAQ_QSO_TEST_FILENAME,
            format=loaders.SINGLE_SPLIT_LABEL,
            **GAMA_2SLAQ_QSO_CONFIG
        )
        assert len(spectra) == 1

        assert spectra[0].flux.unit == u.count
        assert spectra[0].spectral_axis.unit == u.Angstrom

        assert isinstance(spectra[0].uncertainty, VarianceUncertainty)

        assert spectra[0].meta.get("label") is not None
        assert spectra[0].meta.get("header") is not None

    @pytest.mark.xfail(reason="Format is ambiguous")
    def test_2slaq_qso_guess(self, shared_datadir):
        spectra = SpectrumList.read(
            shared_datadir / GAMA_2SLAQ_QSO_TEST_FILENAME,
        )
        assert len(spectra) == 1

        assert spectra[0].flux.unit == u.count
        assert spectra[0].spectral_axis.unit == u.Angstrom

        assert isinstance(spectra[0].uncertainty, VarianceUncertainty)

        assert spectra[0].meta.get("label") is not None
        assert spectra[0].meta.get("header") is not None

    def test_gama_lt(self, shared_datadir):
        spectra = SpectrumList.read(
            shared_datadir / GAMA_GAMA_LT_TEST_FILENAME,
            format=loaders.SINGLE_SPLIT_LABEL,
            **GAMA_LT_CONFIG
        )
        assert len(spectra) == 1

        assert spectra[0].flux.unit == u.count
        assert spectra[0].spectral_axis.unit == u.Angstrom

        assert spectra[0].uncertainty is None

        assert spectra[0].meta.get("label") is not None
        assert spectra[0].meta.get("header") is not None

    @pytest.mark.xfail(reason="Format is ambiguous")
    def test_gama_lt_guess(self, shared_datadir):
        spectra = SpectrumList.read(
            shared_datadir / GAMA_GAMA_LT_TEST_FILENAME,
        )
        assert len(spectra) == 1

        assert spectra[0].flux.unit == u.count
        assert spectra[0].spectral_axis.unit == u.Angstrom

        assert spectra[0].uncertainty is None

        assert spectra[0].meta.get("label") is not None
        assert spectra[0].meta.get("header") is not None

    def test_gama_wigglez(self, shared_datadir):
        spectra = SpectrumList.read(
            shared_datadir / GAMA_WIGGLEZ_TEST_FILENAME,
            format=loaders.SINGLE_SPLIT_LABEL,
            **GAMA_WIGGLEZ_CONFIG
        )
        assert len(spectra) == 1

        assert spectra[0].flux.unit == u.count
        assert spectra[0].spectral_axis.unit == u.Angstrom

        assert isinstance(spectra[0].uncertainty, VarianceUncertainty)

        assert spectra[0].meta.get("label") is not None
        assert spectra[0].meta.get("header") is not None

    @pytest.mark.xfail(reason="Format is ambiguous")
    def test_gama_wigglez_guess(self, shared_datadir):
        spectra = SpectrumList.read(
            shared_datadir / GAMA_WIGGLEZ_TEST_FILENAME,
        )
        assert len(spectra) == 1

        assert spectra[0].flux.unit == u.count
        assert spectra[0].spectral_axis.unit == u.Angstrom

        assert isinstance(spectra[0].uncertainty, VarianceUncertainty)

        assert spectra[0].meta.get("label") is not None
        assert spectra[0].meta.get("header") is not None


class TestMultilineSingle:

    def test_gama_2dfgrs(self, shared_datadir):
        spectra = SpectrumList.read(
            shared_datadir / GAMA_2DFGRS_TEST_FILENAME,
            format=loaders.MULTILINE_SINGLE_LABEL,
            **GAMA_2DFGRS_CONFIG
        )
        assert len(spectra) == 2

        assert spectra[0].flux.unit == u.count
        assert spectra[0].spectral_axis.unit == u.Angstrom
        assert spectra[1].flux.unit == u.count
        assert spectra[1].spectral_axis.unit == u.Angstrom

        assert isinstance(spectra[0].uncertainty, VarianceUncertainty)
        assert spectra[1].uncertainty is None

        assert spectra[0].meta.get("label") is not None
        assert spectra[0].meta.get("header") is not None
        assert spectra[1].meta.get("label") is not None
        assert spectra[1].meta.get("header") is not None

    @pytest.mark.xfail(reason="Format is ambiguous")
    def test_gama_2dfgrs_guess(self, shared_datadir):
        spectra = SpectrumList.read(
            shared_datadir / GAMA_2DFGRS_TEST_FILENAME,
        )
        assert len(spectra) == 2

        assert spectra[0].flux.unit == u.count
        assert spectra[0].spectral_axis.unit == u.Angstrom
        assert spectra[1].flux.unit == u.count
        assert spectra[1].spectral_axis.unit == u.Angstrom

        assert isinstance(spectra[0].uncertainty, VarianceUncertainty)
        assert spectra[1].uncertainty is None

        assert spectra[0].meta.get("label") is not None
        assert spectra[0].meta.get("header") is not None
        assert spectra[1].meta.get("label") is not None
        assert spectra[1].meta.get("header") is not None

    def test_gama_gama(self, shared_datadir):
        spectra = SpectrumList.read(
            shared_datadir / GAMA_GAMA_TEST_FILENAME,
            format=loaders.MULTILINE_SINGLE_LABEL,
            **GAMA_GAMA_CONFIG
        )
        assert len(spectra) == 3

        assert spectra[0].flux.unit == u.Unit("10^-17 erg/s/cm^2/A")
        assert spectra[0].spectral_axis.unit == u.Angstrom
        assert spectra[1].flux.unit == u.count
        assert spectra[1].spectral_axis.unit == u.Angstrom
        assert spectra[2].flux.unit == u.count
        assert spectra[2].spectral_axis.unit == u.Angstrom

        assert isinstance(spectra[0].uncertainty, StdDevUncertainty)
        assert isinstance(spectra[1].uncertainty, StdDevUncertainty)
        assert spectra[2].uncertainty is None

        assert spectra[0].meta.get("label") is not None
        assert spectra[0].meta.get("header") is not None
        assert spectra[1].meta.get("label") is not None
        assert spectra[1].meta.get("header") is not None
        assert spectra[2].meta.get("label") is not None
        assert spectra[2].meta.get("header") is not None

    def test_gama_gama_by_label(self, shared_datadir):
        spectra = SpectrumList.read(
            shared_datadir / GAMA_GAMA_TEST_FILENAME,
            format="GAMA"
        )
        assert len(spectra) == 3

        assert spectra[0].flux.unit == u.Unit("10^-17 erg/s/cm^2/A")
        assert spectra[0].spectral_axis.unit == u.Angstrom
        assert spectra[1].flux.unit == u.count
        assert spectra[1].spectral_axis.unit == u.Angstrom
        assert spectra[2].flux.unit == u.count
        assert spectra[2].spectral_axis.unit == u.Angstrom

        assert isinstance(spectra[0].uncertainty, StdDevUncertainty)
        assert isinstance(spectra[1].uncertainty, StdDevUncertainty)
        assert spectra[2].uncertainty is None

        assert spectra[0].meta.get("label") is not None
        assert spectra[0].meta.get("header") is not None
        assert spectra[1].meta.get("label") is not None
        assert spectra[1].meta.get("header") is not None
        assert spectra[2].meta.get("label") is not None
        assert spectra[2].meta.get("header") is not None

    def test_gama_gama_guess(self, shared_datadir):
        spectra = SpectrumList.read(shared_datadir / GAMA_GAMA_TEST_FILENAME)
        assert len(spectra) == 3

        assert spectra[0].flux.unit == u.Unit("10^-17 erg/s/cm^2/A")
        assert spectra[0].spectral_axis.unit == u.Angstrom
        assert spectra[1].flux.unit == u.count
        assert spectra[1].spectral_axis.unit == u.Angstrom
        assert spectra[2].flux.unit == u.count
        assert spectra[2].spectral_axis.unit == u.Angstrom

        assert isinstance(spectra[0].uncertainty, StdDevUncertainty)
        assert isinstance(spectra[1].uncertainty, StdDevUncertainty)
        assert spectra[2].uncertainty is None

        assert spectra[0].meta.get("label") is not None
        assert spectra[0].meta.get("header") is not None
        assert spectra[1].meta.get("label") is not None
        assert spectra[1].meta.get("header") is not None
        assert spectra[2].meta.get("label") is not None
        assert spectra[2].meta.get("header") is not None

    def test_gama_mgc(self, shared_datadir):
        spectra = SpectrumList.read(
            shared_datadir / GAMA_MGC_TEST_FILENAME,
            format=loaders.MULTILINE_SINGLE_LABEL,
            **GAMA_MGC_CONFIG
        )
        assert len(spectra) == 2

        assert spectra[0].flux.unit == u.count
        assert spectra[0].spectral_axis.unit == u.Angstrom
        assert spectra[1].flux.unit == u.count
        assert spectra[1].spectral_axis.unit == u.Angstrom

        assert isinstance(spectra[0].uncertainty, StdDevUncertainty)
        assert spectra[1].uncertainty is None

        assert spectra[0].meta.get("label") is not None
        assert spectra[0].meta.get("header") is not None
        assert spectra[1].meta.get("label") is not None
        assert spectra[1].meta.get("header") is not None

    def test_gama_mgc_guess(self, shared_datadir):
        spectra = SpectrumList.read(shared_datadir / GAMA_MGC_TEST_FILENAME)
        assert len(spectra) == 2

        assert spectra[0].flux.unit == u.count
        assert spectra[0].spectral_axis.unit == u.Angstrom
        assert spectra[1].flux.unit == u.count
        assert spectra[1].spectral_axis.unit == u.Angstrom

        assert isinstance(spectra[0].uncertainty, StdDevUncertainty)
        assert spectra[1].uncertainty is None

        assert spectra[0].meta.get("label") is not None
        assert spectra[0].meta.get("header") is not None
        assert spectra[1].meta.get("label") is not None
        assert spectra[1].meta.get("header") is not None


class TestAAOmega2dF:
    def test_with_rwss(self, shared_datadir):
        spectra = SpectrumList.read(
            shared_datadir / AAOMEGA_WITH_RWSS,
            format="Data Central AAOmega",
        )
        assert len(spectra) == 139
        for spec in spectra:
            assert spec.meta.get("label") is not None
            assert spec.meta.get("header") is not None
            assert spec.meta.get("purpose") is not None
            assert spec.meta.get("fibre_index") is not None

    def test_without_rwss(self, shared_datadir):
        spectra = SpectrumList.read(
            shared_datadir / AAOMEGA_WITHOUT_RWSS,
            format="Data Central AAOmega",
        )
        assert len(spectra) == 153
        for spec in spectra:
            assert spec.meta.get("label") is not None
            assert spec.meta.get("header") is not None
            assert spec.meta.get("purpose") is not None
            assert spec.meta.get("fibre_index") is not None

    def test_with_rwss_guess(self, shared_datadir):
        spectra = SpectrumList.read(
            shared_datadir / AAOMEGA_WITH_RWSS,
        )
        assert len(spectra) == 139
        for spec in spectra:
            assert spec.meta.get("label") is not None
            assert spec.meta.get("header") is not None
            assert spec.meta.get("purpose") is not None
            assert spec.meta.get("fibre_index") is not None

    def test_without_rwss_guess(self, shared_datadir):
        spectra = SpectrumList.read(
            shared_datadir / AAOMEGA_WITHOUT_RWSS,
        )
        assert len(spectra) == 153
        for spec in spectra:
            assert spec.meta.get("label") is not None
            assert spec.meta.get("header") is not None
            assert spec.meta.get("purpose") is not None
            assert spec.meta.get("fibre_index") is not None


class TestObsCore:
    @pytest.mark.skip(reason="Missing test file")
    def test_gama_sdss(self, shared_datadir):
        spectra = SpectrumList.read(
            shared_datadir / GAMA_SDSS_TEST_FILENAME,
            format="GAMA-SDSS obscore",
        )
        spec = spectra[0]
        obscore = spec.meta.get("obscore")
        assert obscore is not None
        assert obscore.get("calib_level") == 2
        assert obscore.get("t_xel") == 1

    def test_2dfgrs(self, shared_datadir):
        spectra = SpectrumList.read(
            shared_datadir / TWODFGRS_TEST_FILENAME, format="2dFGRS obscore",
        )
        spec = spectra[0]
        obscore = spec.meta.get("obscore")
        assert obscore is not None
        assert obscore.get("calib_level") == 2
        assert obscore.get("t_xel") == 1

    @pytest.mark.skip(reason="Galah needs fixing")
    def test_galah(self, shared_datadir):
        spectra = SpectrumList.read(
            shared_datadir / GALAH_TEST_FILENAME, format="GALAH obscore",
        )
        spec = spectra[0]
        obscore = spec.meta.get("obscore")
        assert obscore is not None
        assert obscore.get("calib_level") == 2
        assert obscore.get("t_xel") == 1

    def test_ozdes(self, shared_datadir):
        spectra = SpectrumList.read(
            shared_datadir / OZDES_TEST_FILENAME, format="OzDES obscore",
        )
        spec = spectra[0]
        obscore = spec.meta.get("obscore")
        assert obscore is not None
        assert obscore.get("calib_level") == 3
        assert obscore.get("t_xel") == 4

    def test_aaomega_with_rwss(self, shared_datadir):
        spectra = SpectrumList.read(
            shared_datadir / AAOMEGA_WITH_RWSS,
            format="Data Central AAOmega obscore",
        )
        for spec in spectra:
            if spec.meta["purpose"] == "reduced":
                obscore = spec.meta.get("obscore")
                assert obscore is not None
                assert obscore.get("calib_level") == 2
                assert obscore.get("t_xel") == 1

    def test_aaomega_without_rwss(self, shared_datadir):
        spectra = SpectrumList.read(
            shared_datadir / AAOMEGA_WITHOUT_RWSS,
            format="Data Central AAOmega obscore",
        )
        for spec in spectra:
            if spec.meta["purpose"] == "reduced":
                obscore = spec.meta.get("obscore")
                assert obscore is not None
                assert obscore.get("calib_level") == 2
                assert obscore.get("t_xel") == 1

    def test_gama_2qz(self, shared_datadir):
        spectra = SpectrumList.read(
            shared_datadir / GAMA_2QZ_TEST_FILENAME,
            format="GAMA-2QZ obscore",
        )
        spec = spectra[0]
        obscore = spec.meta.get("obscore")
        assert obscore is not None
        assert obscore.get("calib_level") == 2
        assert obscore.get("t_xel") == 1

    @pytest.mark.skip(reason="Bad file")
    def test_2slaq_qso(self, shared_datadir):
        spectra = SpectrumList.read(
            shared_datadir / GAMA_2SLAQ_QSO_TEST_FILENAME,
            format="GAMA-2SLAQ-QSO obscore",
        )
        spec = spectra[0]
        obscore = spec.meta.get("obscore")
        assert obscore is not None
        assert obscore.get("calib_level") == 2
        assert obscore.get("t_xel") == 1

    def test_gama_lt(self, shared_datadir):
        spectra = SpectrumList.read(
            shared_datadir / GAMA_GAMA_LT_TEST_FILENAME,
            format="GAMA-LT obscore",
        )
        spec = spectra[0]
        obscore = spec.meta.get("obscore")
        assert obscore is not None
        assert obscore.get("calib_level") == 2
        assert obscore.get("t_xel") == 1

    def test_gama(self, shared_datadir):
        spectra = SpectrumList.read(
            shared_datadir / GAMA_GAMA_TEST_FILENAME,
            format="GAMA obscore",
        )
        spec = spectra[0]
        obscore = spec.meta.get("obscore")
        assert obscore is not None
        assert obscore.get("calib_level") == 2
        assert obscore.get("t_xel") == 1

    def test_gama_mgc(self, shared_datadir):
        spectra = SpectrumList.read(
            shared_datadir / GAMA_MGC_TEST_FILENAME,
            format="GAMA-MGC obscore",
        )
        spec = spectra[0]
        obscore = spec.meta.get("obscore")
        assert obscore is not None
        assert obscore.get("calib_level") == 2
        assert obscore.get("t_xel") == 1

    def test_gama_wigglez(self, shared_datadir):
        spectra = SpectrumList.read(
            shared_datadir / GAMA_WIGGLEZ_TEST_FILENAME,
            format="GAMA-WiggleZ obscore",
        )
        spec = spectra[0]
        obscore = spec.meta.get("obscore")
        assert obscore is not None
        assert obscore.get("calib_level") == 2
        assert obscore.get("t_xel") == 1

    def test_6dfgs_table(self, shared_datadir):
        spectra = SpectrumList.read(
            shared_datadir / SIXDFGS_TABLE_TEST_FILENAME,
            format="6dFGS-tabular obscore",
        )
        spec = spectra[0]
        obscore = spec.meta.get("obscore")
        assert obscore is not None
        assert obscore.get("calib_level") == 2
        assert obscore.get("t_xel") == 1

    def test_6dfgs_combined(self, shared_datadir):
        spectra = SpectrumList.read(
            shared_datadir / SIXDFGS_COMBINED_TEST_FILENAME,
            format="6dFGS-combined obscore",
        )
        spec = spectra[0]
        obscore = spec.meta.get("obscore")
        assert obscore is not None
        assert obscore.get("calib_level") == 2
        assert obscore.get("t_xel") == 4

    def test_gama_2dfgrs(self, shared_datadir):
        spectra = SpectrumList.read(
            shared_datadir / GAMA_2DFGRS_TEST_FILENAME,
            format="GAMA-2dFGRS obscore",
        )
        spec = spectra[0]
        obscore = spec.meta.get("obscore")
        assert obscore is not None
        assert obscore.get("calib_level") == 2
        assert obscore.get("t_xel") == 1

DC_TEST_FILENAMES = [
    "000001.fits",
    #"1704070052012671.fits",
    "1D-c0022498-344732.fits",
    #"299869.fit",
    "G12_Y2_009_044.fit",
    #"J091726.21+003424.0_a14_040423.fit",
    #"J113606.3+001155a.fit",
    #"LTF_09_1128_0555.fit",
    "MGC23320.fit",
    "OBJ0032red.fits",
    "OBJ0039red.fits",
    #"OzDES-DR2_00394.fits",
    #"OzDES-DR2_04720.fits",
    #"Spectrum-195254.fit",
    "all-c0022498-344732r_spectrum0.fits",
    "all-c0022498-344732v_spectrum0.fits",
    "all-g2359555-392832.fits",
    "spSpec-51609-0304-087.fit"
]

class TestingGuesses:
    @pytest.mark.parametrize("filename", DC_TEST_FILENAMES)
    def test_guess_all(self, shared_datadir, filename):
        formats = loaders.whatformat(shared_datadir / filename)
        spectra = SpectrumList.read(
            shared_datadir / filename
        )
        print("format for {0} is {1} and number of spectra is {2}".format(filename, formats, len(spectra)))
        for spectrum in spectra:
            print(spectrum)

MARZ_TEST_FILENAMES = [
    "marz/emlLinearSkyAirHelio.fits",
    "marz/emlLinearSkyAirHelioCMB.fits",
    "marz/emlLinearSkyAirNoHelio.fits",
    "marz/emlLinearVacuumNoHelio.fits",
    "marz/emlLogVacuumHelioCMB.fits",
    "marz/emlLogVacuumNoHelio.fits",
    "marz/quasarLinearSkyAirNoHelio.fits",
    #"marz/emlLogVacuumHelioMultiple.fits"
    #"spec-4444-55538-1000.fits"
]
MARZ_CONFIG = {
    "hdus": {
        "0": {
            "purpose": "science",
            "units": {"flux_unit": "10^-17 erg/s/cm^2/A"},
        },
        "1": {"purpose": "error_variance"},
        "2": {"purpose": "sky"},
        "3": {"purpose": "skip"},
    },
    "wcs": {
        "pixel_reference_point_keyword": "CRPIX1",
        "pixel_reference_point_value_keyword": "CRVAL1",
        "pixel_width_keyword": "CDELT1",
        "wavelength_unit": "Angstrom",
    },
    "units": {"flux_unit": "count"},
    "all_standard_units": False,
    "all_keywords": False,
    "valid_wcs": False,
}
class TestSingleSplitMarz:

    @pytest.mark.parametrize("filename", MARZ_TEST_FILENAMES)
    def test_marz(self, shared_datadir, filename):
        spectra = SpectrumList.read(
            shared_datadir / filename,
            format=loaders.SINGLE_SPLIT_LABEL,
            **MARZ_CONFIG
        )
        # Should be main spectra, with sky, and normalised
        assert len(spectra) == 2

        assert spectra[1].flux.unit == u.count

        assert spectra[0].spectral_axis.unit == u.Angstrom
        assert spectra[1].spectral_axis.unit == u.Angstrom

        assert isinstance(spectra[0].uncertainty, VarianceUncertainty)
        assert spectra[1].uncertainty is None

        #assert spectra[0].meta.get("label") is None
        assert spectra[0].meta.get("header") is not None

        #assert spectra[1].meta.get("label") is None
        assert spectra[1].meta.get("header") is not None

    def test_marz_guess_1(self, shared_datadir):
        formats = loaders.whatformat(shared_datadir / "marz/quasarLinearSkyAirNoHelio.fits")
        if len(formats) > 1:
            loaders.unregister(formats[0])
        spectra = SpectrumList.read(
            shared_datadir / "marz/quasarLinearSkyAirNoHelio.fits"
        )
        # Should be main spectra, with sky, and normalised
        assert len(spectra) == 2
    @pytest.mark.xfail(reason="data loader is not up to it yet")
    def test_marz_guess_2(self, shared_datadir):
        formats = loaders.whatformat(shared_datadir / "marz/alldata_combined_runz_x12_b02.fits")
        if len(formats) > 1:
            loaders.unregister(formats[0])
        spectra = SpectrumList.read(
            shared_datadir / "marz/alldata_combined_runz_x12_b02.fits"
        )
        # Should be main spectra, with sky, and normalised
        assert len(spectra) == 2

    def test_not_marz_guess_1(self, shared_datadir):
        formats = loaders.whatformat(shared_datadir / "OBJ0032red.fits")
        if len(formats) > 1:
            loaders.unregister(formats[0])
        spectra = SpectrumList.read(
            shared_datadir / "OBJ0032red.fits"
        )
        # Should be main spectra, with sky, and normalised
        assert len(spectra) > 0
        assert spectra[0].spectral_axis.unit == u.Angstrom
        assert spectra[0].meta.get("header") is not None
    
    def test_not_marz_guess_2(self, shared_datadir):
        formats = loaders.whatformat(shared_datadir / "J091726.21+003424.0_a14_040423.fit")
        if len(formats) > 1:
            loaders.unregister(formats[0])
        spectra = SpectrumList.read(
            shared_datadir / "J091726.21+003424.0_a14_040423.fit"
        )
        loaders.restore_registered_loaders()
        # Should be main spectra, with sky, and normalised
        assert len(spectra) > 0
        assert spectra[0].spectral_axis.unit == u.Angstrom
        assert spectra[0].meta.get("header") is not None
    #def test_marz_as_ozdes(self, shared_datadir):
    #    spectra = SpectrumList.read(
    #        shared_datadir / "marz/quasarLinearSkyAirNoHelio.fits"
    #    )
    #def test_other_marz(self, shared_datadir):
    #    spectra = SpectrumList.read(
    #        shared_datadir / "marz/alldata_combined_runz_x12_b02.fits",
    #        format="OzDES"
    #    )
