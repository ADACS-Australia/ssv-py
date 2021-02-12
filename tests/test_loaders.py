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
import ssv
import ssv.utils

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
'''
Testing the following API contract:
1.  Can read the spectra of a fits file by correctly guessing its format
2.  If the format is ambiguous then we provide a way to resolve this by "unregistering" one of the format loaders
'''
class TestSingleSplit:


    # This picks up the wcs1d-fits format that needs to be removed to ensure the Data Central drivers can be used
    def test_ozdes_guess(self, shared_datadir):
        formats = ssv.ssvloaders.whatformat(shared_datadir / OZDES_TEST_FILENAME)
        if len(formats) > 1:
            ssv.ssvloaders.unregister(formats[0])
        spectra = SpectrumList.read(shared_datadir / OZDES_TEST_FILENAME)
        ssv.ssvloaders.restore_registered_loaders()

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

    

    def test_gama_2qz_guess(self, shared_datadir):
        formats = ssv.ssvloaders.whatformat(shared_datadir / GAMA_2QZ_TEST_FILENAME)
        if len(formats) > 1:
            ssv.ssvloaders.unregister(formats[0])
        spectra = SpectrumList.read(
            shared_datadir / GAMA_2QZ_TEST_FILENAME,
        )
        ssv.ssvloaders.restore_registered_loaders()

        assert len(spectra) == 1

        assert spectra[0].flux.unit == u.count
        assert spectra[0].spectral_axis.unit == u.Angstrom

        assert isinstance(spectra[0].uncertainty, VarianceUncertainty)

        assert spectra[0].meta.get("label") is not None
        assert spectra[0].meta.get("header") is not None

    def test_2slaq_qso_guess(self, shared_datadir):
        formats = ssv.ssvloaders.whatformat(shared_datadir / GAMA_2SLAQ_QSO_TEST_FILENAME)
        if len(formats) > 1:
            ssv.ssvloaders.unregister(formats[0])
        spectra = SpectrumList.read(
            shared_datadir / GAMA_2SLAQ_QSO_TEST_FILENAME,
        )
        ssv.ssvloaders.restore_registered_loaders()

        assert len(spectra) == 1

        assert spectra[0].flux.unit == u.count
        assert spectra[0].spectral_axis.unit == u.Angstrom

        assert isinstance(spectra[0].uncertainty, VarianceUncertainty)

        assert spectra[0].meta.get("label") is not None
        assert spectra[0].meta.get("header") is not None

    #@pytest.mark.xfail(reason="Format is ambiguous")
    def test_gama_lt_guess(self, shared_datadir):
        formats = ssv.ssvloaders.whatformat(shared_datadir / GAMA_GAMA_LT_TEST_FILENAME)
        if len(formats) > 1:
            ssv.ssvloaders.unregister(formats[0])
        spectra = SpectrumList.read(
            shared_datadir / GAMA_GAMA_LT_TEST_FILENAME,
        )
        ssv.ssvloaders.restore_registered_loaders()

        assert len(spectra) == 1

        assert spectra[0].flux.unit == u.count
        assert spectra[0].spectral_axis.unit == u.Angstrom

        assert spectra[0].uncertainty is None

        assert spectra[0].meta.get("label") is not None
        assert spectra[0].meta.get("header") is not None

    #@pytest.mark.xfail(reason="Format is ambiguous")
    def test_gama_wigglez_guess(self, shared_datadir):
        formats = ssv.ssvloaders.whatformat(shared_datadir / GAMA_WIGGLEZ_TEST_FILENAME)
        if len(formats) > 1:
            ssv.ssvloaders.unregister(formats[0])
        spectra = SpectrumList.read(
            shared_datadir / GAMA_WIGGLEZ_TEST_FILENAME,
        )
        ssv.ssvloaders.restore_registered_loaders()
        
        assert len(spectra) == 1

        assert spectra[0].flux.unit == u.count
        assert spectra[0].spectral_axis.unit == u.Angstrom

        assert isinstance(spectra[0].uncertainty, VarianceUncertainty)

        assert spectra[0].meta.get("label") is not None
        assert spectra[0].meta.get("header") is not None


class TestMultilineSingle:

    def test_gama_2dfgrs_guess(self, shared_datadir):
        formats = ssv.ssvloaders.whatformat(shared_datadir / GAMA_2DFGRS_TEST_FILENAME)
        if len(formats) > 1:
            for i in range(len(formats)-1):
                ssv.ssvloaders.unregister(formats[i])
        spectra = SpectrumList.read(
            shared_datadir / GAMA_2DFGRS_TEST_FILENAME,
        )
        ssv.ssvloaders.restore_registered_loaders()

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
        formats = ssv.ssvloaders.whatformat(shared_datadir / filename)
        spectra = SpectrumList.read(
            shared_datadir / filename
        )
        print("format for {0} is {1} and number of spectra is {2}".format(filename, formats, len(spectra)))
        for spectrum in spectra:
            print(spectrum)

    @pytest.mark.parametrize("filename", DC_TEST_FILENAMES)
    def test_guess_all_utils(self, shared_datadir, filename):
        spectra = ssv.utils.read_spectra_file(
            shared_datadir / filename
        )
        assert spectra is not None

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
class TestSingleSplitMarz:

    @pytest.mark.parametrize("filename", MARZ_TEST_FILENAMES)
    def test_marz(self, shared_datadir, filename):
        formats = ssv.ssvloaders.whatformat(shared_datadir / filename)
        if len(formats) > 1:
            ssv.ssvloaders.unregister(formats[0])
        spectra = SpectrumList.read(
            shared_datadir / filename
        )
        ssv.ssvloaders.restore_registered_loaders()
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

    @pytest.mark.parametrize("filename", MARZ_TEST_FILENAMES)
    def test_marz_spectra(self, shared_datadir, filename):
        spectra = ssv.utils.read_spectra_file(
            shared_datadir / filename
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
        formats = ssv.ssvloaders.whatformat(shared_datadir / "marz/quasarLinearSkyAirNoHelio.fits")
        if len(formats) > 1:
            ssv.ssvloaders.unregister(formats[0])
        spectra = SpectrumList.read(
            shared_datadir / "marz/quasarLinearSkyAirNoHelio.fits"
        )
        ssv.ssvloaders.restore_registered_loaders()

        # Should be main spectra, with sky, and normalised
        assert len(spectra) == 2
    @pytest.mark.xfail(reason="data loader is not up to it yet")
    def test_marz_guess_2(self, shared_datadir):
        formats = ssv.ssvloaders.whatformat(shared_datadir / "marz/alldata_combined_runz_x12_b02.fits")
        if len(formats) > 1:
            ssv.ssvloaders.unregister(formats[0])
        spectra = SpectrumList.read(
            shared_datadir / "marz/alldata_combined_runz_x12_b02.fits"
        )
        # Should be main spectra, with sky, and normalised
        assert len(spectra) == 2

    def test_not_marz_guess_1(self, shared_datadir):
        formats = ssv.ssvloaders.whatformat(shared_datadir / "OBJ0032red.fits")
        if len(formats) > 1:
            ssv.ssvloaders.unregister(formats[0])
        spectra = SpectrumList.read(
            shared_datadir / "OBJ0032red.fits"
        )
        # Should be main spectra, with sky, and normalised
        assert len(spectra) > 0
        assert spectra[0].spectral_axis.unit == u.Angstrom
        assert spectra[0].meta.get("header") is not None
    
    def test_not_marz_guess_2(self, shared_datadir):
        formats = ssv.ssvloaders.whatformat(shared_datadir / "J091726.21+003424.0_a14_040423.fit")
        if len(formats) > 1:
            ssv.ssvloaders.unregister(formats[0])
        spectra = SpectrumList.read(
            shared_datadir / "J091726.21+003424.0_a14_040423.fit"
        )
        ssv.ssvloaders.restore_registered_loaders()
        # Should be main spectra, with sky, and normalised
        assert len(spectra) > 0
        assert spectra[0].spectral_axis.unit == u.Angstrom
        assert spectra[0].meta.get("header") is not None

class TestLoaderRegistryWorkaround:
    def test_unregister_and_register(self, shared_datadir):
        ssv.ssvloaders.restore_registered_loaders()
        formats = ssv.ssvloaders.whatformat(shared_datadir / GAMA_2SLAQ_QSO_TEST_FILENAME)
        assert len(formats) == 2

        ssv.ssvloaders.unregister(formats[0])
        formats = ssv.ssvloaders.whatformat(shared_datadir / GAMA_2SLAQ_QSO_TEST_FILENAME)
        assert len(formats) == 1

        ssv.ssvloaders.restore_registered_loaders()
        formats = ssv.ssvloaders.whatformat(shared_datadir / GAMA_2SLAQ_QSO_TEST_FILENAME)
        assert len(formats) == 2
