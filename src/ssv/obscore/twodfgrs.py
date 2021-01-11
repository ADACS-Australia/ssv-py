import numpy as np

from astropy.coordinates import FK4
from astropy.coordinates import SkyCoord

from specutils import SpectrumList
from specutils.io.registers import data_loader

from .common import get_times
from ..loaders import FITS_FILE_EXTS, no_auto_identify


@data_loader(
    label="2dFGRS obscore", extensions=FITS_FILE_EXTS, dtype=SpectrumList,
    identifier=no_auto_identify,
)
def obscore_loader_2dfgrs(fname):
    spectra = SpectrumList.read(fname, format="2dFGRS")

    # primary_header is for the image in first extension...
    # There may be multiple observations of one spectrum; taken on multiple
    # dates
    for spec in spectra:
        spec.meta["obscore"] = {}
        obscore = spec.meta["obscore"]
        hdr = spec.meta["header"]
        primary_hdr = spec.meta["primary_header"]

        c = SkyCoord(
            ra=hdr["OBSRA"] * 180.0 / np.pi,
            dec=hdr["OBSDEC"] * 180.0 / np.pi,
            unit="deg",
            frame=FK4,
            equinox="B1950",
        )
        c2 = c.transform_to("icrs")

        t1, t2 = get_times(hdr)
        if t1 is not None:
            obscore["t_min"] = t1.to_value('mjd',subfmt='float')
        if t2 is not None:
            obscore["t_max"] = t2.to_value('mjd',subfmt='float')

        # These are in RADIANS in the header!
        # OBJRA and OBJDEC are the most accurate.
        obscore["s_ra"] = c2.ra.deg
        obscore["s_dec"] = c2.dec.deg
        obscore["s_fov"] = 2.1 / 3600
        # obscore["obs_collection"] = "2dfgrs_fdr"
        # obscore["facility_name"] = "AAT"
        # obscore["dataproduct_type"] = "spectrum"
        obscore["dataproduct_subtype"] = "science"
        obscore["calib_level"] = 2

        if "REFEXP" in hdr:
            obscore["t_exptime"] = hdr["REFEXP"]

        nspecpix = len(spec.spectral_axis)
        obscore["em_xel"] = nspecpix
        obscore["em_ucd"] = "em.wl"
        obscore["em_unit"] = "angstrom"
        obscore["s_xel1"] = nspecpix
        obscore["s_xel2"] = 1
        obscore["t_xel"] = 1
        obscore["em_min"] = spec.spectral_axis[0].meter
        obscore["em_max"] = spec.spectral_axis[-1].meter
        obscore["em_res_power"] = 644
        obscore["em_res_power_min"] = 400
        obscore["em_res_power_max"] = 889
        obscore["em_resolution"] = 9.0 * 1e-10

        obscore["o_ucd"] = "phot.count"
        # spectra are calibrated in flux: ROW1 hdr comment: Flux-calibrated
        # spectrum in 10^-17 erg/s/cm^2/A
        # obscore['o_ucd'] = 'phot.flux'
        # obscore['o_unit'] = '1.0E-17 erg/s/cm^2/A'
        # obscore['o_calib_status'] = 'absolute'
        # or perhaps 'relative' since these are fibre spectra?
        obscore["instrument_name"] = "2dF"
        obscore["em_calib_status"] = "calibrated"

        if "NAME" in primary_hdr:
            obscore["target_name"] = primary_hdr["NAME"]
        if "OBSNAME" in hdr:
            obscore["alt_target_name"] = hdr["OBSNAME"]
        # if('WAVRESOL' in hdr):
        #    obscore['em_resolution'] = hdr['WAVRESOL']*1e-10

        if "SNR" in hdr:
            obscore["em_snr"] = hdr["SNR"]

        # alternative name: OBJECT
        if "SPECID" in hdr:
            obscore["obs_id"] = hdr["SPECID"]
        if "Z" in hdr:
            obscore["redshift"] = hdr["Z"]
    return spectra
