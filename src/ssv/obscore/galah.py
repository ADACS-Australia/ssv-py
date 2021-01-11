import re

from specutils import SpectrumList
from specutils.io.registers import data_loader

from .common import get_times
from ..loaders import FITS_FILE_EXTS, no_auto_identify


@data_loader(
    label="GALAH obscore", extensions=FITS_FILE_EXTS, dtype=SpectrumList,
    identifier=no_auto_identify,
)
def galah_obscore_loader(fname):
    spectra = SpectrumList.read(fname, format="GALAH")
    spec = spectra[0]
    spec.meta["obscore"] = {}
    obscore = spec.meta["obscore"]
    hdr = spec.meta["header"]

    # band_name = B, V, R or I
    if hdr["SPECTID"] == "BL":
        obscore["band_name"] = "B"
    if hdr["SPECTID"] == "RR":
        obscore["band_name"] = "I"
    if hdr["SPECTID"] == "RD":
        obscore["band_name"] = "R"
    if hdr["SPECTID"] == "GN":
        obscore["band_name"] = "V"
    t1, t2 = get_times(hdr, duration_kw="EXPOSED")
    if t1 is not None:
        obscore["t_min"] = t1.to_value('mjd',subfmt='float')
    if t2 is not None:
        obscore["t_max"] = t2.to_value('mjd',subfmt='float')

    obscore["s_ra"] = hdr["MEANRA"]
    obscore["s_dec"] = hdr["MEANDEC"]
    obscore["s_fov"] = 2.1 / 3600
    obscore["obs_collection"] = "galah_dr3"
    obscore["facility_name"] = "AAT"
    # obscore["dataproduct_type"] = "spectrum"
    if hdr["NCOMBINE"] > 1:
        obscore["dataproduct_subtype"] = "combined"
    else:
        obscore["dataproduct_subtype"] = "science"
    obscore["calib_level"] = 2

    obscore["t_exptime"] = hdr["TOTALEXP"]
    obscore["t_resolution"] = hdr["EXPOSED"]
    obscore["t_xel"] = hdr["NCOMBINE"]

    nspecpix = len(spec.spectral_axis)
    obscore["em_xel"] = nspecpix
    obscore["em_ucd"] = "em.wl"
    obscore["em_unit"] = "angstrom"
    obscore["s_xel1"] = nspecpix
    obscore["s_xel2"] = 1
    obscore["em_min"] = spec.spectral_axis[0].meter
    obscore["em_max"] = spec.spectral_axis[-1].meter

    # some info in aao_itc/includes/itc_forms.inc , but not resolving powers -
    # in another program?
    obscore["em_res_power"] = 26000
    # get these from the function....
    # if(cen_res == None):
    #    (cen_res,cen_rp,cen_rp_min,cen_rp_max) = CalcAAOmegaResolutions(hdr)
    # if(cen_res != None):
    #    obscore['em_res_power'] = cen_rp
    #    obscore['em_res_power_min'] = cen_rp_min
    #    obscore['em_res_power_max'] = cen_rp_max
    #    obscore['em_resolution'] = cen_res*1E-10

    obscore["o_ucd"] = "phot.count"
    # spectra are calibrated in flux: ROW1 hdr comment: Flux-calibrated
    # spectrum in 10^-17 erg/s/cm^2/A
    # obscore['o_ucd'] = 'phot.flux'
    # obscore['o_unit'] = '1.0E-17 erg/s/cm^2/A'
    # obscore['o_calib_status'] = 'absolute'
    # or perhaps 'relative' since these are fibre spectra?
    obscore["instrument_name"] = "2dF-HERMES"
    obscore["proposal_id"] = re.sub("\\s+", "", hdr["AAOPRGID"])
    obscore["em_calib_status"] = "calibrated"
    # Should see if we can extract this from somewhere else than the file name,
    # as we may not always have access to it
    objid = re.sub("[0-9].fits", "", re.sub(".*/", "", str(fname)))
    obscore["target_name"] = objid
    obscore["obs_id"] = objid
    # if('OBJCOM' in hdr):
    #    obscore['alt_target_name'] = hdr['OBJCOM']
    # if('WAVRESOL' in hdr):
    #    obscore['em_resolution'] = hdr['WAVRESOL']*1e-10

    # if('SN' in hdr):
    #    obscore['em_snr'] = hdr['SN']

    # alternative name: OBJECT
    # if('Z' in hdr):
    #    obscore['redshift'] = hdr['Z']
    return spectra
