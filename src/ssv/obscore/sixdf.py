import re

from astropy.time import Time, TimeDelta

from specutils import SpectrumList
from specutils.io.registers import data_loader

from ..loaders import no_auto_identify


@data_loader(
    label="6dFGS-tabular obscore", dtype=SpectrumList,
    identifier=no_auto_identify,
)
def sixdf_1d_obscore_loader(fname):
    spectra = SpectrumList.read(fname, format="6dFGS-tabular")
    spec = spectra[0]

    spec.meta["obscore"] = {}
    obscore = spec.meta["obscore"]
    hdr = spec.meta["header"]
    # No date or time info in 6dF spectra

    obscore["s_ra"] = hdr["OBSRA"]
    obscore["s_dec"] = hdr["OBSDEC"]
    obscore["s_fov"] = 6.7 / 3600
    obscore["facility_name"] = "UKST"
    # obscore["dataproduct_type"] = "spectrum"
    obscore["dataproduct_subtype"] = "science"
    obscore["calib_level"] = 2

    # No exposure time

    nspecpix = len(spec.spectral_axis)
    obscore["em_xel"] = nspecpix
    obscore["em_ucd"] = "em.wl"
    obscore["em_unit"] = "angstrom"
    obscore["s_xel1"] = nspecpix
    obscore["s_xel2"] = 1
    obscore["t_xel"] = 1
    obscore["em_min"] = spec.spectral_axis[0].meter
    obscore["em_max"] = spec.spectral_axis[-1].meter
    obscore["em_res_power"] = 1000
    obscore["em_resolution"] = 8.5 * 1e-10

    obscore["o_ucd"] = "phot.flux"
    obscore["o_unit"] = "counts/s"
    obscore["o_calib_status"] = "absolute"
    obscore["instrument_name"] = "6dF"
    obscore["em_calib_status"] = "calibrated"

    obscore["target_name"] = hdr.get("NAME_V")
    obscore["alt_target_name"] = hdr.get("TARGET")
    obscore["redshift"] = hdr.get("Z")

    return spectra


@data_loader(
    label="6dFGS-combined obscore", dtype=SpectrumList,
    identifier=no_auto_identify,
)
def sixdf_combined_obscore_loader(fname):
    spectra = SpectrumList.read(fname, format="6dFGS-combined")

    # number of spectra in the SpectrumList
    ridx = None
    vidx = None
    vridx = None
    for idx, spec in enumerate(spectra):
        hdr = spec.meta["header"]
        band = re.sub("SPECTRUM ", "", hdr["EXTNAME"])
        if band == "V":
            vidx = idx
        if band == "R":
            ridx = idx
        if band == "VR":
            vridx = idx

    # calculate date/time info here for convenience
    # as they depend on each other a little
    vhdr = spectra[vidx].meta["header"]
    rhdr = spectra[ridx].meta["header"]

    rexp = rhdr["EXP_R"]
    rnexp = rhdr["NCOMB_R"]
    vexp = vhdr["EXP_V"]
    vnexp = vhdr["NCOMB_V"]

    rdate = re.sub("/", "-", rhdr["UTDATE_R"])
    rstart = rhdr["UTSTRT_R"]
    rt1 = Time(rdate + "T" + rstart)
    rt2 = rt1 + TimeDelta(rnexp * rexp, format="sec")

    vdate = re.sub("/", "-", vhdr["UTDATE_V"])
    vstart = vhdr["UTSTRT_V"]
    vt1 = Time(vdate + "T" + vstart)
    vt2 = vt1 + TimeDelta(vnexp * vexp, format="sec")

    vrt1 = min(rt1, vt1)
    vrt2 = max(rt2, vt2)

    for spec in spectra:
        spec.meta["obscore"] = {}
        obscore = spec.meta["obscore"]
        hdr = spec.meta["header"]
        band = re.sub("SPECTRUM ", "", hdr["EXTNAME"])

        # no date or time info in 6dF spectra
        if band == "V":
            obscore["t_min"] = vt1.to_value('mjd',subfmt='float')
            obscore["t_max"] = vt2.to_value('mjd',subfmt='float')
            obscore["t_resolution"] = vexp
            obscore["t_xel"] = vnexp
            obscore["t_exptime"] = vexp * vnexp
            obscore["calib_level"] = 2
            obscore["obs_id"] = hdr["OBSID_V"]
            obscore["band_name"] = "V"
            if vnexp > 1:
                obscore["dataproduct_subtype"] = "combined"
            else:
                obscore["dataproduct_subtype"] = "science"

        if band == "R":
            obscore["t_min"] = rt1.to_value('mjd',subfmt='float')
            obscore["t_max"] = rt2.to_value('mjd',subfmt='float')
            obscore["t_xel"] = rnexp
            obscore["t_resolution"] = rexp
            obscore["t_exptime"] = rexp * rnexp
            obscore["obs_id"] = hdr["OBSID_R"]
            obscore["calib_level"] = 2
            obscore["band_name"] = "R"
            if rnexp > 1:
                obscore["dataproduct_subtype"] = "combined"
            else:
                obscore["dataproduct_subtype"] = "science"
        if band == "VR":
            obscore["t_min"] = vrt1.to_value('mjd',subfmt='float')
            obscore["t_max"] = vrt2.to_value('mjd',subfmt='float')
            obscore["t_resolution"] = min(rexp, vexp)
            obscore["t_exptime"] = rexp * rnexp + vexp * vnexp
            obscore["t_xel"] = rnexp + vnexp
            obscore["obs_id"] = vhdr["OBSID_V"] + "-" + rhdr["OBSID_R"]
            obscore["calib_level"] = 3
            obscore["band_name"] = "VR"
            obscore["dataproduct_subtype"] = "combined"

        obscore["s_ra"] = hdr["OBSRA"]
        obscore["s_dec"] = hdr["OBSDEC"]
        obscore["s_fov"] = 6.7 / 3600
        obscore["facility_name"] = "UKST"
        # obscore["dataproduct_type"] = "spectrum"

        nspecpix = len(spec.spectral_axis)
        obscore["em_xel"] = nspecpix
        obscore["em_ucd"] = "em.wl"
        obscore["em_unit"] = "angstrom"
        obscore["s_xel1"] = nspecpix
        obscore["s_xel2"] = 1
        obscore["em_min"] = spec.spectral_axis[0].meter
        obscore["em_max"] = spec.spectral_axis[-1].meter
        # obscore['em_res_power'] = 1000
        # obscore['em_res_power_min'] = 1050
        # obscore['em_res_power_max'] = 1681
        # obscore['em_resolution'] = 8.5*1E-10

        obscore["o_ucd"] = "phot.count"
        obscore["instrument_name"] = "6dF"
        obscore["em_calib_status"] = "calibrated"

        if "NAME_V" in vhdr:
            obscore["target_name"] = vhdr["NAME_V"]
        if "TARGET" in hdr and hdr["TARGET"] != "not__in__config":
            obscore["alt_target_name"] = hdr["TARGET"]
        # if('WAVRESOL' in hdr):
        #    obscore['em_resolution'] = hdr['WAVRESOL']*1e-10

        # if('SN' in hdr):
        # obscore['em_snr'] = 10

        if "Z_HELIO" in hdr:
            obscore["redshift"] = hdr["Z_HELIO"]

    return spectra
