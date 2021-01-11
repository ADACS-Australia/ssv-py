from specutils import SpectrumList
from specutils.io.registers import data_loader

from .common import get_times
from ..loaders import FITS_FILE_EXTS, no_auto_identify


@data_loader(
    label="OzDES obscore", extensions=FITS_FILE_EXTS, dtype=SpectrumList,
    identifier=no_auto_identify,
)
def ozdes_obscore_loader(fname):
    spectra = SpectrumList.read(fname, format="OzDES")

    # a list to store the times from the spectra that have the info in the
    # headers
    times = []
    # used to sum up the exposure time of all spectra
    combined_time = 0
    # index of the combined spectrum; should always be zero, but we assign it
    # when purpose == "combined"
    combined_idx = 0
    # number of individual spectra
    nepochs = 0
    mid_mjd = []
    redshift = None
    for idx, spec in enumerate(spectra):
        # dict to store all the obscore params
        # easier to pass one obscore dict to other functions, than remembering
        # which obscore params were assigned
        spec.meta["obscore"] = {}
        obscore = spec.meta["obscore"]

        hdr = spec.meta["header"]
        t1, t2 = get_times(hdr, duration_kw="EXPOSED")
        if t1 is not None:
            times.append(t1)
            obscore["t_min"] = t1.to_value('mjd',subfmt='float')
        if t2 is not None:
            times.append(t2)
            obscore["t_max"] = t2.to_value('mjd',subfmt='float')
        if t1 is not None and t2 is not None:
            mid_mjd.append(t1.to_value('mjd',subfmt='float') + (t2.to_value('mjd',subfmt='float') - t1.to_value('mjd',subfmt='float')) * 0.5)
        obscore["s_ra"] = hdr["RA"]
        obscore["s_dec"] = hdr["DEC"]
        obscore["s_fov"] = 2.1 / 3600

        # obscore["obs_collection"] = "ozdes_dr2"
        # obscore['facility_name'] = None
        # 'N/A'; does not seem relevant here

        if "EXPOSED" in hdr:
            exptime = hdr["EXPOSED"]
            combined_time = combined_time + exptime
            obscore["t_exptime"] = exptime
            # obscore['t_resolution'] = exptime
        # else:
        #    obscore['t_resolution'] = None

        nspecpix = len(spec.spectral_axis)
        obscore["em_xel"] = nspecpix
        obscore["em_ucd"] = "em.wl"
        obscore["em_unit"] = "angstrom"
        obscore["s_xel1"] = nspecpix
        obscore["s_xel2"] = 1
        obscore["em_min"] = spec.spectral_axis[0].meter
        obscore["em_max"] = spec.spectral_axis[-1].meter
        obscore["em_res_power"] = 1361
        obscore["em_res_power_min"] = 1050
        obscore["em_res_power_max"] = 1681
        obscore["em_resolution"] = 5.318 * 1e-10

        obscore["instrument_name"] = "2dF-AAOmega"
        obscore["facility_name"] = "AAT"
        obscore["em_calib_status"] = "calibrated"

        obscore["t_xel"] = 1

        # obscore["dataproduct_type"] = "spectrum"
        if spec.meta["purpose"] == "combined":
            obscore["dataproduct_subtype"] = "combined"
            # since it is a combined spectrum, it is slightly higher
            # calib_level than ordinary spectra (2)
            obscore["calib_level"] = 3
            if "SOURCE" in hdr:
                obscore["target_name"] = hdr["SOURCE"]
                obscore["obs_id"] = hdr["SOURCE"]
            if "Z" in hdr:
                redshift = hdr["Z"]
                obscore["redshift"] = redshift

            combined_idx = idx
        else:
            obscore["dataproduct_subtype"] = "science"
            obscore["calib_level"] = 2
            nepochs = nepochs + 1
            if "ORIGTARG" in hdr:
                obscore["target_name"] = hdr["ORIGTARG"]
                obscore["obs_id"] = hdr["ORIGTARG"]

    combined_obscore = spectra[combined_idx].meta["obscore"]
    if len(mid_mjd) > 1:
        delta = mid_mjd[1] - mid_mjd[0]
        for idx in range(2, len(mid_mjd)):
            dm = mid_mjd[idx] - mid_mjd[idx - 1]
            if dm < delta:
                delta = dm
        combined_obscore["t_resolution"] = delta * 24.0 * 3600

    # set some values for the combined frame...
    combined_obscore["t_min"] = min(times).to_value('mjd',subfmt='float')
    combined_obscore["t_max"] = max(times).to_value('mjd',subfmt='float')
    combined_obscore["t_exptime"] = combined_time
    combined_obscore["t_xel"] = nepochs

    for spec in spectra:
        if redshift is not None:
            obscore["redshift"] = redshift

    return spectra
