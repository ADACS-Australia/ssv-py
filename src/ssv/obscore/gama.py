from astropy.time import Time, TimeDelta

from specutils import SpectrumList, Spectrum1D
from specutils.io.registers import data_loader

from ..loaders import FITS_FILE_EXTS, no_auto_identify
from .common import get_times


def gama_sdss_time_to_mjd(sdss_time):
    """
    SDSS stores time in seconds since 0 MJD using TAI as the scale, this
    converts the value to a more normal mjd
    """
    if sdss_time is None:
        return None
    base_time = Time(0.0, format="mjd")
    return base_time + TimeDelta(sdss_time, format="sec", scale="tai")


@data_loader(
    label="GAMA-SDSS obscore", extensions=FITS_FILE_EXTS, dtype=Spectrum1D,
    identifier=no_auto_identify,
)
def obscore_loader_gama_sdss(fname):
    """
    This is a specutils loader which adds obscore information to a SDSS spectra
    in the GAMA survey
    """
    spec = Spectrum1D.read(fname, format="SDSS-I/II spSpec")
    spec.meta["obscore"] = {}
    obscore = spec.meta["obscore"]

    hdr = spec.meta["header"]

    obscore["t_min"] = gama_sdss_time_to_mjd(hdr.get("TAI-BEG")).to_value('mjd',subfmt='float')
    obscore["t_max"] = gama_sdss_time_to_mjd(hdr.get("TAI-END")).to_value('mjd',subfmt='float')
    obscore["s_ra"] = hdr["RA"]
    obscore["s_dec"] = hdr["DEC"]
    obscore["s_seeing"] = hdr.get("SEEING50")
    # obscore["obs_collection"] = "gama_dr2"
    # obscore["facility_name"] = "sdss"
    obscore["dataproduct_subtype"] = "science"
    obscore["calib_level"] = 2

    obscore["t_exptime"] = hdr.get("EXPTIME")

    nspecpix = len(spec.spectral_axis)
    obscore["em_xel"] = nspecpix
    obscore["em_ucd"] = "em.wl"
    obscore["em_unit"] = spec.spectral_axis.unit
    obscore["s_xel1"] = nspecpix
    obscore["s_xel2"] = 1
    obscore["t_xel"] = 1
    obscore["em_min"] = spec.spectral_axis[0].meter
    obscore["em_max"] = spec.spectral_axis[-1].meter
    obscore["em_res_power"] = 2000
    # http://classic.sdss.org/dr7/instruments/spectrographs/index.html
    obscore["em_res_power_min"] = 1850
    obscore["em_res_power_max"] = 2200
    # lambda_cen=(3800+9200)/2=6500 A; Delta Lambda = 6500/2000 = 3.25
    obscore["em_resolution"] = 3.25 * 1e-10

    # obscore['o_ucd'] = 'phot.count'
    # spectra are calibrated in flux:
    # BUNIT   = '1.0E-17 erg/cm/s/Ang' / units
    obscore["o_ucd"] = "phot.flux"
    obscore["o_unit"] = "1.0E-17 erg/s/cm^2/A"
    obscore["o_calib_status"] = "absolute"
    # or perhaps 'relative' since these are fibre spectra?

    obscore["instrument_name"] = "SDSS"
    # 3 arcsec diameter fibres
    obscore["s_fov"] = 3.0 / 3600
    obscore["em_calib_status"] = "calibrated"

    obscore["target_name"] = hdr.get("GAMANAME")
    # if('OBJECT' in hdr):
    #    obscore['alt_target_name'] = hdr['OBJECT']
    # if('WAVRESOL' in hdr):
    #    obscore['em_resolution'] = hdr['WAVRESOL']*1e-10

    # SN_G and SN_I also available, choosing R-band as good compromise/average
    obscore["em_snr"] = hdr.get("SN_R")

    # alternative name: OBJECT
    obscore["obs_id"] = hdr.get("SPECID")
    obscore["redshift"] = hdr.get("Z")

    return spec


@data_loader(
    label="GAMA-WiggleZ obscore", dtype=SpectrumList,
    identifier=no_auto_identify,
)
def wigglez_obscore_loader(fname):
    spectra = SpectrumList.read(fname, format="GAMA-WiggleZ")
    for spec in spectra:
        spec.meta["obscore"] = {}
        obscore = spec.meta["obscore"]
        if spec.meta["purpose"] == "reduced":

            hdr = spec.meta["header"]
            # NOTE: It's unclear whether DATE in header (Image generation time)
            # corresponds to date of observation or reduction
            # Missing start and end times of exposures; Missing exposure time
            t1 = hdr.get("DATE")
            if t1:
                tval = Time(t1, format="isot", scale="utc")
                obscore["t_min"] = tval.to_value('mjd',subfmt='float')
            # if t2 is not None:
            #    obscore["t_max"] = t2.to_value('mjd',subfmt='float')

            obscore["s_ra"] = hdr["RA"]
            obscore["s_dec"] = hdr["DEC"]
            obscore["s_fov"] = 2.1 / 3600
            # obscore["obs_collection"] = "gama_dr2"
            # obscore["facility_name"] = "wigglez"

            obscore["dataproduct_subtype"] = "science"
            obscore["calib_level"] = 2
            # Missing exptime from Wigglez spectra! :-(
            # if('T_EXP' in hdr):
            obscore["t_exptime"] = 3600

            nspecpix = len(spec.spectral_axis)
            obscore["em_xel"] = nspecpix
            obscore["em_ucd"] = "em.wl"
            obscore["em_unit"] = "angstrom"
            obscore["s_xel1"] = nspecpix
            obscore["s_xel2"] = 1
            obscore["t_xel"] = 1
            obscore["em_min"] = spec.spectral_axis[0].meter
            obscore["em_max"] = (
                spec.spectral_axis[nspecpix - 1].meter
            )

            obscore["o_ucd"] = "phot.count"
            # spectra are calibrated in flux: ROW1 hdr comment: Flux-calibrated
            # spectrum in 10^-17 erg/s/cm^2/A
            # obscore['o_ucd'] = 'phot.flux'
            # or perhaps 'relative' since these are fibre spectra?
            # obscore['o_calib_status'] = 'absolute'
            obscore["instrument_name"] = "2dF-AAOmega"
            obscore["em_calib_status"] = "calibrated"

            if "GAMANAME" in hdr:
                obscore["target_name"] = hdr["GAMANAME"]
            if "OBJECT" in hdr:
                obscore["alt_target_name"] = hdr["OBJECT"]
            # if('WAVRESOL' in hdr):
            #    obscore['em_resolution'] = hdr['WAVRESOL']*1e-10

            # if('SN' in hdr):
            #    obscore['em_snr'] = hdr['SN']

            # alternative name: OBJECT
            if "SPECID" in hdr:
                obscore["obs_id"] = hdr["SPECID"]
            if "Z" in hdr:
                obscore["redshift"] = hdr["Z"]

    return spectra


@data_loader(
    label="GAMA-6dFGS obscore", extensions=FITS_FILE_EXTS, dtype=SpectrumList,
    identifier=no_auto_identify,
)
def gama_6dfgs_obscore_loader(fname):
    spec = Spectrum1D.read(fname, format="6dFGS-split")
    obscore = {}
    hdr = spec.meta["header"]
    # no date or time info in 6dF spectra
    #   (t1,t2) = GetTimes_GAMA_GAMA(hdr)
    #   if(t1 != None):
    #       obscore['t_min'] = t1.to_value('mjd',subfmt='float')
    #   if(t2 != None):
    #       obscore['t_max'] = t2.to_value('mjd',subfmt='float')

    obscore["s_ra"] = hdr["RA"]
    obscore["s_dec"] = hdr["DEC"]
    obscore["s_fov"] = 6.7 / 3600
    # obscore["obs_collection"] = "gama_dr2"
    # obscore["facility_name"] = "6dfgs"
    obscore["dataproduct_subtype"] = "science"
    obscore["calib_level"] = 2

    # No exposure time
    # if('T_EXP' in hdr):
    #    obscore['t_exptime'] = hdr['T_EXP']
    # exposure times are 3x10 min (red arm) and 3x20 min (visual arm); Jones et
    # al. 2004
    obscore["t_exptime"] = 3600

    nspecpix = len(spec.spectral_axis)
    obscore["em_xel"] = nspecpix
    obscore["em_ucd"] = "em.wl"
    obscore["em_unit"] = "angstrom"
    obscore["s_xel1"] = nspecpix
    obscore["s_xel2"] = 1
    obscore["t_xel"] = 1
    obscore["em_min"] = spec.spectral_axis[0].meter
    obscore["em_max"] = spec.spectral_axis[nspecpix - 1].meter
    obscore["em_res_power"] = 1000
    # obscore['em_res_power_min'] = 1050
    # obscore['em_res_power_max'] = 1681
    obscore["em_resolution"] = 8.5 * 1e-10

    obscore["o_ucd"] = "phot.count"
    # spectra are calibrated in flux: ROW1 hdr comment: Flux-calibrated
    # spectrum in 10^-17 erg/s/cm^2/A
    # obscore['o_ucd'] = 'phot.flux'
    # obscore['o_unit'] = '1.0E-17 erg/s/cm^2/A'
    # or perhaps 'relative' since these are fibre spectra?
    # obscore['o_calib_status'] = 'absolute'
    obscore["instrument_name"] = "6dF"
    obscore["em_calib_status"] = "calibrated"

    if "GAMANAME" in hdr:
        obscore["target_name"] = hdr["GAMANAME"]
    if "TARGET" in hdr:
        obscore["alt_target_name"] = hdr["TARGET"]
    # if('WAVRESOL' in hdr):
    #    obscore['em_resolution'] = hdr['WAVRESOL']*1e-10

    # if('SN' in hdr):
    obscore["em_snr"] = 10

    # alternative name: OBJECT
    if "SPECID" in hdr:
        obscore["obs_id"] = hdr["SPECID"]
    if "Z" in hdr:
        obscore["redshift"] = hdr["Z"]

    spec.meta["obscore"] = obscore
    return spec


@data_loader(
    label="GAMA-2dFGRS obscore", extensions=FITS_FILE_EXTS, dtype=SpectrumList,
    identifier=no_auto_identify,
)
def gama_2dfgrs_obscore_loader(fname):
    spectra = SpectrumList.read(fname, format="GAMA-2dFGRS")
    # number of spectra in the SpectrumList
    for spec in spectra:
        spec.meta["obscore"] = {}
        obscore = spec.meta["obscore"]
        if spec.meta["purpose"] == "reduced":

            hdr = spec.meta["header"]
            t1, t2 = get_times(hdr)
            if t1 is not None:
                obscore["t_min"] = t1.to_value('mjd',subfmt='float')
            if t2 is not None:
                obscore["t_max"] = t2.to_value('mjd',subfmt='float')

            obscore["s_ra"] = hdr["RA"]
            obscore["s_dec"] = hdr["DEC"]
            obscore["s_fov"] = 2.1 / 3600
            if "SEEING" in hdr:
                obscore["s_seeing"] = hdr["SEEING"]
            # obscore["obs_collection"] = "gama_dr2"
            # obscore["facility_name"] = "2dfgrs"
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
            obscore["em_max"] = (
                spec.spectral_axis[nspecpix - 1].meter
            )
            obscore["em_res_power"] = 644
            obscore["em_res_power_min"] = 400
            obscore["em_res_power_max"] = 889
            obscore["em_resolution"] = 9.0 * 1e-10

            obscore["o_ucd"] = "phot.count"
            # spectra are calibrated in flux: ROW1 hdr comment: Flux-calibrated
            # spectrum in 10^-17 erg/s/cm^2/A
            # obscore['o_ucd'] = 'phot.flux'
            # obscore['o_unit'] = '1.0E-17 erg/s/cm^2/A'
            # or perhaps 'relative' since these are fibre spectra?
            # obscore['o_calib_status'] = 'absolute'
            obscore["instrument_name"] = "2dF"
            obscore["em_calib_status"] = "calibrated"

            if "GAMANAME" in hdr:
                obscore["target_name"] = hdr["GAMANAME"]
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


@data_loader(
    label="GAMA obscore", extensions=FITS_FILE_EXTS, dtype=SpectrumList,
    identifier=no_auto_identify,
)
def gama_obscore_loader(fname):
    spectra = SpectrumList.read(fname, format="GAMA")
    for spec in spectra:
        spec.meta["obscore"] = {}
        obscore = spec.meta["obscore"]
        if spec.meta["purpose"] == "reduced":

            hdr = spec.meta["header"]
            t1, t2 = get_times(
                hdr, date_kw="DATE-OBS", duration_kw="T_EXP"
            )
            if t1 is not None:
                obscore["t_min"] = t1.to_value('mjd',subfmt='float')
            if t2 is not None:
                obscore["t_max"] = t2.to_value('mjd',subfmt='float')

            obscore["s_ra"] = hdr["RA"]
            obscore["s_dec"] = hdr["DEC"]
            obscore["s_fov"] = 2.1 / 3600
            # obscore["obs_collection"] = "gama_dr2"
            # obscore["facility_name"] = "gama"
            obscore["dataproduct_subtype"] = "science"
            obscore["calib_level"] = 2

            if "T_EXP" in hdr:
                obscore["t_exptime"] = hdr["T_EXP"]

            nspecpix = len(spec.spectral_axis)
            obscore["em_xel"] = nspecpix
            obscore["em_ucd"] = "em.wl"
            obscore["em_unit"] = "angstrom"
            obscore["s_xel1"] = nspecpix
            obscore["s_xel2"] = 1
            obscore["t_xel"] = 1
            obscore["em_min"] = spec.spectral_axis[0].meter
            obscore["em_max"] = (
                spec.spectral_axis[nspecpix - 1].meter
            )
            obscore["em_res_power"] = 1361
            obscore["em_res_power_min"] = 1050
            obscore["em_res_power_max"] = 1681
            obscore["em_resolution"] = 5.318 * 1e-10

            # obscore['o_ucd'] = 'phot.count'
            # spectra are calibrated in flux: ROW1 hdr comment: Flux-calibrated
            # spectrum in 10^-17 erg/s/cm^2/A
            obscore["o_ucd"] = "phot.flux"
            obscore["o_unit"] = "1.0E-17 erg/s/cm^2/A"
            # or perhaps 'relative' since these are fibre spectra?
            obscore["o_calib_status"] = "absolute"
            obscore["instrument_name"] = "2dF-AAOmega"
            obscore["em_calib_status"] = "calibrated"

            if "GAMANAME" in hdr:
                obscore["target_name"] = hdr["GAMANAME"]
            # if('OBJECT' in hdr):
            #    obscore['alt_target_name'] = hdr['OBJECT']
            # if('WAVRESOL' in hdr):
            #    obscore['em_resolution'] = hdr['WAVRESOL']*1e-10

            if "SN" in hdr:
                obscore["em_snr"] = hdr["SN"]

            # alternative name: OBJECT
            if "SPECID" in hdr:
                obscore["obs_id"] = hdr["SPECID"]
            if "Z" in hdr:
                obscore["redshift"] = hdr["Z"]

    return spectra


@data_loader(
    label="GAMA-2QZ obscore", extensions=FITS_FILE_EXTS, dtype=SpectrumList,
    identifier=no_auto_identify,
)
def gama_2qz_obscore_loader(fname):
    spectra = SpectrumList.read(fname, format="GAMA-2QZ")
    for spec in spectra:
        spec.meta["obscore"] = {}
        obscore = spec.meta["obscore"]
        if spec.meta["purpose"] == "reduced":

            hdr = spec.meta["header"]
            # NOTE: It's unclear whether DATE-OBS in header (Image generation
            # time) corresponds to date of observation or reduction Missing
            # start and end times of exposures; Missing exposure time
            t1, t2 = get_times(hdr, date_kw="DATE-OBS")
            if t1 is not None:
                obscore["t_min"] = t1.to_value('mjd',subfmt='float')
            if t2 is not None:
                obscore["t_max"] = t2.to_value('mjd',subfmt='float')

            obscore["s_ra"] = hdr["RA"]
            obscore["s_dec"] = hdr["DEC"]
            # obscore["obs_collection"] = "gama_dr2"
            # obscore["facility_name"] = "2qz"
            obscore["dataproduct_subtype"] = "science"
            obscore["calib_level"] = 2

            # if('EXPTIME' in hdr):
            # Exptime missing from header, assuming 1 hour (Boyle et al. 2000)
            obscore["t_exptime"] = 3600

            nspecpix = len(spec.spectral_axis)
            obscore["em_xel"] = nspecpix
            obscore["em_ucd"] = "em.wl"
            obscore["em_unit"] = "angstrom"
            obscore["s_xel1"] = nspecpix
            obscore["s_xel2"] = 1
            obscore["t_xel"] = 1
            obscore["em_min"] = spec.spectral_axis[0].meter
            obscore["em_max"] = (
                spec.spectral_axis[nspecpix - 1].meter
            )

            obscore["o_ucd"] = "phot.count"
            obscore["instrument_name"] = "2dF"
            obscore["s_fov"] = 2.1 / 3600
            obscore["em_calib_status"] = "calibrated"

            if "GAMANAME" in hdr:
                obscore["target_name"] = hdr["GAMANAME"]
            if "OBJECT" in hdr:
                obscore["alt_target_name"] = hdr["OBJECT"]
            if "WAVRESOL" in hdr:
                obscore["em_resolution"] = 9 * 1e-10

            # if('SN' in hdr):
            #    obscore['em_snr'] = hdr['SN']

            # alternative name: OBJECT
            if "SPECID" in hdr:
                obscore["obs_id"] = hdr["SPECID"]
            if "Z" in hdr:
                obscore["redshift"] = hdr["Z"]

    return spectra


def gama_lt_get_times(hdr):
    # format (start time): yyyy-mm-ddThh:mm:ss
    date = hdr.get("DATE-OBS")
    if not date:
        return None, None
    t1 = Time(date, format="isot", scale="utc")
    exptime = hdr.get("EXPTIME")
    if not exptime:
        return t1, None
    t2 = t1 + TimeDelta(exptime, format="sec")
    return t1, t2


@data_loader(
    label="GAMA-LT obscore", extensions=FITS_FILE_EXTS, dtype=SpectrumList,
    identifier=no_auto_identify,
)
def gama_lt_obscore_loader(fname):
    spectra = SpectrumList.read(fname, format="GAMA-LT")
    for spec in spectra:
        spec.meta["obscore"] = {}
        obscore = spec.meta["obscore"]
        if spec.meta["purpose"] == "reduced":

            hdr = spec.meta["header"]
            t1, t2 = gama_lt_get_times(hdr)
            if t1 is not None:
                obscore["t_min"] = t1.to_value('mjd',subfmt='float')
            if t2 is not None:
                obscore["t_max"] = t2.to_value('mjd',subfmt='float')

            obscore["s_ra"] = hdr["RA"]
            obscore["s_dec"] = hdr["DEC"]
            if "ESTSEE" in hdr:
                obscore["s_seeing"] = hdr["ESTSEE"]
            # obscore["obs_collection"] = "gama_dr2"
            # obscore["facility_name"] = "gama_lt"
            obscore["dataproduct_subtype"] = "science"
            obscore["calib_level"] = 2

            if "EXPTIME" in hdr:
                obscore["t_exptime"] = hdr["EXPTIME"]

            nspecpix = len(spec.spectral_axis)
            obscore["em_xel"] = nspecpix
            obscore["em_ucd"] = "em.wl"
            obscore["em_unit"] = "angstrom"
            obscore["s_xel1"] = nspecpix
            obscore["s_xel2"] = 1
            obscore["t_xel"] = 1
            obscore["em_min"] = spec.spectral_axis[0].meter
            obscore["em_max"] = (
                spec.spectral_axis[nspecpix - 1].meter
            )

            obscore["o_ucd"] = "phot.count"
            obscore["instrument_name"] = "FrodoSpec"
            obscore["s_fov"] = 9.84 / 3600
            obscore["em_calib_status"] = "calibrated"

            if "GAMANAME" in hdr:
                obscore["target_name"] = hdr["GAMANAME"]
            if "OBJECT" in hdr:
                obscore["alt_target_name"] = hdr["OBJECT"]
            if "WAVRESOL" in hdr:
                obscore["em_resolution"] = hdr["WAVRESOL"] * 1e-10

            # if('SN' in hdr):
            #    obscore['em_snr'] = hdr['SN']

            # alternative name: OBJECT
            if "SPECID" in hdr:
                obscore["obs_id"] = hdr["SPECID"]
            if "Z" in hdr:
                obscore["redshift"] = hdr["Z"]

    return spectra


@data_loader(
    label="GAMA-MGC obscore", extensions=FITS_FILE_EXTS, dtype=SpectrumList,
    identifier=no_auto_identify,
)
def gama_mgc_obscore_loader(fname):
    spectra = SpectrumList.read(fname, format="GAMA-MGC")
    for spec in spectra:
        spec.meta["obscore"] = {}
        obscore = spec.meta["obscore"]
        if spec.meta["purpose"] == "reduced":

            hdr = spec.meta["header"]
            t1, t2 = get_times(
                hdr, date_kw="DATE-OBS", start_kw="UT", duration_kw="EXPTIME"
            )
            if t1 is not None:
                obscore["t_min"] = t1.to_value('mjd',subfmt='float')
            if t2 is not None:
                obscore["t_max"] = t2.to_value('mjd',subfmt='float')

            obscore["s_ra"] = hdr["RA"]
            obscore["s_dec"] = hdr["DEC"]
            obscore["s_fov"] = 2.1 / 3600
            if "SEEING" in hdr:
                obscore["s_seeing"] = hdr["SEEING"]
            # obscore["obs_collection"] = "gama_dr2"
            # obscore["facility_name"] = "mgc"
            obscore["dataproduct_subtype"] = "science"
            obscore["calib_level"] = 2

            if "EXPTIME" in hdr:
                obscore["t_exptime"] = hdr["EXPTIME"]

            nspecpix = len(spec.spectral_axis)
            obscore["em_xel"] = nspecpix
            obscore["em_ucd"] = "em.wl"
            obscore["em_unit"] = "angstrom"
            obscore["s_xel1"] = nspecpix
            obscore["s_xel2"] = 1
            obscore["t_xel"] = 1
            obscore["em_min"] = spec.spectral_axis[0].meter
            obscore["em_max"] = (
                spec.spectral_axis[nspecpix - 1].meter
            )

            obscore["o_ucd"] = "phot.count"
            obscore["instrument_name"] = "2dF"
            obscore["em_calib_status"] = "calibrated"

            if "GAMANAME" in hdr:
                obscore["target_name"] = hdr["GAMANAME"]
            if "OBJECT" in hdr:
                obscore["alt_target_name"] = hdr["OBJECT"]
            if "SN" in hdr:
                obscore["em_snr"] = hdr["SN"]

            # alternative name: OBJECT
            if "SPECID" in hdr:
                obscore["obs_id"] = hdr["SPECID"]
            if "Z" in hdr:
                obscore["redshift"] = hdr["Z"]

    return spectra


@data_loader(
    label="GAMA-2SLAQ-LRG obscore", extensions=FITS_FILE_EXTS,
    dtype=SpectrumList, identifier=no_auto_identify,
)
def gama_2slaq_lrg_obscore_loader(fname):
    spectra = SpectrumList.read(fname, format="2SLAQ-LRG")

    # restrict top just first item in SpectrumList, which is the science
    # spectrum second in list is the sky spectrum
    spec = spectra[0]
    spec.meta["obscore"] = {}
    obscore = spec.meta["obscore"]
    hdr = spec.meta["header"]
    t1, t2 = get_times(hdr, duration_kw="EXPOSED")
    if t1 is not None:
        obscore["t_min"] = t1.to_value('mjd',subfmt='float')
    if t2 is not None:
        obscore["t_max"] = t2.to_value('mjd',subfmt='float')

    obscore["s_ra"] = hdr["RA"]
    obscore["s_dec"] = hdr["DEC"]
    obscore["s_fov"] = 2.1 / 3600
    # obscore["obs_collection"] = "gama_dr2"
    # obscore["facility_name"] = "2slaq-lrg"
    obscore["dataproduct_subtype"] = "science"
    obscore["calib_level"] = 2

    if "EXPOSED" in hdr:
        obscore["t_exptime"] = hdr["EXPOSED"]

    nspecpix = len(spec.spectral_axis)
    obscore["em_xel"] = nspecpix
    obscore["em_ucd"] = "em.wl"
    obscore["em_unit"] = "angstrom"
    obscore["s_xel1"] = nspecpix
    obscore["s_xel2"] = 1
    obscore["t_xel"] = 1
    obscore["em_min"] = spec.spectral_axis[0].meter
    obscore["em_max"] = (
        spec.spectral_axis[nspecpix - 1].meter
    )

    obscore["o_ucd"] = "phot.count"
    obscore["instrument_name"] = "2dF"
    obscore["em_calib_status"] = "calibrated"

    if "GAMANAME" in hdr:
        obscore["target_name"] = hdr["GAMANAME"]
    if "OBJECT" in hdr:
        obscore["alt_target_name"] = hdr["OBJECT"]
    if "SPECID" in hdr:
        obscore["obs_id"] = hdr["SPECID"]
    if "Z" in hdr:
        obscore["redshift"] = hdr["Z"]
    return spectra


@data_loader(
    label="GAMA-2SLAQ-QSO obscore", extensions=FITS_FILE_EXTS,
    dtype=SpectrumList, identifier=no_auto_identify,
)
def gama_2slaq_qso_obscore_loader(fname):
    spectra = SpectrumList.read(fname, format="GAMA-2SLAQ-QSO")
    for spec in spectra:
        spec.meta["obscore"] = {}
        obscore = spec.meta["obscore"]

        hdr = spec.meta["header"]
        t1, t2 = get_times(hdr, duration_kw="EXPOSED")
        if t1 is not None:
            obscore["t_min"] = t1.to_value('mjd',subfmt='float')
        if t2 is not None:
            obscore["t_max"] = t2.to_value('mjd',subfmt='float')

        obscore["s_ra"] = hdr["RA"]
        obscore["s_dec"] = hdr["DEC"]
        obscore["s_fov"] = 2.1 / 3600

        # obscore["obs_collection"] = "gama_dr2"
        # obscore["facility_name"] = "2slaq-qso"
        obscore["dataproduct_subtype"] = "science"
        obscore["calib_level"] = 2

        if "EXPOSED" in hdr:
            obscore["t_exptime"] = hdr["EXPOSED"]

        nspecpix = len(spec.spectral_axis)
        obscore["em_xel"] = nspecpix
        obscore["em_ucd"] = "em.wl"
        obscore["em_unit"] = "angstrom"
        obscore["s_xel1"] = nspecpix
        obscore["s_xel2"] = 1
        obscore["t_xel"] = 1
        obscore["em_min"] = spec.spectral_axis[0].meter
        obscore["em_max"] = (
            spec.spectral_axis[nspecpix - 1].meter
        )

        obscore["o_ucd"] = "phot.count"
        obscore["instrument_name"] = "2dF"
        obscore["em_calib_status"] = "calibrated"

        if "GAMANAME" in hdr:
            obscore["target_name"] = hdr["GAMANAME"]
        if "SPECID" in hdr:
            obscore["obs_id"] = hdr["SPECID"]
        if "Z" in hdr:
            obscore["redshift"] = hdr["Z"]

    return spectra
