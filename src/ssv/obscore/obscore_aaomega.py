import re

import numpy as np

from specutils import SpectrumList
from specutils.io.registers import data_loader

from .common import get_times
from ..loaders import FITS_FILE_EXTS, no_auto_identify


def calc_aaomega_resolutions(hdr):
    # TODO: Add MODE = MOS / IFU
    grat = hdr.get("GRATID")
    gang = hdr.get("GRATANGL")
    cang = hdr.get("CAMANGL")
    order = hdr.get("ORDER")
    if not (grat and gang and cang and order):
        return None, None, None, None
    # perhaps I should get this from somewhere in header?
    npixels = 2048
    # use value for MOS
    resolutionpix = 3.4
    # check if hdr['INSTRUME'] contains KOALA or SPIRAL ???
    # resolutionpix = 2.1

    rad = 180.0 / np.pi
    flcam = 247
    pix = 0.015
    hwid = (npixels * pix / flcam) / 2.0
    ddisp = np.cos(72.5 / (3.15 * 190))

    slant = {
        "580V": 0.7,
        "385R": 0.6,
        "1700B": 0.2,
        "1500V": 0.0,
        "1000R": 1.2,
        "1000I": 1.8,
        "3200R": 0.0,
        "2500V": 0.0,
        "2000R": 0.0,
        "1700I": 0.7,
        "1700D": 0.2,
    }

    slantr = slant[grat] / rad
    linespmm = int(grat[:-1])

    gangr = gang / rad
    cangr = cang / rad - gangr
    lcen = 1e7 * (np.sin(gangr) + np.sin(cangr)) / (linespmm * order)
    lblaze = 1e7 * 2 * np.sin(gangr + slantr) / (linespmm * order)
    # Get central and blaze wavelengths
    lcen = int(lcen + 0.5)
    lblaze = int(lblaze + 0.5)

    dispc = 1e7 * pix / flcam * np.cos(cangr) / (order * linespmm)
    resolutionpix = resolutionpix * np.cos(gangr) / np.cos(cangr)

    resa = resolutionpix * dispc
    res = lcen / resa
    lcb = 1e7 * (np.sin(gangr) + np.sin(cangr - hwid)) / (order * linespmm)
    lcr = 1e7 * (np.sin(gangr) + np.sin(cangr + hwid)) / (order * linespmm)
    leb = ddisp * lcb
    ler = ddisp * lcr

    dcb = 1e7 * pix / flcam * np.cos(cangr - hwid) / (order * linespmm)
    dcr = 1e7 * pix / flcam * np.cos(cangr + hwid) / (order * linespmm)
    deb = ddisp * dcb
    der = ddisp * dcr

    racb = resolutionpix * dcb
    racr = resolutionpix * dcr
    raeb = racb * ddisp
    raer = racr * ddisp
    rcb = lcb / (resolutionpix * dcb)
    rcr = lcr / (resolutionpix * dcr)
    reb = rcb / ddisp
    rer = rcr / ddisp

    dispc = int((1000 * dispc) + 0.5) / 1000
    resa = int((1000 * resa) + 0.5) / 1000
    res = int(res + 0.5)
    resa = int((1000 * resa) + 0.5) / 1000
    res = int(res + 0.5)
    lcb = int(lcb + 0.5)
    lcr = int(lcr + 0.5)
    leb = int(leb + 0.5)
    ler = int(ler + 0.5)
    dcb = int((1000 * dcb) + 0.5) / 1000
    dcr = int((1000 * dcr) + 0.5) / 1000
    deb = int((1000 * deb) + 0.5) / 1000
    der = int((1000 * der) + 0.5) / 1000
    racb = int((1000 * racb) + 0.5) / 1000
    racr = int((1000 * racr) + 0.5) / 1000
    raeb = int((1000 * raeb) + 0.5) / 1000
    raer = int((1000 * raer) + 0.5) / 1000
    rcb = int(rcb + 0.5)
    rcr = int(rcr + 0.5)
    reb = int(reb + 0.5)
    rer = int(rer + 0.5)
    covc = lcr - lcb
    cove = ler - leb
    cov = ler - lcb

    cen_res = resa
    cen_rp = res
    cen_rp_min = rcb
    cen_rp_max = rcr
    # cen_res = FWHM in Angstrom;
    # cen_rp = resolving power at central wavelength
    # cen_rp_min = min resolving power
    # cen_rp_max = max resolving power
    return cen_res, cen_rp, cen_rp_min, cen_rp_max


@data_loader(
    label="Data Central AAOmega obscore", extensions=FITS_FILE_EXTS,
    dtype=SpectrumList, identifier=no_auto_identify,
)
def aaomega_obscore_loader(fname):
    spectra = SpectrumList.read(fname, format="Data Central AAOmega")

    for spec in spectra:
        cen_res, cen_rp, cen_rp_min, cen_rp_max = calc_aaomega_resolutions(
            spec.meta["header"]
        )
        if cen_res is not None:
            break

    for spec in spectra:
        # Don't produce obscore values for sky
        if spec.meta["purpose"] == "reduced":
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
            obscore["s_seeing"] = hdr.get("SEEING")
            obscore["obs_collection"] = "aat_archive"
            obscore["facility_name"] = "AAT"
            # obscore["dataproduct_type"] = "spectrum"
            obscore["dataproduct_subtype"] = "science"
            obscore["calib_level"] = 2

            obscore["t_exptime"] = hdr.get("EXPOSED")

            nspecpix = len(spec.spectral_axis)
            obscore["em_xel"] = nspecpix
            obscore["em_ucd"] = "em.wl"
            obscore["em_unit"] = "angstrom"
            obscore["s_xel1"] = nspecpix
            obscore["s_xel2"] = 1
            obscore["t_xel"] = 1
            obscore["em_min"] = spec.spectral_axis[0].meter
            obscore["em_max"] = spec.spectral_axis[-1].meter
            obscore["em_res_power"] = cen_rp
            obscore["em_res_power_min"] = cen_rp_min
            obscore["em_res_power_max"] = cen_rp_max
            obscore["em_resolution"] = cen_res * 1e-10

            obscore["o_ucd"] = "phot.count"
            # spectra are calibrated in flux: ROW1 hdr comment: Flux-calibrated
            # spectrum in 10^-17 erg/s/cm^2/A
            # obscore['o_ucd'] = 'phot.flux'
            # obscore['o_unit'] = '1.0E-17 erg/s/cm^2/A'
            # obscore['o_calib_status'] = 'absolute'
            # or perhaps 'relative' since these are fibre spectra?

            obscore["instrument_name"] = "2dF-AAOmega"
            obscore["em_calib_status"] = "calibrated"

            if "OBJECT" in hdr:
                obscore["target_name"] = hdr["OBJECT"]
            if "OBJCOM" in hdr:
                obscore["alt_target_name"] = hdr["OBJCOM"]

            # alternative name: OBJECT
            if "OBJPIV" in hdr:
                obscore["obs_id"] = "%s-%s" % (
                    re.sub(".sds", "", hdr["CFG_FILE"]),
                    hdr["OBJPIV"],
                )

    return spectra
