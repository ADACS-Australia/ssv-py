import re

from astropy.time import Time, TimeDelta


def get_times(
    hdr, *, date_kw="UTDATE", start_kw="UTSTART", end_kw="UTEND",
    duration_kw="REFEXP"
):
    obs_date = hdr.get(date_kw)
    if not obs_date:
        return None, None

    date = re.sub(":", "-", obs_date)
    if start_kw in hdr:
        tstart = hdr[start_kw]
        tstr = date + "T" + tstart
    # we are missing UTSTART
    else:
        tstr = date
    t1 = Time(tstr, format="isot", scale="utc")

    if end_kw in hdr:
        tend = hdr[end_kw]
        tstr = date + "T" + tend
        t2 = Time(tstr, format="isot", scale="utc")
    elif duration_kw in hdr:
        exptime = hdr[duration_kw]
        t2 = t1 + TimeDelta(exptime, format="sec")
    else:
        return t1, None
    return t1, t2
