# flake8: noqa
uselocal = True
if uselocal:
    from . import loaders
    from . import aaomega_2df
    from . import gama
    from . import galah
    from . import marz
    from . import ozdes
    from . import obscore
else:
    import dc_loaders as loaders
from .plotting import plot_spectra
from .helpers import (
    specutils_spectra_to_table_spectra, specutils_spectrum_to_table_spectrum,
)
from ._version import get_versions
__version__ = get_versions()['version']
del get_versions
