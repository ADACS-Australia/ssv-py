import ssv.plugin_collection as plugin_collection
import ssv.utils as utils
from ssv.viewer import SimpleSpectrum
from Naked.toolshed.shell import execute_js, muterun_js
import sys
import os
import json
from pathlib import Path
from astropy.io import fits

class MarzCLI(plugin_collection.Plugin):
    """Marz CLI plugin; can be used to find the best fit of redshift and template to the input spectrum
    """
    def __init__(self):
        super().__init__()
        self.description = 'MarzCLI'

    def reduce(self, argument):
        """Implementation of the Marz CLI

        Parameters
        ----------
        argument : str
            Path to the FITS file containing the spectrum to be fit

        Returns
        -------
        Float
            The fitted redshift of the input spectrum
        String
            The name of the fitted template of the input spectrum
        """

        import os
        from subprocess import run, Popen, PIPE, STDOUT

        spectrum = SimpleSpectrum('fits2JSON', argument)
        spectrum_json = json.dumps(utils.toMarzJSON(spectrum), indent=2)

        plugin_dir = os.path.dirname(os.path.abspath(__file__))

        marzcli_path = plugin_dir + '/marzcli.js'
        cmd = plugin_dir + '/marzcli.sh ' + marzcli_path

        proc = run(cmd, shell=True, input=spectrum_json, text=True, stdout=PIPE, stderr=STDOUT, close_fds=True)
        lines = proc.stdout.split('\n')


        print("="*20 + " Marz CLI Plugin " + "="*20)
        for line in lines[:-6:-1]:
            print(line)
        
        tokens = lines[-5].split(',')
        bestfit_redshift = float(tokens[8])
        bestfit_template = tokens[7].strip()
        print(f"BEST REDSHIFT = {bestfit_redshift}")
        print(f"BEST TEMPLATE = {bestfit_template}")
        return bestfit_redshift, bestfit_template

