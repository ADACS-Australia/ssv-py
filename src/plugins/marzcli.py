import ssv.plugin_collection as plugin_collection
from Naked.toolshed.shell import execute_js, muterun_js
import sys
import os
from pathlib import Path

class MarzCLI(plugin_collection.Plugin):
    """This plugin is just the identity function: it returns the argument
    """
    def __init__(self):
        super().__init__()
        self.description = 'MarzCLI'

    def perform_operation(self, argument):
        """The actual implementation of the identity plugin is to just return the
        argument
        """
        return argument
    
    def reduce_old(self, argument):
        """The actual implementation of the identity plugin is to just return the
        argument
        """
        import os
        climarzv2 = str(os.path.dirname(os.path.abspath(__file__)))+'/marzcli.js'
        infile = argument
        filename, file_extension = os.path.splitext(infile)
        result = muterun_js(climarzv2, arguments="--numCPUs=0 "+argument)
        Lines = result.stdout.decode('utf-8').split('\n')
        
        header = Lines[-3][1:].split(',')
        tokens = Lines[-3].split(',')
        bestfit_redshift = float(tokens[8])
        bestfit_template = tokens[7].strip()
        return (bestfit_redshift, bestfit_template)

    def reduce(self, argument):
        """The actual implementation of the identity plugin is to just return the
        argument
        """
        import os
        from subprocess import Popen, PIPE, STDOUT

        cmd = str(os.path.dirname(os.path.abspath(__file__)))+'/marzcli.sh '+argument
        p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
        output = p.stdout.read()
        Lines = output.decode('utf-8').split('\n')
        print("last line"+str(Lines[-1]))
        print("second last line"+str(Lines[-2]))
        print("third last line"+str(Lines[-3]))
        print("forth last line"+str(Lines[-4]))
        print("fifth last line"+str(Lines[-5]))
        header = Lines[-6][1:].split(',')
        tokens = Lines[-5].split(',')
        bestfit_redshift = float(tokens[8])
        bestfit_template = tokens[7].strip()
        return (bestfit_redshift, bestfit_template)

