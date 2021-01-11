#!/bin/bash
python fits2json.py | node ../../src/plugins/marzcli.js --verbose  --numCPUs=0  /dev/stdin 
