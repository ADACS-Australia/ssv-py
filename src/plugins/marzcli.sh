#!/bin/bash
python fits2json.py --input $2 | node $1 --verbose  --numCPUs=0  /dev/stdin 
