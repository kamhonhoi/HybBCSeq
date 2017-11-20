#!/usr/bin/python

# Kam Hon Hoi
# 11/20/17

# This script is to copy over the optimized_pwm package that was compiled without the 
# fpectl option

# When the ''' undefined symbol: PyFPE_jbuf ''' occured, it means an incompatible package was installed.

# Run this script to overwrite the exisiting optimized_pwm package with the one not compiled with the fpectl option

import os

os.system("cp optimized_pwm.so ../../HybBCSeq-venv/lib/python2.7/site-packages/")

