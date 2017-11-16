#!/usr/bin/python

# Kam Hon Hoi
# 11/16/2017

# This script is used to setup the virtual environment for use with the Hybridoma Barcoded Sequencing Workflow

import os
import virtualenv

# Creating the virtualenv
venv_dir="HybBCSeq-venv"
virtualenv.create_environment(venv_dir)
# Installing custom Python package
os.system("cp ./HybBCSeq-bin/*.so ./HybBCSeq-venv/lib/python2.7/site-packages/")
